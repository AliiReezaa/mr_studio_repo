//! Starts and stops the bundled R engine (contracts C10).
//!
//! The secret and one-time login token are generated here and passed on the
//! engine's standard input, never as arguments or environment variables. The
//! engine answers with one JSON line on standard output once it is listening.

use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::mpsc;
use std::time::Duration;

pub struct Engine {
    pub child: Child,
    pub port: u16,
    pub login_token: String,
}

fn random_hex(bytes: usize) -> Result<String, String> {
    let mut buf = vec![0u8; bytes];
    getrandom::getrandom(&mut buf).map_err(|e| format!("random generator failed: {e}"))?;
    Ok(buf.iter().map(|b| format!("{b:02x}")).collect())
}

/// Locations of the R home and engine scripts. Development builds may point
/// elsewhere with MRGUI_R_HOME and MRGUI_ENGINE_DIR.
pub fn locate(resource_dir: &Path) -> (PathBuf, PathBuf) {
    let r_home = std::env::var_os("MRGUI_R_HOME").map(PathBuf::from).unwrap_or_else(|| {
        if cfg!(windows) { resource_dir.join("runtime").join("Lib").join("R") }
        else { resource_dir.join("runtime").join("lib").join("R") }
    });
    let engine_dir = std::env::var_os("MRGUI_ENGINE_DIR").map(PathBuf::from)
        .unwrap_or_else(|| resource_dir.join("engine"));
    (r_home, engine_dir)
}

/// The R executable started directly (the `Rscript` wrapper is not relocatable).
fn r_binary(r_home: &Path) -> PathBuf {
    if cfg!(windows) { r_home.join("bin").join("x64").join("Rterm.exe") }
    else { r_home.join("bin").join("exec").join("R") }
}

pub fn start(r_home: &Path, engine_dir: &Path, log_path: &Path, timeout: Duration) -> Result<Engine, String> {
    let secret = random_hex(32)?;
    let login_token = random_hex(32)?;
    let bin = r_binary(r_home);
    if !bin.exists() {
        return Err(format!("R engine not found at {}", bin.display()));
    }
    if let Some(dir) = log_path.parent() { let _ = std::fs::create_dir_all(dir); }
    let log = std::fs::File::create(log_path).map_err(|e| format!("cannot create log: {e}"))?;

    let mut cmd = Command::new(&bin);
    cmd.arg("--no-echo").arg("--no-restore").arg("--vanilla")
        .arg(format!("--file={}", engine_dir.join("main.R").display()))
        .env("R_HOME", r_home)
        .env("R_SHARE_DIR", r_home.join("share"))
        .env("R_INCLUDE_DIR", r_home.join("include"))
        .env("R_DOC_DIR", r_home.join("doc"))
        .env("OMP_NUM_THREADS", "1")            // contracts C8.4: one thread per task
        .env("OPENBLAS_NUM_THREADS", "1")
        .env("MKL_NUM_THREADS", "1")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::from(log));
    if !cfg!(windows) {
        cmd.env("LANG", "en_US.UTF-8").env("LC_ALL", "en_US.UTF-8"); // contracts C10.1
    }
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        const CREATE_NO_WINDOW: u32 = 0x0800_0000;
        cmd.creation_flags(CREATE_NO_WINDOW);
    }

    let mut child = cmd.spawn().map_err(|e| format!("cannot start engine: {e}"))?;
    #[cfg(windows)]
    if let Err(e) = kill_with_shell(&child) {
        let _ = child.kill();
        return Err(e);
    }
    {
        let mut stdin = child.stdin.take().ok_or("engine stdin unavailable")?;
        let init = serde_json::json!({ "secret": secret, "login_token": login_token,
                                       "parent_pid": std::process::id() });
        writeln!(stdin, "{init}").map_err(|e| format!("cannot send secrets: {e}"))?;
    } // stdin dropped here: the engine sees end of input

    let stdout = child.stdout.take().ok_or("engine stdout unavailable")?;
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let mut line = String::new();
        let _ = BufReader::new(stdout).read_line(&mut line);
        let _ = tx.send(line);
    });
    let line = match rx.recv_timeout(timeout) {
        Ok(l) => l,
        Err(_) => { let _ = child.kill(); return Err("engine did not become ready in time".into()); }
    };
    let ready: serde_json::Value = serde_json::from_str(line.trim())
        .map_err(|_| { let _ = child.kill(); format!("unexpected engine output: {}", line.trim()) })?;
    let port = ready.get("port").and_then(|p| p.as_u64()).filter(|p| *p > 0 && *p < 65536)
        .ok_or_else(|| { let _ = child.kill(); "engine did not report a port".to_string() })? as u16;
    Ok(Engine { child, port, login_token })
}

/// Windows: put the engine in a job object that is killed when the last handle
/// to it closes. The handle is deliberately kept open for the life of the shell,
/// so the OS stops the engine if the shell exits or crashes.
#[cfg(windows)]
fn kill_with_shell(child: &Child) -> Result<(), String> {
    use std::os::windows::io::AsRawHandle;
    use windows_sys::Win32::System::JobObjects::{
        AssignProcessToJobObject, CreateJobObjectW, JobObjectExtendedLimitInformation,
        SetInformationJobObject, JOBOBJECT_EXTENDED_LIMIT_INFORMATION, JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE,
    };
    unsafe {
        let job = CreateJobObjectW(std::ptr::null(), std::ptr::null());
        if job.is_null() {
            return Err("cannot create job object".into());
        }
        let mut info: JOBOBJECT_EXTENDED_LIMIT_INFORMATION = std::mem::zeroed();
        info.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
        let ok = SetInformationJobObject(
            job,
            JobObjectExtendedLimitInformation,
            &info as *const _ as *const core::ffi::c_void,
            std::mem::size_of::<JOBOBJECT_EXTENDED_LIMIT_INFORMATION>() as u32,
        );
        if ok == 0 {
            return Err("cannot configure job object".into());
        }
        if AssignProcessToJobObject(job, child.as_raw_handle() as _) == 0 {
            return Err("cannot assign engine to job object".into());
        }
        // `job` is intentionally never closed: closing happens when the shell exits.
    }
    Ok(())
}
