// MR Studio desktop shell — week-1 prototype.
// Opens a local loading page, starts the R engine, then points the window at
// the engine's one-time login URL. Kills the engine when the app exits.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod engine;
mod smoke;

use std::sync::Mutex;
use std::time::Duration;
use tauri::{Manager, RunEvent, WebviewUrl, WebviewWindowBuilder};

struct EngineState(Mutex<Option<engine::Engine>>);

fn main() {
    let app = tauri::Builder::default()
        .manage(EngineState(Mutex::new(None)))
        .setup(|app| {
            let win = WebviewWindowBuilder::new(app, "main", WebviewUrl::App("index.html".into()))
                .title("MR Studio")
                .inner_size(1440.0, 900.0)
                .min_inner_size(1024.0, 700.0)
                .build()?;
            let loading_url = win.url()?;
            let resource_dir = app.path().resource_dir()?;
            let log_path = app.path().app_log_dir()?.join("engine.log");
            let handle = app.handle().clone();

            std::thread::spawn(move || {
                let (r_home, engine_dir) = engine::locate(&resource_dir);
                match engine::start(&r_home, &engine_dir, &log_path, Duration::from_secs(90)) {
                    Ok(eng) => {
                        if std::env::var_os("MRGUI_SMOKE").is_some() {
                            let passed = smoke::run(eng.port, &eng.login_token);
                            *handle.state::<EngineState>().0.lock().unwrap() = Some(eng);
                            handle.exit(if passed { 0 } else { 1 });
                            return;
                        }
                        let url = format!("http://127.0.0.1:{}/auth?t={}", eng.port, eng.login_token);
                        *handle.state::<EngineState>().0.lock().unwrap() = Some(eng);
                        if let Ok(u) = url.parse() { let _ = win.navigate(u); }
                    }
                    Err(msg) => {
                        if std::env::var_os("MRGUI_SMOKE").is_some() {
                            println!("FAIL engine start: {msg}");
                            handle.exit(1);
                            return;
                        }
                        let mut u = loading_url.clone();
                        u.query_pairs_mut().append_pair("error", &msg);
                        let _ = win.navigate(u);
                    }
                }
            });
            Ok(())
        })
        .build(tauri::generate_context!())
        .expect("error while building MR Studio");

    app.run(|handle, event| {
        if let RunEvent::Exit = event {
            if let Some(mut eng) = handle.state::<EngineState>().0.lock().unwrap().take() {
                let _ = eng.child.kill();
                let _ = eng.child.wait();
            }
        }
    });
}
