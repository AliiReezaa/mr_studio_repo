//! Smoke test run by the build pipeline on each OS (MRGUI_SMOKE=1): checks the
//! bundled engine through the real C10 handshake, then exits with 0 or 1.

use std::io::{Read, Write};
use std::net::TcpStream;
use std::time::Duration;

struct Resp { status: u16, set_cookie: Option<String>, body: String }

fn get(port: u16, path: &str, cookie: Option<&str>) -> Result<Resp, String> {
    let mut s = TcpStream::connect(("127.0.0.1", port)).map_err(|e| e.to_string())?;
    s.set_read_timeout(Some(Duration::from_secs(20))).ok();
    let mut req = format!("GET {path} HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\nConnection: close\r\n");
    if let Some(c) = cookie { req.push_str(&format!("Cookie: {c}\r\n")); }
    req.push_str("\r\n");
    s.write_all(req.as_bytes()).map_err(|e| e.to_string())?;
    let mut raw = Vec::new();
    s.read_to_end(&mut raw).map_err(|e| e.to_string())?;
    let text = String::from_utf8_lossy(&raw).to_string();
    let (head, body) = text.split_once("\r\n\r\n").unwrap_or((&text, ""));
    let status = head.split_whitespace().nth(1).and_then(|c| c.parse().ok()).unwrap_or(0);
    let set_cookie = head.lines()
        .find(|l| l.to_ascii_lowercase().starts_with("set-cookie:"))
        .map(|l| l[11..].trim().split(';').next().unwrap_or("").to_string());
    Ok(Resp { status, set_cookie, body: body.to_string() })
}

pub fn run(port: u16, token: &str) -> bool {
    let mut ok = true;
    let mut check = |name: &str, pass: bool| { println!("{} {name}", if pass { "PASS" } else { "FAIL" }); ok &= pass; };
    match get(port, "/", None) { Ok(r) => check("page without cookie is refused (403)", r.status == 403), Err(e) => check(&e, false) }
    match get(port, &format!("/auth?t={token}"), None) {
        Ok(r) => {
            check("one-time login gives 302 and a cookie", r.status == 302 && r.set_cookie.is_some());
            let cookie = r.set_cookie.unwrap_or_default();
            match get(port, "/", Some(&cookie)) {
                Ok(p) => check("engine page served with cookie (200)", p.status == 200 && p.body.contains("MR GUI engine prototype")),
                Err(e) => check(&e, false),
            }
            match get(port, &format!("/auth?t={token}"), None) {
                Ok(p) => check("login token cannot be reused (403)", p.status == 403),
                Err(e) => check(&e, false),
            }
        }
        Err(e) => check(&e, false),
    }
    ok
}
