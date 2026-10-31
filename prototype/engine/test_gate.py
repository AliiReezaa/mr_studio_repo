"""Gate test for the engine connection (contracts C10). Acts as the desktop shell:
starts the engine, passes secrets on stdin, then attacks and uses the connection.
Usage: python3 test_gate.py"""
import asyncio, http.client as httpc, json, os, secrets, subprocess, sys, time
import websockets

HERE = os.path.dirname(os.path.abspath(__file__))
results = []

def check(desc, ok):
    results.append(ok)
    print(("PASS " if ok else "FAIL ") + desc)

def http(port, path, headers=None):
    c = httpc.HTTPConnection("127.0.0.1", port, timeout=10)
    c.request("GET", path, headers=headers or {})
    r = c.getresponse(); body = r.read(); c.close()
    return r.status, dict(r.getheaders()), body

async def ws_try(port, origin, cookie, host=None):
    """Returns 'open' if the server keeps the socket open, 'closed' if it closes it."""
    hdrs = {}
    if cookie: hdrs["Cookie"] = cookie
    if host: hdrs["Host"] = host
    try:
        async with websockets.connect(f"ws://127.0.0.1:{port}/websocket/", origin=origin,
                                      additional_headers=hdrs, open_timeout=5) as ws:
            try:
                await asyncio.wait_for(ws.recv(), timeout=1.5)
                return "open"            # server sent something: socket accepted
            except asyncio.TimeoutError:
                return "open"            # still open, waiting for client init
            except websockets.ConnectionClosed:
                return "closed"
    except Exception as e:
        return "closed"

async def ws_session(port, origin, cookie):
    async with websockets.connect(f"ws://127.0.0.1:{port}/websocket/", origin=origin,
                                  additional_headers={"Cookie": cookie}, open_timeout=5) as ws:
        await ws.send(json.dumps({"method": "init", "data": {".clientdata_output_hello_hidden": False}}))
        deadline = time.time() + 10
        while time.time() < deadline:
            msg = await asyncio.wait_for(ws.recv(), timeout=10)
            if "connected" in msg:
                return True
    return False

def main():
    secret, token = secrets.token_hex(32), secrets.token_hex(32)
    # MRGUI_R_CMD: the engine command the shell uses, e.g. '<runtime>/lib/R/bin/exec/R --no-echo --no-restore'
    rcmd = os.environ.get("MRGUI_R_CMD", "Rscript").split("|")
    script = os.path.join(HERE, "main.R")
    argv = rcmd + ([script] if rcmd == ["Rscript"] else ["--file=" + script])
    p = subprocess.Popen(argv, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    p.stdin.write(json.dumps({"secret": secret, "login_token": token}) + "\n"); p.stdin.flush(); p.stdin.close()
    t0 = time.time(); line = p.stdout.readline()
    try:
        ready = json.loads(line)
    except Exception:
        print("engine did not start:", line, p.stderr.read()[-2000:]); sys.exit(1)
    port = ready["port"]
    check(f"engine reports ready on stdout ({time.time()-t0:.1f}s, port {port})", ready.get("ready") is True)
    if sys.platform.startswith("linux"):
        check("secret not visible in the engine command line", secret not in open(f"/proc/{p.pid}/cmdline").read())
    host = f"127.0.0.1:{port}"
    try:
        s, h, _ = http(port, "/")
        check("page without cookie -> 403", s == 403)
        s, h, _ = http(port, "/auth?t=" + "0" * 64)
        check("wrong login token -> 403", s == 403)
        s, h, _ = http(port, "/auth?t=" + token, {"Host": "evil.example:80"})
        check("login with foreign Host header -> 403 (DNS rebinding)", s == 403)
        s, h, _ = http(port, "/auth?t=" + token)
        sc = h.get("Set-Cookie", "")
        check("correct one-time token -> 302 with HttpOnly SameSite=Strict cookie",
              s == 302 and "HttpOnly" in sc and "SameSite=Strict" in sc and h.get("Location") == "/")
        cookie = sc.split(";")[0]
        s, h, _ = http(port, "/auth?t=" + token)
        check("login token cannot be reused -> 403", s == 403)
        s, h, b = http(port, "/", {"Cookie": cookie})
        check("page with cookie -> 200 and Shiny UI", s == 200 and b"MR GUI engine prototype" in b)
        s, h, _ = http(port, "/", {"Cookie": cookie, "Host": "evil.example"})
        check("valid cookie but foreign Host -> 403", s == 403)
        s, h, _ = http(port, "/", {"Cookie": "mrgui_s=" + "f" * 64})
        check("forged cookie -> 403", s == 403)
        s, h, _ = http(port, "/session/abc/download/x", {})
        check("session download route without cookie -> 403", s == 403)
        s, h, _ = http(port, "/shared/shiny.min.js")
        print(f"INFO static library file /shared/shiny.min.js without cookie -> {s} (served by httpuv outside R; public file only)")

        good = f"http://127.0.0.1:{port}"
        check("WebSocket with cookie and own Origin stays open", asyncio.run(ws_try(port, good, cookie)) == "open")
        check("WebSocket from a foreign web page Origin is closed", asyncio.run(ws_try(port, "http://evil.example", cookie)) == "closed")
        check("full Shiny session over the gated WebSocket renders server output", asyncio.run(ws_session(port, good, cookie)))
        check("WebSocket without cookie is closed", asyncio.run(ws_try(port, good, None)) == "closed")
        check("WebSocket with foreign Host is closed", asyncio.run(ws_try(port, good, cookie, host="evil.example")) == "closed")
    finally:
        p.terminate()
        try: p.wait(timeout=10)
        except subprocess.TimeoutExpired: p.kill()
    print(f"\n{sum(results)} passed, {len(results) - sum(results)} failed")
    sys.exit(0 if all(results) else 1)

main()
