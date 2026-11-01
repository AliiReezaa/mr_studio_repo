#!/usr/bin/env python3
"""Make a conda-built macOS runtime relocatable.

conda writes the build prefix into some Mach-O load commands (for example
libncurses.6.dylib -> /<build prefix>/lib/libtinfo.6.dylib). After the runtime
folder is moved, or installed inside an app bundle, those absolute paths no
longer exist and the program aborts. This script rewrites every load command,
install name and rpath that starts with the build prefix to a path relative to
the file itself (@loader_path/...), then re-signs each changed file ad hoc
(required on Apple Silicon).

It also makes R's own launcher script (lib/R/bin/R) and etc/ldpaths relocatable. conda
writes the build prefix into them as plain text; `R CMD INSTALL`, `Rscript` and `R CMD`
start R through that launcher, so it must find its own folder wherever the runtime lives
(including a folder whose name contains a space, such as /Applications/MR Studio.app).

Usage: relocate_macos.py OLD_PREFIX NEW_PREFIX
  OLD_PREFIX  the folder the runtime was built in (as it appears in otool -L)
  NEW_PREFIX  the folder the runtime lives in now
Exit code 1 if any reference to OLD_PREFIX remains in a Mach-O file.
"""
import mmap
import os
import re
import subprocess
import sys

MAGICS = {b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe", b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca"}


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True)


def is_macho(path):
    try:
        with open(path, "rb") as fh:
            return fh.read(4) in MAGICS
    except OSError:
        return False


def mentions(path, needle):
    try:
        if os.path.getsize(path) == 0:
            return False
        with open(path, "rb") as fh, mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            return mm.find(needle) != -1
    except (OSError, ValueError):
        return False


def deps_and_id(path):
    res = run(["otool", "-L", path])
    if res.returncode != 0:
        return None, None
    lines = res.stdout.splitlines()[1:]
    deps = [ln.strip().split(" (")[0] for ln in lines if ln.strip()]
    res_id = run(["otool", "-D", path])
    ident = None
    if res_id.returncode == 0:
        out = res_id.stdout.strip().splitlines()
        if len(out) >= 2:
            ident = out[1].strip()
    return deps, ident


def rpaths(path):
    res = run(["otool", "-l", path])
    if res.returncode != 0:
        return []
    return re.findall(r"cmd LC_RPATH\n\s+cmdsize \d+\n\s+path (.+?) \(offset \d+\)", res.stdout)


LAUNCHER_HEAD = """# --- MR Studio: relocatable launcher. The install prefix is derived from this file's own location.
MRGUI_SELF="$0"
while [ -L "$MRGUI_SELF" ]; do
  MRGUI_LINK="$(readlink "$MRGUI_SELF")"
  case "$MRGUI_LINK" in /*) MRGUI_SELF="$MRGUI_LINK" ;; *) MRGUI_SELF="$(dirname "$MRGUI_SELF")/$MRGUI_LINK" ;; esac
done
MRGUI_PREFIX="$(cd "$(dirname "$MRGUI_SELF")/../../.." && pwd -P)"
"""


def relocate_text_files(old, new):
    """Rewrite the build prefix in R's launcher script and etc/ldpaths. Returns a list of problems."""
    problems = []
    r_home = os.path.join(new, "lib", "R")
    launcher = os.path.join(r_home, "bin", "R")
    if os.path.isfile(launcher):
        text = open(launcher, encoding="utf-8").read()
        if "MRGUI_PREFIX" not in text:
            lines = text.split("\n", 1)           # keep the #! line first
            body = lines[1].replace(old, "${MRGUI_PREFIX}")
            # the launcher uses ${SED} unquoted, which breaks when the folder name has a space
            # (for example /Applications/MR Studio.app); the system sed is always present
            body = body.replace("SED=${MRGUI_PREFIX}/bin/sed", "SED=sed")
            text = lines[0] + "\n" + LAUNCHER_HEAD + body
            open(launcher, "w", encoding="utf-8").write(text)
        if old in open(launcher, encoding="utf-8").read():
            problems.append(launcher)
    ldpaths = os.path.join(r_home, "etc", "ldpaths")
    if os.path.isfile(ldpaths):
        text = open(ldpaths, encoding="utf-8").read()
        text = text.replace(old, "${R_HOME}/../..")
        open(ldpaths, "w", encoding="utf-8").write(text)
        if old in text:
            problems.append(ldpaths)
    return problems


def main():
    old, new = sys.argv[1].rstrip("/"), os.path.realpath(sys.argv[2]).rstrip("/")
    needle = old.encode()
    changed = failed = 0
    for root, _dirs, files in os.walk(new, followlinks=False):
        for name in files:
            path = os.path.join(root, name)
            if os.path.islink(path) or not is_macho(path) or not mentions(path, needle):
                continue
            deps, ident = deps_and_id(path)
            if deps is None:
                continue
            here = os.path.dirname(path)
            cmds = []
            for dep in deps:
                if dep == ident:
                    continue
                if dep.startswith(old + "/"):
                    target = os.path.join(new, os.path.relpath(dep, old))
                    cmds.append(["install_name_tool", "-change", dep,
                                 "@loader_path/" + os.path.relpath(target, here), path])
            if ident and ident.startswith(old + "/"):
                cmds.append(["install_name_tool", "-id", "@rpath/" + os.path.basename(ident), path])
            for rp in rpaths(path):
                if rp.startswith(old):
                    target = os.path.join(new, os.path.relpath(rp, old))
                    rel = "@loader_path/" + os.path.relpath(target, here)
                    cmds.append(["install_name_tool", "-rpath", rp, rel, path])
            if not cmds:
                continue
            ok = True
            for cmd in cmds:
                res = run(cmd)
                if res.returncode != 0:
                    # a duplicate rpath cannot be rewritten; delete it instead
                    if cmd[1] == "-rpath":
                        res = run(["install_name_tool", "-delete_rpath", cmd[2], path])
                    if res.returncode != 0:
                        ok = False
                        print("WARN", " ".join(cmd), res.stderr.strip(), file=sys.stderr)
            sign = run(["codesign", "--force", "--sign", "-", path])
            if sign.returncode != 0:
                ok = False
                print("WARN codesign", path, sign.stderr.strip(), file=sys.stderr)
            changed += 1
            failed += 0 if ok else 1
    text_problems = relocate_text_files(old, new)
    for path in text_problems:
        print("  STILL (text):", path)
    # verification: no Mach-O load command may still point at the old prefix
    remaining = []
    for root, _dirs, files in os.walk(new, followlinks=False):
        for name in files:
            path = os.path.join(root, name)
            if os.path.islink(path) or not is_macho(path) or not mentions(path, needle):
                continue
            deps, ident = deps_and_id(path)
            bad = [d for d in (deps or []) if d.startswith(old)] + [r for r in rpaths(path) if r.startswith(old)]
            if bad:
                remaining.append((path, bad))
    print(f"relocate_macos: rewrote {changed} files, {failed} with warnings, {len(remaining)} still reference {old}")
    for path, bad in remaining[:20]:
        print("  STILL:", path, bad[:3])
    sys.exit(1 if (remaining or text_problems) else 0)


if __name__ == "__main__":
    main()
