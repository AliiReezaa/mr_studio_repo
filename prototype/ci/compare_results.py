"""Cross-OS reproducibility check (contracts C8.6).
Compares estimates and SEs of every method result across operating systems.
Usage: python3 compare_results.py <reports_dir>"""
import glob, json, os, sys

TOL = 1e-8
root = sys.argv[1]
by_os = {}
for d in sorted(glob.glob(os.path.join(root, "report-*"))):
    os_id = os.path.basename(d).replace("report-", "")
    files = glob.glob(os.path.join(d, "**", "result_*.json"), recursive=True)
    by_os[os_id] = {os.path.basename(f): json.load(open(f)) for f in files}

if len(by_os) < 2:
    print("fewer than two operating systems reported; nothing to compare"); sys.exit(1)

ref_os = sorted(by_os)[0]
failures = 0
for name, ref in by_os[ref_os].items():
    for other, results in by_os.items():
        if other == ref_os:
            continue
        res = results.get(name)
        if res is None:
            print(f"FAIL {name}: missing on {other}"); failures += 1; continue
        if res["fingerprint"] != ref["fingerprint"] and res["provenance"]["os"] == ref["provenance"]["os"]:
            print(f"FAIL {name}: fingerprint differs on the same OS family")
            failures += 1
        for a, b in zip(ref.get("estimates") or [], res.get("estimates") or []):
            for key in ("estimate", "se", "ci_lower", "ci_upper"):
                x, y = a[key], b[key]
                if x is None and y is None:
                    continue
                rel = abs(x - y) / max(1.0, abs(x))
                status = "PASS" if rel <= TOL else "FAIL"
                failures += status == "FAIL"
                print(f"{status} {name} {a['estimand']} {key}: {ref_os}={x!r} {other}={y!r} rel_diff={rel:.2e}")
print(f"\n{failures} failure(s) across {len(by_os)} operating systems")
sys.exit(1 if failures else 0)
