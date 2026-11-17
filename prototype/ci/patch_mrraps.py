"""Remove the rsnps import from the pinned mr.raps source before it is installed into the bundled runtime (Round 38).

mr.raps imports rsnps only for the optional online annotation of its own scatter plot
(mr.raps.scatterplot with annotate = TRUE, which queries NCBI). MR Studio never calls it.
rsnps cannot be installed from conda-forge on Apple silicon (its dependency RJSONIO has no
osx-arm64 build), so the import is removed on every system alike. The estimator code is unchanged.
Usage: python patch_mrraps.py <mr.raps source folder>
"""
import pathlib
import re
import sys

src = pathlib.Path(sys.argv[1])
desc = (src / "DESCRIPTION").read_text(encoding="utf-8")
new = re.sub(r"\n[ \t]*rsnps,", "", desc)
new = re.sub(r"\nRemotes:\n[ \t]*ropensci/rsnps[ \t]*", "", new)
ns = (src / "NAMESPACE").read_text(encoding="utf-8")
new_ns = ns.replace("import(rsnps)\n", "")
if new == desc or new_ns == ns or "rsnps" in new or "rsnps" in new_ns:
    sys.exit("patch_mrraps: the rsnps import was not found as expected; check the pinned mr.raps source")
(src / "DESCRIPTION").write_text(new, encoding="utf-8")
(src / "NAMESPACE").write_text(new_ns, encoding="utf-8")
print("patch_mrraps: rsnps import removed (used only by the optional plot annotation)")
