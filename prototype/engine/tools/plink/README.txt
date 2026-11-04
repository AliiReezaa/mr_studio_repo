PLINK 1.9 for clumping on this computer (Round 19).

The build pipeline (.github/workflows/prototype.yml) downloads the official build pinned in prototype/engine/ld_core.R
(PLINK_PIN: version, URL and SHA-256 per operating system), checks the SHA-256, and puts the program and its LICENSE here.
On macOS it is signed ad hoc, because Apple Silicon refuses to run unsigned code.

PLINK 1.9: www.cog-genomics.org/plink/1.9/ (Chang et al. 2015, GigaScience 4:7; Purcell et al. 2007, Am J Hum Genet 81:559).
Licence: GNU General Public License v3 (LICENSE, shipped next to the program). Source code: github.com/chrchang/plink-ng
