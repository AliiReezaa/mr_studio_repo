# Week-1 prototype

Goal: prove, before building the real app, that a one-click desktop app can carry its own R engine on macOS and Windows, that the window-to-engine connection is secure, and that results are identical across operating systems.

## What is here

| Folder | Contents |
|---|---|
| `engine/` | R engine entry point (`main.R`), connection gate (`gate.R`) and its attack/use test (`test_gate.py`) |
| `shell/` | Tauri 2 desktop shell (Rust): starts the engine, passes secrets on stdin, opens the one-time login URL, stops the engine on exit; `MRGUI_SMOKE=1` runs a self-test and exits |
| `runtime/` | Empty in the repository; the pipeline fills it with the relocatable conda-forge R runtime |
| `ci/` | Runtime definition (`runtime.yml`) and the cross-OS result comparison |
| `../.github/workflows/prototype.yml` | Pipeline for macOS arm64, macOS Intel and Windows |

## Results so far (Linux, 2026-10-29)

| Exit criterion | Status |
|---|---|
| Relocatable R runtime: built, moved to a folder with a space in its name, still runs | Passed on Linux. `Rscript` is not relocatable; the shell starts `lib/R/bin/exec/R` with `R_HOME` set. |
| Contract tests (37) run with the moved runtime | 37/37 pass |
| Connection gate (C10): one-time token, HttpOnly cookie, Host and Origin checks, WebSocket gating, full Shiny session | 16/16 pass |
| Desktop shell compiles and opens the engine page | Passed on Linux (Tauri 2.12, Rust 1.97); window showed the live Shiny session |
| Engine stops when the shell exits normally and when the shell is force-killed | Passed on Linux (engine polls its parent; Windows uses a job object, type-checked only) |
| Smoke test through the shell (`MRGUI_SMOKE=1`) | 4/4 pass, exit code 0 |
| macOS arm64, macOS Intel, Windows: all of the above plus installers | **Pending: needs the pipeline run** |
| Same estimates on every OS within 1e-8 (C8.6) | **Pending: compared automatically by the pipeline** |

Findings recorded as decisions: bundle **R 4.5.3** (R 4.6 lacks conda-forge builds of some packages); the runtime is about 1.5 GB before trimming, so trimming is a v1 task to meet the installer size target.

Bugs found and fixed during the prototype: a token comparison that accepted any token (`xor()` is logical in R; now `bitwXor()`), a WebSocket gate that blocked Shiny (handlers must return `NULL` to pass), and secrets lost to lazy evaluation (`force()`).

## How to run the pipeline

1. Create a GitHub repository and upload this whole folder, including the hidden `.github` folder.
2. Open the repository's **Actions** tab, choose **prototype**, then **Run workflow**.
3. When it finishes (about 30 to 60 minutes), download the `installer-*` artifacts and test them: on a Mac, open the `.dmg` and drag the app to Applications (first launch: right-click, then Open); on Windows, run the `.exe` (SmartScreen: More info, then Run anyway).

Cost: build minutes are free for public repositories. For private repositories on GitHub Free the monthly allowance covers roughly 190 macOS minutes, so a private repository can afford only one or two full runs a month.
