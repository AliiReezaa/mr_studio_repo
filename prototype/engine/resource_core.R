# ---- Round 66: memory-aware heavy analyses -----------------------------------------------------------------------------------------
# The long genome-wide analyses (LHC-MR, LD score regression, LD score computation) can exhaust the memory of an ordinary laptop.
# Before they start, the app reads the computer's memory and offers only settings whose estimated peak fits in a safe share of the
# memory that is free at that moment; while LHC-MR runs, a guard process stops it cleanly if free memory falls too low.

RES_SAFE_SHARE <- 0.5    # an option is offered only if its estimated peak fits in half of the memory free now
RES_GUARD_FLOOR <- 0.10  # the guard stops a run when free memory falls below 10% of the total (or below 1 GB)

.num1 <- function(x) suppressWarnings(as.numeric(x[1]))

#' Total and available physical memory in GB (NA when the system does not tell). Available = what can be used without swapping:
#' MemAvailable on Linux; on macOS the free percentage of memory_pressure (else free + inactive + speculative + purgeable pages);
#' on Windows the Available memory of the performance counters (as Task Manager shows), else FreePhysicalMemory.
mem_info <- function() {
  tot <- NA_real_; av <- NA_real_
  tryCatch({
    if (.Platform$OS.type == "windows") {
      o <- system2("powershell", c("-NoProfile", "-Command", shQuote("$o = Get-CimInstance Win32_OperatingSystem; $m = Get-CimInstance Win32_PerfFormattedData_PerfOS_Memory; \"$($o.TotalVisibleMemorySize) $($o.FreePhysicalMemory) $($m.AvailableBytes)\"")), stdout = TRUE, stderr = FALSE)
      v <- as.numeric(strsplit(trimws(o[length(o)]), "\\s+")[[1]]); tot <- v[1] / 1048576; av <- if (length(v) >= 3 && is.finite(v[3])) v[3] / 1073741824 else v[2] / 1048576   # Round 100: Available (free + standby cache, as Task Manager shows); FreePhysicalMemory leaves out the standby cache
    } else if (Sys.info()[["sysname"]] == "Darwin") {
      tot <- .num1(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE)) / 1073741824
      vs <- system2("vm_stat", stdout = TRUE); ps <- .num1(sub(".*page size of ([0-9]+) bytes.*", "\\1", vs[1]))
      pg <- function(k) { l <- grep(paste0("^", k, ":"), vs, value = TRUE); if (length(l)) .num1(gsub("[^0-9]", "", sub(".*:", "", l))) else 0 }
      av <- (pg("Pages free") + pg("Pages inactive") + pg("Pages speculative") + pg("Pages purgeable")) * ps / 1073741824
      mp <- tryCatch(system2("memory_pressure", stdout = TRUE, stderr = FALSE), error = function(e) character()); fp <- grep("free percentage", mp, value = TRUE)
      if (length(fp)) av <- tot * .num1(gsub("[^0-9]", "", sub(".*:", "", fp[1]))) / 100   # Apple's own measure (counts reclaimable file cache)
    } else if (file.exists("/proc/meminfo")) {
      m <- readLines("/proc/meminfo"); g <- function(k) .num1(gsub("[^0-9]", "", grep(paste0("^", k, ":"), m, value = TRUE))) / 1048576
      tot <- g("MemTotal"); av <- g("MemAvailable")
    }
  }, error = function(e) NULL)
  list(total_gb = tot, available_gb = av)
}

#' Whether a run needing `need_gb` at its peak may start now. Returns list(ok, msg) with the numbers in the message.
mem_check <- function(need_gb, what, mi = mem_info(), share = RES_SAFE_SHARE) {   # Round 85: share, for an estimate measured on real data
  if (!is.finite(mi$available_gb)) return(list(ok = TRUE, msg = "", unknown = TRUE))
  if (need_gb <= share * mi$available_gb) return(list(ok = TRUE, msg = ""))
  list(ok = FALSE, msg = sprintf("%s needs about %.1f GB of memory at its peak, and %.1f GB is free now (%.0f GB in total). To keep the computer responsive it starts only when about %.1f GB is free: close other programs, then try again.",
                                  what, need_gb, mi$available_gb, mi$total_gb, need_gb / share))
}

#' Start Rscript in the background with the environment it needs: RHOME, because a moved (installed) runtime's Rscript has no valid
#' built-in location (as run_pair_script does), R_LIBS for the same packages, and `vars` (a named list). Variables are set with
#' Sys.setenv for the call and restored, never put into the shell command, so paths with spaces ("MR Studio") are safe.
rscript_background <- function(args, log = NULL, vars = list()) {
  vars <- c(list(RHOME = R.home(), R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)), vars)
  old <- Sys.getenv(names(vars), unset = NA, names = TRUE); do.call(Sys.setenv, vars)
  on.exit(for (k in names(old)) if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, stats::setNames(list(old[[k]]), k)), add = TRUE)
  system2(default_rscript(), args, stdout = if (is.null(log)) FALSE else log, stderr = if (is.null(log)) FALSE else log, wait = FALSE)
}

# ---- LHC-MR estimates (calibrated on the validation fit, see validation/darrous2021_lhcmr/README.md) ---------------------------------
# Peak memory grows with the variants and the workers; time falls with the workers, which share the 250 fits.
# Peak memory = base + main process data + one data copy per worker (forked workers' resident memory, counted in full: conservative).
# Time does not depend on the number of variants: one likelihood evaluation takes about 1 to 2 s on 3,000 to 100,000 variants
# (it works on binned distributions). Measured on the BMI and height validation run (2026-11-20, 2019 Intel Mac, 16 GB): the 50
# starting-point fits took 7.6 worker-minutes each, the 200 jackknife fits 2.0 (they start from the best solution), so a run is about
# 13 worker-hours (6.5 hours on 2 workers, 3.3 on 4), plus about 10 minutes of preparation. Memory: peak 4.46 GB with 2 workers
# (main process 2.2 GB, about 1.1 GB per worker) and 6.0 GB with 4, on about 1 million variants.
LHCMR_CAL <- list(base_gb = 1.2, main_gb_per_mvar = 1.0, per_worker_gb_per_mvar = 1.1, sp_fits = 50L, sp_minutes = 7.6, jk_fits = 200L, jk_minutes = 2.0, prep_hours = 0.2, provisional = FALSE)

lhcmr_estimate <- function(variants, workers, cal = LHCMR_CAL) {
  mv <- variants / 1e6
  list(workers = as.integer(workers), memory_gb = cal$base_gb + cal$main_gb_per_mvar * mv + workers * cal$per_worker_gb_per_mvar * mv, hours = cal$prep_hours + (cal$sp_fits * cal$sp_minutes + cal$jk_fits * cal$jk_minutes) / 60 / workers)
}

#' Settings offered for LHC-MR on this computer: worker counts whose estimated peak fits in the safe share of free memory.
#' Windows runs one worker (the package parallelises by forking). Returns a data frame (possibly empty) and the memory read.
lhcmr_options <- function(variants, mi = mem_info(), cal = LHCMR_CAL, max_workers = max(1L, parallel::detectCores() %/% 2L)) {
  w <- if (.Platform$OS.type == "windows") 1L else unique(pmin(c(1L, 2L, 4L, 6L, 8L), max_workers))
  o <- do.call(rbind, lapply(w, function(k) as.data.frame(lhcmr_estimate(variants, k, cal))))
  o$label <- ifelse(o$workers == 1L, "Gentle", ifelse(o$workers == 2L, "Balanced", "Faster"))
  if (is.finite(mi$available_gb)) o <- o[o$memory_gb <= RES_SAFE_SHARE * mi$available_gb, , drop = FALSE] else o <- o[o$workers == 1L, , drop = FALSE]
  rownames(o) <- NULL; attr(o, "memory") <- mi; o
}

#' Guard: a small separate R process that checks free memory every `every` seconds while process `pid` lives, and stops it (and
#' its forked workers) when free memory falls below the floor, writing the reason into the status file.
mem_guard_start <- function(pid, status_file, every = 30, engine = getOption("mrgui.engine_dir", ""), log = NULL) {
  code <- sprintf(paste0("source(file.path(%s, 'resource_core.R')); pid <- %dL; sf <- %s; alive <- function() if (.Platform$OS.type == 'windows') any(grepl(paste0(' ', pid, ' '), system2('tasklist', c('/FI', shQuote(paste('PID eq', pid))), stdout = TRUE))) else isTRUE(tryCatch(tools::pskill(pid, 0L), error = function(e) FALSE))",
    "; repeat { Sys.sleep(%d); if (!alive()) break; mi <- mem_info(); fl <- suppressWarnings(as.numeric(Sys.getenv('MRSTUDIO_GUARD_FLOOR', RES_GUARD_FLOOR))); if (is.finite(mi$available_gb) && mi$available_gb < max(1, fl * mi$total_gb)) {",
    " if (.Platform$OS.type != 'windows') system(sprintf('pkill -TERM -P %%d', pid)); tools::pskill(pid)",
    "; writeLines(jsonlite::toJSON(list(state = 'stopped', msg = sprintf('Stopped to protect the computer: free memory fell to %%.1f GB of %%.0f GB. Run it again with fewer workers, or close other programs first.', mi$available_gb, mi$total_gb), time = format(Sys.time(), '%%Y-%%m-%%d %%H:%%M:%%S')), auto_unbox = TRUE), sf); break } }"),
    deparse(engine), as.integer(pid), deparse(status_file), as.integer(every))
  rscript_background(c("--vanilla", "-e", shQuote(code)), log = log); invisible(TRUE)
}

# ---- LDSC and LD score computation (measured 2026-11-20 on a 16 GB Intel Mac): GenomicSEM::ldsc on two 1.1 million variant files
# peaked at 1.34 GB; reading the Data library files before it adds about 0.12 GB per million rows of the larger dataset. The LD score
# computation peaked at 0.51 GB on chromosome 2 (one of the largest); PLINK's extraction is a small separate process.
ldsc_memory_gb <- function(rows) 1.5 + 0.12 * max(rows, 0, na.rm = TRUE) / 1e6
LDSCORE_MEMORY_GB <- 1.0

# ---- Round 72c: keep the computer awake during a run, and a system notification when it finishes -----------------------------------
# Both follow the Settings switches (awake, notify; on by default). Nothing is sent online. The keep-awake helper is tied to the engine
# process, so it ends by itself if the app closes; keep_awake_stop() ends it when the run finishes or is stopped.
.awake <- new.env()
keep_awake_start <- function() {
  keep_awake_stop()
  sys <- Sys.info()[["sysname"]]; pid <- Sys.getpid()
  r <- tryCatch(.with_user_env(
    if (sys == "Darwin") {   # caffeinate: no idle sleep while the engine process lives (-w), released when stopped
      as.integer(system2("/bin/sh", c("-c", shQuote(sprintf("/usr/bin/caffeinate -i -w %d >/dev/null 2>&1 & echo $!", pid))), stdout = TRUE))
    } else if (.Platform$OS.type == "windows") {   # SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED) until the engine process ends
      pf <- tempfile(fileext = ".pid")
      ps <- sprintf("$PID | Out-File -Encoding ascii '%s'; Add-Type -Namespace MRS -Name P -MemberDefinition '[DllImport(\"kernel32.dll\")] public static extern uint SetThreadExecutionState(uint f);'; [MRS.P]::SetThreadExecutionState(0x80000001) | Out-Null; Wait-Process -Id %d", normalizePath(pf, winslash = "/", mustWork = FALSE), pid)
      system2("powershell.exe", c("-NoProfile", "-WindowStyle", "Hidden", "-Command", shQuote(ps, type = "cmd")), wait = FALSE, stdout = FALSE, stderr = FALSE)
      for (i in 1:20) { if (file.exists(pf) && length(x <- suppressWarnings(as.integer(readLines(pf, warn = FALSE)))) && is.finite(x[1])) break; Sys.sleep(0.25) }
      if (file.exists(pf)) suppressWarnings(as.integer(readLines(pf, warn = FALSE))[1]) else NA_integer_
    } else NA_integer_), error = function(e) NA_integer_)
  .awake$pid <- if (length(r) && is.finite(r[1])) r[1] else NA_integer_
  invisible(is.finite(.awake$pid))
}
keep_awake_stop <- function() {
  p <- .awake$pid; .awake$pid <- NULL
  if (!is.null(p) && is.finite(p)) tryCatch(tools::pskill(p), error = function(e) NULL)
  invisible(TRUE)
}
#' A system notification (macOS Notification Centre, Windows notification area). Text is MR Studio's own, never from a page.
notify_done <- function(title, msg) {
  # Round 72e: inside the app the desktop shell shows it as MR Studio (Tauri notification); this line on standard output asks it to.
  if (isTRUE(getOption("mrstudio.shell_notify"))) { cat("MRSTUDIO_NOTIFY ", jsonlite::toJSON(list(title = title, body = msg), auto_unbox = TRUE), "\n", sep = ""); flush(stdout()); return(invisible(TRUE)) }
  sys <- Sys.info()[["sysname"]]; clean <- function(x) gsub("[\"'`$\\\\]", "", x)
  tryCatch(.with_user_env(
    if (sys == "Darwin") {
      system2("/usr/bin/osascript", c("-e", shQuote(sprintf("display notification \"%s\" with title \"%s\"", clean(msg), clean(title)))), wait = FALSE, stdout = FALSE, stderr = FALSE)
    } else if (.Platform$OS.type == "windows") {
      ps <- sprintf("Add-Type -AssemblyName System.Windows.Forms; $n = New-Object System.Windows.Forms.NotifyIcon; $n.Icon = [System.Drawing.SystemIcons]::Information; $n.Visible = $true; $n.ShowBalloonTip(8000, '%s', '%s', 'Info'); Start-Sleep -Seconds 9; $n.Dispose()", clean(title), clean(msg))
      system2("powershell.exe", c("-NoProfile", "-WindowStyle", "Hidden", "-Command", shQuote(ps, type = "cmd")), wait = FALSE, stdout = FALSE, stderr = FALSE)
    }), error = function(e) NULL)
  invisible(TRUE)
}
