# ---- Round 66: memory-aware heavy analyses -----------------------------------------------------------------------------------------
# The long genome-wide analyses (LHC-MR, LD score regression, LD score computation) can exhaust the memory of an ordinary laptop.
# Before they start, the app reads the computer's memory and offers only settings whose estimated peak fits in a safe share of the
# memory that is free at that moment; while LHC-MR runs, a guard process stops it cleanly if free memory falls too low.

RES_SAFE_SHARE <- 0.5    # an option is offered only if its estimated peak fits in half of the memory free now
RES_GUARD_FLOOR <- 0.10  # the guard stops a run when free memory falls below 10% of the total (or below 1 GB)

.num1 <- function(x) suppressWarnings(as.numeric(x[1]))

#' Total and available physical memory in GB (NA when the system does not tell). Available = what can be used without swapping:
#' MemAvailable on Linux; on macOS the free percentage of memory_pressure (else free + inactive + speculative + purgeable pages);
#' FreePhysicalMemory on Windows.
mem_info <- function() {
  tot <- NA_real_; av <- NA_real_
  tryCatch({
    if (.Platform$OS.type == "windows") {
      o <- system2("powershell", c("-NoProfile", "-Command", shQuote("$o = Get-CimInstance Win32_OperatingSystem; \"$($o.TotalVisibleMemorySize) $($o.FreePhysicalMemory)\"")), stdout = TRUE, stderr = FALSE)
      v <- as.numeric(strsplit(trimws(o[length(o)]), "\\s+")[[1]]); tot <- v[1] / 1048576; av <- v[2] / 1048576
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
mem_check <- function(need_gb, what, mi = mem_info()) {
  if (!is.finite(mi$available_gb)) return(list(ok = TRUE, msg = "", unknown = TRUE))
  if (need_gb <= RES_SAFE_SHARE * mi$available_gb) return(list(ok = TRUE, msg = ""))
  list(ok = FALSE, msg = sprintf("%s needs about %.1f GB of memory at its peak, and %.1f GB is free now (%.0f GB in total). To keep the computer responsive it starts only when about %.1f GB is free: close other programs, then try again.",
                                  what, need_gb, mi$available_gb, mi$total_gb, need_gb / RES_SAFE_SHARE))
}

# ---- LHC-MR estimates (calibrated on the validation fit, see validation/darrous2021_lhcmr/README.md) ---------------------------------
# Peak memory grows with the variants and the workers; time falls with the workers, which share the 250 fits.
# Peak memory = base + main process data + one data copy per worker (forked workers' resident memory, counted in full: conservative).
# Time does not depend on the number of variants: one likelihood evaluation takes about 1 to 2 s on 3,000 to 100,000 variants
# (it works on binned distributions; measured 2026-11-20 on a 2019 Intel Mac), and each fit needs a few hundred to 5,000 of them
# (Nelder-Mead, maxit 5000), so about 20 to 100 minutes per fit. A run is 50 starting-point fits and 200 jackknife fits.
# Measured 2026-11-20 (2019 Intel Mac, 16 GB): the package fit on the 1-million-variant validation pair peaked at 4.46 GB with 2 workers
# (main process 2.2 GB, about 1.1 GB per worker; 188 samples); fits took about 22 minutes each (run_equivalence.R). The time per fit is
# confirmed on the full validation run (validation/darrous2021_lhcmr/README.md).
LHCMR_CAL <- list(base_gb = 1.2, main_gb_per_mvar = 1.0, per_worker_gb_per_mvar = 1.1, minutes_per_fit = 22, fits = 250L, provisional = TRUE)

lhcmr_estimate <- function(variants, workers, cal = LHCMR_CAL) {
  mv <- variants / 1e6
  list(workers = as.integer(workers), memory_gb = cal$base_gb + cal$main_gb_per_mvar * mv + workers * cal$per_worker_gb_per_mvar * mv, hours = cal$fits * cal$minutes_per_fit / 60 / workers)
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
mem_guard_start <- function(pid, status_file, every = 30, engine = getOption("mrgui.engine_dir", "")) {
  code <- sprintf(paste0("source(file.path(%s, 'resource_core.R')); pid <- %dL; sf <- %s; alive <- function() if (.Platform$OS.type == 'windows') any(grepl(paste0(' ', pid, ' '), system2('tasklist', c('/FI', shQuote(paste('PID eq', pid))), stdout = TRUE))) else isTRUE(tryCatch(tools::pskill(pid, 0L), error = function(e) FALSE))",
    "; repeat { Sys.sleep(%d); if (!alive()) break; mi <- mem_info(); fl <- suppressWarnings(as.numeric(Sys.getenv('MRSTUDIO_GUARD_FLOOR', RES_GUARD_FLOOR))); if (is.finite(mi$available_gb) && mi$available_gb < max(1, fl * mi$total_gb)) {",
    " if (.Platform$OS.type != 'windows') system(sprintf('pkill -TERM -P %%d', pid)); tools::pskill(pid)",
    "; writeLines(jsonlite::toJSON(list(state = 'stopped', msg = sprintf('Stopped to protect the computer: free memory fell to %%.1f GB of %%.0f GB. Run it again with fewer workers, or close other programs first.', mi$available_gb, mi$total_gb), time = format(Sys.time(), '%%Y-%%m-%%d %%H:%%M:%%S')), auto_unbox = TRUE), sf); break } }"),
    deparse(engine), as.integer(pid), deparse(status_file), as.integer(every))
  system2(default_rscript(), c("--vanilla", "-e", shQuote(code)), stdout = FALSE, stderr = FALSE, wait = FALSE); invisible(TRUE)
}

# ---- LDSC and LD score computation (measured 2026-11-20 on a 16 GB Intel Mac): GenomicSEM::ldsc on two 1.1 million variant files
# peaked at 1.34 GB; reading the Data library files before it adds about 0.12 GB per million rows of the larger dataset. The LD score
# computation peaked at 0.51 GB on chromosome 2 (one of the largest); PLINK's extraction is a small separate process.
ldsc_memory_gb <- function(rows) 1.5 + 0.12 * max(rows, 0, na.rm = TRUE) / 1e6
LDSCORE_MEMORY_GB <- 1.0
