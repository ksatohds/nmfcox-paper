# Porting recipe: paper scripts → nmfkc (@develop) reproduction scripts

The paper's figures/tables were produced by scripts that `source("nmf_cox3.R")` (a frozen core).
Those functions are now in the R package **nmfkc** (develop branch). Verified: the package
`nmf.cox` reproduces `nmf_cox3.R` **bit-identically** (|Δγ|,|Δβ(t)|,|Δobjfunc| < 1e-6).

## The mechanical edits (apply to each script)

1. Replace the source lines
   `source("nmf_cox2.R"); source("nmf_cox3.R")`  (or any `source("nmf_cox3.R")`)
   with
   `suppressMessages(library(nmfkc))`
   Keep `library(survival)` and any other `library(...)` / `source("coxvc.r")` / `source("coxvcf.r")`.

2. Rename function calls:
   - `nmf.cox3(`       → `nmf.cox(`
   - `nmf.cox3.cv(`    → `nmf.cox.cv(`
   - `nmf.cox3.cf(`    → `nmf.cox.cf(`
   - `nmf.cox3.phtest(`→ `nmf.cox.phtest(`

3. Rename the smoothing argument in those calls:
   - `X.smooth=`   → `X.L2.smooth=`
   - in `nmf.cox.cv` only: `smooth.grid=` → `X.L2.smooth=`  and  `rank.grid=` → `rank=`

3b. In scripts that read `nmf.cox.cv` RESULT fields, rename the accessors (the package renamed them):
   - `$best.smooth`     → `$X.L2.smooth.best`
   - `$best.rank`       → `$rank.best`
   - `$best.smooth.1se` → `$X.L2.smooth.best.1se`
   - `$best.rank.1se`   → `$rank.best.1se`
   (`$cvpl` and `$se` are unchanged.)

4. Leave EVERYTHING ELSE unchanged: all other arguments (`ties`, `x.update`, `seed`,
   `inference`, `maxit`, `X.restriction`, `X.init`, `nstart`, `A.normalize`, `nfolds`,
   `mc.cores`, `aic`, ...) pass through `...` and reproduce exactly. Do NOT change data
   handling, seeds, λ_X values, output filenames, or paths. Keep the `setwd(...)` line.

## Exceptions

- `nmf.cox2(...)` (the v2 reference used only in the solver-benchmark `verify_newton_*`):
  nmfkc has no `nmf.cox2`. Keep `source("nmf_cox2.R")` for that ONE v2-reference call, or drop
  the cox2 column with a note. All `nmf.cox3*` calls still port as above.
- `coxvc.r` / `coxvcf.r` (Perperoglou reduced-rank regression, used by the RRR-comparison
  figures) are external code, copied into this folder; keep the `source(...)` lines.

## Verification each ported script must pass

- It must `Rscript` without error.
- Figure/redraw scripts (read an `.rds`, draw a `.pdf/.png`): run to completion; the output is
  bit-identical to the paper figure (deterministic).
- Heavy Monte-Carlo / server scripts (`mc_*`, `verify_newton_*`): do NOT re-run the full 500-rep
  job. Verify the ported `nmf.cox*` call works on a tiny smoke (few reps / one dataset), then
  note "server job; rds provided" in the provenance row.

## Install (for the record / README)

```r
remotes::install_github("ksatohds/nmfkc@develop")
```
