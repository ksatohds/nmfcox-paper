# NMF-COX: Reproducibility Materials

Reproducibility code for:

> **Separating proportional and non-proportional hazards with a non-negative time basis: interpretable time-varying effects and cross-fit oracle inference in the Cox model**
> Kenichi Satoh (Shiga University), 2026

The method separates a constant proportional-hazards effect **γ** for covariates *z* from a
low-rank, non-negative time-varying structure **β(t) = XΘ** for non-proportional covariates *a*
(a semi-NMF: the time basis *X ≥ 0* with column sums one; the loading *Θ* sign-free), and
estimates γ on the partial likelihood with cross-fit (Neyman-orthogonal) oracle inference.

All figures/tables were produced with the R package **`nmfkc`** (develop branch). The frozen
core `nmf_cox3.R` used during development is reproduced bit-for-bit by `nmfkc::nmf.cox`
(see `PORTING_RECIPE.md`), so these scripts and the paper give identical numbers.

## Requirements

### R
- R ≥ 4.3.0 (developed on R 4.4.1).
- Linux recommended for the Monte-Carlo scripts (`parallel::mclapply` uses fork); on Windows they
  run serially with identical results. Run heavy MC on a compute server with `OPENBLAS_NUM_THREADS=1`.

### R packages
```r
# Core: NMF-COX is integrated into nmfkc (nmf.cox / nmf.cox.cv / nmf.cox.cf / nmf.cox.phtest / nmf.cox.inference)
remotes::install_github("ksatohds/nmfkc@develop")   # CRAN release planned Aug 2026
install.packages("survival")                        # datasets (veteran, gbsg, mgus2, lung, pbc) + Cox
# optional: install.packages("asaur")               # pharmacoSmoking, used only by tab:newton
```

### Third-party code (RRR comparison, Fig 7 & 9 only)
The reduced-rank regression (RRR) baseline uses A. Perperoglou's `coxvc.r` / `coxvcf.r`. These are
**not redistributed here**. Download them from <https://github.com/drperpo> and place them in the
repository root. Scripts that need them stop with a message if they are missing.

## How to run
Set the working directory to the repository root (all paths are relative to it), then run any script:
```bash
Rscript generate_gamma_forest_20260714.R
Rscript verification/verify_newton_alldata_20260714.R   # server job
```
Figures write to `v3_fig/`. Pre-computed Monte-Carlo results and multi-start fits are shipped as
`verification/*.rds` and `v3_fig/*.rds`; delete them to regenerate from scratch (the `mc_*` and
`verify_*` scripts are 500-replication server jobs).

## Repository structure

| File | Paper element | Description |
|------|---------------|-------------|
| `generate_convergence_fig_20260714r2.R` | Fig `convergence` | Objective-function convergence (`nmf.cox`, veteran/gbsg/mgus2) |
| `verification/mc_beta_20260712.R` → `generate_beta_sim_20260713B.R` | Fig `betasim`, Table `betasim` | β̂ simulation MC (500 reps) + redraw |
| `generate_confound_fig_20260714.R` | Fig `confound` | Confounding simulation forest (self-contained: redraw or regenerate) |
| `generate_reassigned_basis_beta_20260714.R` → `generate_mono_figs_20260712.R` → `make_pdf_figs_20260712.R` | Fig `v3basis` | Shared non-negative basis *X* + β(t) of non-PH covariates |
| `generate_mgus2hr_20260714.R` | Fig `mgus2hr` | mgus2 value-specific hazard ratios |
| `generate_mono_figs_20260712.R` | Fig `dummyhr` | Dummy-covariate HR panel (mono) |
| `compare_rrr_basis_count_20260714.R` → `regen_rrrfit_fig_20260714.R` | Fig `rrrfit` | NMF-COX vs RRR by basis count *(needs coxvc)* |
| `generate_gamma_table_20260714.R` → `generate_gamma_forest_20260714.R` | Fig `gamma` | Cross-fit γ̂ forest (γ preserved vs full-PH Cox) |
| `generate_rank_compare_gbsg_20260713r3.R` → `generate_rank_compare_gbsg_mono_20260714.R` | Fig `rankcmp` | Rank comparison on gbsg, NMF-COX vs RRR *(needs coxvc)* |
| `verification/mc_cf_oracle_20260712.R`, `verification/mc_overrank_R2_20260714.R` | Table `mcsim` | Cross-fit oracle MC (baseline + over-rank rows) |
| `verification/mc_select_R3_20260713.R` | Table `select` | Rank / λ_X selection MC (AIC/BIC/CVPL) |
| `verification/verify_newton_alldata_20260714.R` | Table `newton` | Block-Newton solver robustness across 7 datasets |
| `verification/cvpl_mgus2_allcc_20260714.R` | §illustration | CVPL selection of λ_X on mgus2 (all complete cases) |
| `nmf_cox2.R` | — | Frozen v2 reference implementation (used only for the cox2 column of Table `newton`) |
| `PORTING_RECIPE.md` | — | How the `nmf_cox3.R` scripts map to `nmfkc` |
| `論文図表引用元nmfkc版.md` | — | Detailed figure/table ↔ script provenance (Japanese) |

`nmfkc::nmf.cox` returns the PH coefficient as `$gamma` (the frozen core's back-compat alias
`$alpha` was dropped); the scripts here read `$gamma` accordingly. The v2 reference `nmf_cox2.R`
still returns `$alpha`.

## Citation
If you use this code, please cite the paper above and the `nmfkc` package
(<https://github.com/ksatohds/nmfkc>).

## License
MIT (see `LICENSE`). Note: the third-party `coxvc.r` / `coxvcf.r` are **not** included and carry
their own terms from <https://github.com/drperpo>.
