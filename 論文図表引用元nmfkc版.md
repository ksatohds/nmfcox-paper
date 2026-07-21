# 論文の図表 ↔ 再現スクリプト対応（nmfkc 版）

本フォルダは、投稿論文 **NMF-COX**（`SiM20260722.tex`）の図表を、
論文用の凍結コア `nmf_cox3.R` ではなく **R パッケージ nmfkc（@develop 版）** から
再現するためのスクリプト一式です。移植は機械的置換のみで、数値は凍結コアと一致します
（`nmf.cox` が `nmf_cox3.R` をビット一致で再現することを検証済み。移植手順は `PORTING_RECIPE.md`）。

```r
remotes::install_github("ksatohds/nmfkc@develop")   # 8/15 頃 CRAN 公開予定
```

- **実行ディレクトリ = このリポジトリの直下**（各スクリプトのパスはすべて直下からの相対）。
- 図の再描画に必要な事前計算結果（server での Monte Carlo 500 反復や多スタート当てはめの出力）は
  `verification/*.rds`・`v3_fig/*.rds` として同梱。これらを削除すれば各生成スクリプトが再計算します
  （重い MC はサーバ `OPENBLAS_NUM_THREADS=1 Rscript ...` 前提）。
- **RRR 比較図（Fig 7・9）** は Perperoglou の縮小ランク回帰コード `coxvc.r` / `coxvcf.r` を使います。
  第三者コードのため本リポジトリには**同梱していません**。<https://github.com/drperpo> から入手し
  リポジトリ直下に置いてください（未配置ならスクリプトが入手先を案内して停止します）。

---

## 図（Figures）

| 論文 | 図ファイル | データ生成スクリプト | 清書・再描画 | データ／設定 | 備考 |
|---|---|---|---|---|---|
| **fig:convergence**（収束） | `v3_fig/convergence.pdf` | `generate_convergence_fig_20260714r2.R`（自己完結：`nmf.cox` を当てはめ `objfunc.iter` を描画） | — | veteran/gbsg/mgus2 | 港：`library(nmfkc)`＋`nmf.cox`＋`X.L2.smooth`。 |
| **fig:betasim**（β̂ シミュ） | `v3_fig/beta_sim.pdf` | `verification/mc_beta_20260712.R`（server MC 500 反復・`nmf.cox.cf`／`nmf.cox(inference=TRUE)` → `mc_beta_beta500B57.rds`） | `generate_beta_sim_20260713B.R`（rds を読んで描画＋LaTeX 行出力） | n=150/300/600/1200, βA=.4, rank1, λ_X=300 | 再描画は純描画（`nmf.cox` 呼び出しなし）。**tab:betasim** も同 rds。 |
| **fig:confound**（交絡） | `v3_fig/confound_forest.pdf` | `generate_confound_fig_20260714.R`（rds があれば再描画、無ければ server MC を実行 → `mc_conf_conf500B58.rds`） | 同左（自己完結） | 500 反復・4 設定（n300/600, ρ=0/.5/.7） | 同梱 rds で即再描画。 |
| **fig:v3basis**（共有非負基底＋β(t)） | `v3_fig/reassigned_basis_beta.pdf` | `generate_reassigned_basis_beta_20260714.R`（`nmf.cox` 多スタート当てはめ → `reassigned_basis_beta_results.rds`） | `generate_mono_figs_20260712.R`（白黒版 PNG）→ `make_pdf_figs_20260712.R`（PNG→PDF） | veteran/gbsg/mgus2、cox.zph 準拠 split、λ_X=1000/3000/3000、多スタート NSTART=8 | mgus2 は全 complete-case N=1338, A={age,hgb,creat}（査読 r1 Major#2）。診断ブロック(i)は旧 rds と次元不一致ならスキップ（図生成には無関係）。 |
| **fig:mgus2hr**（値別 HR） | `v3_fig/mgus2_HR_by_value.pdf` | `generate_mgus2hr_20260714.R`（`reassigned_basis_beta_results.rds` を読む） | `make_pdf_figs_20260712.R` | mgus2 の β(t) を代表値の HR に変換 | 純描画（rds 依存）。 |
| **fig:dummyhr**（ダミー別 HR） | `v3_fig/dummy_hr_by_covariate.pdf` | `generate_mono_figs_20260712.R`（`reassigned_basis_beta_results.rds` の非PH共変量 HR パネル） | `make_pdf_figs_20260712.R` | veteran/gbsg/mgus2 | 純描画（rds 依存）。 |
| **fig:rrrfit**（RRR 基底本数比較） | `v3_fig/compare_rrr_basis_count.pdf` | `compare_rrr_basis_count_20260714.R`（`nmf.cox`＋RRR〔coxvcf〕→ `compare_rrr_basis_count_results.rds`）**［coxvc 必要］** | `regen_rrrfit_fig_20260714.R`（rds を読んで描画） | veteran/gbsg/mgus2, `nmf.cox` rank2/λ_X, RRR reduced-rank r=2, q=2,3,4,6 | 港：`fN$alpha`→`fN$gamma`（nmfkc は `$gamma` のみ）。 |
| **fig:gamma**（γ フォレスト） | `v3_fig/gamma_forest.pdf` | `generate_gamma_table_20260714.R`（`nmf.cox.cf` 3 データ → `gamma_table_20260714.rds/.csv`） | `generate_gamma_forest_20260714.R`（rds を読んで描画） | veteran/gbsg/mgus2, cf rank2, λ_X=1000/3000/3000, nfolds=5, seed=1 | γ 保存（③≈② full-PH, 平均\|Δ\|≈0.027）。 |
| **fig:rankcmp**（rank 別 gbsg） | `v3_fig/rank_compare_gbsg.pdf` | `generate_rank_compare_gbsg_20260713r3.R`（NMF は rds から、RRR〔coxvcf〕を M=5 で再当てはめ → `rank_compare_gbsg_results.rds`）**［coxvc 必要］** | `generate_rank_compare_gbsg_mono_20260714.R`（白黒版） | gbsg, NMF Q=2,3,4／RRR M=5 anchored＋M=6, reduced rank r=2,3,4 | r3 は NMF 呼び出しなし（rds 読込）ゆえ移植不要。 |

## 表（Tables）

| 論文 | 生成スクリプト | データ／設定 | 備考 |
|---|---|---|---|
| **tab:mcsim**（クロスフィット MC） | `verification/mc_cf_oracle_20260712.R`（基準行）＋ `verification/mc_overrank_R2_20260714.R`（over-rank 行） | server MC 500 反復。基準：4 設定（n300/600, βA=.4/.8, rank1/2）, λ_X=300, nfolds=5。over-rank：R=2/true-rank-1/Q=2, n=300, λ_X=300 | `nmf.cox.cf` 移植済み。バイアス≈0, カバレッジ 0.93–0.94, 対 stratified oracle 差 0.004–0.011。 |
| **tab:betasim**（β̂ MC） | `verification/mc_beta_20260712.R` → `mc_beta_beta500B57.rds`（fig:betasim と共有） | 上記 fig:betasim と同 | `nmf.cox(inference=TRUE)`。 |
| **tab:select**（rank/λ_X 選択） | `verification/mc_select_R3_20260713.R` | server MC。R=3, rank grid 1:3, λ_X grid c(30,100,300,1000), nfolds=5, aic/bic/cvpl | `nmf.cox.cv` 移植＋結果フィールド改名（`$rank.best`, `$X.L2.smooth.best`, `.1se`）。 |
| **tab:newton**（ソルバ堅牢性） | `verification/verify_newton_alldata_20260714.R` | 7 データ（sim＋veteran/gbsg/pharma/lung/pbc/mgus2）, rank2, λ_X 300–3000。newton vs optim vs cox2 | v2 参照 `nmf_cox2.R` を同梱（cox2 列のみ）。港：`fo/fn$alpha`→`$gamma`, `f2`（=nmf.cox2）は `$alpha` 維持。asaur 未導入なら pharma 行はスキップ。 |
| **λ_X 選択（mgus2, 本文）** | `verification/cvpl_mgus2_allcc_20260714.R` | mgus2 全 complete-case, rank2, λ_X grid, cvpl, nfolds=5 | λ_X=3000 が CVPL 1-SE 帯内であることを確認。 |
| tab:methods / tab:data / tab:perp | —（本文の定性表・データ記述） | — | 計算を伴わないため再現スクリプトなし。 |

## 補足：移植で必要だった調整（`PORTING_RECIPE.md` に加えて）

- **`$alpha` 別名の廃止**：nmfkc の `nmf.cox` は PH 係数を `$gamma` のみで返します（凍結コアの後方互換別名 `$alpha` は無し）。
  `compare_rrr_basis_count`（`fN$alpha`→`fN$gamma`）と `verify_newton`（`fo/fn$alpha`→`$gamma`; ただし `nmf.cox2` の `f2$alpha` は維持）で対応。
  `nmf.cox.cf` はもとから `$gamma` を返すためフォレスト/表側は変更不要。
- **cv 結果フィールドの改名**：`$best.smooth`→`$X.L2.smooth.best` 等（`PORTING_RECIPE.md` §3b）。`mc_select`・`cvpl_mgus2` に適用。
- **重い MC/検証（`mc_*` / `verify_newton`）は本ローカルでは smoke のみ確認**。公表数値は server 500 反復の出力（同梱 rds／再計算はサーバ推奨）。
