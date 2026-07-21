# 論文用にベクター(PDF)図を生成。既存の rds ベース図スクリプトを、png() を pdf() に差し替えて source する
# （元スクリプトは無編集＝単体実行では従来どおり PNG を出力）。全て rds から再描画＝再フィット不要。
# Working directory = repository root (all paths below are relative to it). If needed: setwd("/path/to/nmfcox-paper")

## png(file.png, width=W_px, height=H_px, res=R, type=...) を、同じ物理寸法の cairo_pdf(file.pdf, W/R in, H/R in) に転送。
## cairo_pdf はフォントをサブセット埋め込みする（base pdf() は非埋め込み＝雑誌不可）。
png <- function(filename, width, height, res=72, type=NULL, ...)
  grDevices::cairo_pdf(sub("\\.png$", ".pdf", filename), width = width/res, height = height/res)

for (s in c("generate_mono_figs_20260712.R",     # fig 1-4: reassigned_basis_beta / mgus2_HR / dummy_hr / compare_rrr_basis_count
            "generate_gamma_forest_20260712.R",  # fig gamma_forest
            "generate_confound_fig_20260712.R")) # fig confound_forest
  source(s, local = FALSE)

cat("=== PDF figures written to v3_fig/*.pdf ===\n")
