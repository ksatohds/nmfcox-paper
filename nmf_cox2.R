# =====================================================================
# nmf_cox2.R  —  NMF-COX v2 (Path B redesign, session prototype 2026-07-10)
#
# 目的: v1 の二重計上（W=exp(z'α)/S を分解 → offset が z'α を再計上 → α/2 縮小）を
#   根本解決する。offset を「z から作った W」ではなく、「a 由来の時間変動修飾」として
#   直接モデル化し、部分尤度上で推定する。
#
# モデル（log-link, 時間変動係数の非負低ランク表現）:
#   λ_i(t) = h0(t) · exp( z_i'α + o_i(t) ),   o_i(t_j) = x_j' Θ a_i = (XΘA)_{ji}
#   X ∈ R_{≥0}^{K×Q}（時間基底, 非負）,  Θ ∈ R^{Q×L}（符号自由）,  A: L×n
#   ⇒ 時間変動係数 β(t_j) = (XΘ)_{j·} ∈ R^L（covariate a の時間変動効果）
#
# 推定（交互最適化, いずれも部分尤度 = Breslow）:
#   (α,Θ)-block: coxph(Surv(start,stop,event) ~ z + F) ,  F_{row,(q,l)} = X[j,q]·a[l,i]
#                （o = z'α + XΘA が (α,Θ) に線形 → 1本の凸 coxph で同時推定, Θ符号自由）
#   X-block:     given (α,Θ), s_i = Θ a_i (Q×n) ;  maximize
#                Σ_j [ x_j's_(j) − log Σ_{l∈R_j} exp(z_l'α + x_j's_l) ] − (λ/2)Σ||x_j−x_{j-1}||²
#                over vec(X) ≥ 0  by L-BFGS-B（解析勾配 = Schoenfeld 型）。
#
# ★identity-link KL-NMF ではなく log-link 部分尤度なので X≥0・Θ符号自由が両立する。
#
# 返り値: alpha, X, Theta, beta.t(=XΘ, K×L), cox.fit, offset, alpha.history, loglik ほか。
# 依存: survival。日付サフィックスなし・上書き可（コア関数）。
# =====================================================================

suppressMessages(require(survival))

nmf.cox2 <- function(formula, data, A, rank = 2,
                     X.smooth = 0, smooth.time = "gap", ties = "breslow",
                     maxit = 50, tol = 1e-4, relax = 1,
                     robust = TRUE, verbose = TRUE, X.init = NULL,
                     inference = TRUE, patience = 5, nonneg = TRUE) {
  ## nonneg=FALSE は X 非負制約を外した sign-free reduced-rank 版（#4 比較用; 通常の NMF-COX は nonneg=TRUE）。
  ## X-block は Breslow 部分尤度（D_j 全体）のみ実装。Efron は未実装ゆえ拒否する（#1 査読 r2 対応 2026-07-11）。
  if (!identical(ties, "breslow"))
    stop("nmf.cox2: only ties='breslow' is supported (Efron is not implemented in the X-block).")
  Q <- rank
  .eps <- 1e-10

  ## ---- 1. 応答と Cox 共変量 z ----
  mf <- model.frame(formula, data)
  y  <- model.response(mf)
  if (!inherits(y, "Surv")) stop("LHS must be a Surv() object.")
  time <- as.numeric(y[, 1]); status <- as.integer(y[, 2])
  n <- length(time)
  Z <- model.matrix(attr(mf, "terms"), mf)
  Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  p <- ncol(Z); znames <- colnames(Z)

  ## ---- 2. NMF 共変量 A（L×n; 符号自由）----
  if (inherits(A, "formula")) {
    Amf <- model.frame(A, data)
    Amat <- model.matrix(A, Amf)
    Amat <- Amat[, colnames(Amat) != "(Intercept)", drop = FALSE]
    A <- t(Amat)                       # L×n
  } else {
    A <- as.matrix(A)
    if (ncol(A) != n && nrow(A) == n) A <- t(A)   # 受理: n×L も
    if (ncol(A) != n) stop("A must have n columns (L×n).")
  }
  L <- nrow(A); anames <- rownames(A); if (is.null(anames)) anames <- paste0("A", 1:L)
  ## A 行を標準化（中心化＝時点定数で不変, スケール＝Θ が吸収; 数値安定化）
  A.center <- rowMeans(A); A.scale <- apply(A, 1, sd); A.scale[A.scale < .eps] <- 1
  A <- (A - A.center) / A.scale

  ## ---- 3. イベント時刻・リスク集合展開（counting process）----
  et <- sort(unique(time[status == 1])); K <- length(et)
  s_prev <- c(0, et[-K])
  rid <- rj <- rstart <- rstop <- rev <- integer(0)
  for (j in seq_len(K)) {
    risk <- which(time >= et[j])
    m <- length(risk)
    rid   <- c(rid, risk)
    rj    <- c(rj, rep(j, m))
    rstart<- c(rstart, rep(s_prev[j], m))
    rstop <- c(rstop, rep(et[j], m))
    rev   <- c(rev, as.integer(status[risk] == 1 & time[risk] == et[j]))
  }
  rstart <- as.numeric(rstart); rstop <- as.numeric(rstop)
  Zexp <- Z[rid, , drop = FALSE]                    # (Nexp × p)
  ## 各イベント時刻 j のリスク集合行インデックス
  rows_by_j <- split(seq_along(rj), rj)

  ## ---- 4. 初期化 ----
  # X: 非負・列は相異なる基底（定数 + 時間ランプ）。β(t)=XΘ を柔軟に。
  if (is.null(X.init)) {
    tt <- (et - min(et)) / (max(et) - min(et) + .eps)   # [0,1] 時間
    X <- matrix(0, K, Q)
    X[, 1] <- 1
    if (Q >= 2) for (q in 2:Q) X[, q] <- tt^(q - 1)
    ## 退化回避（列を相異させる微小摂動）。決定的＝再現性のため runif ではなく余弦を使用。
    X <- X + 1e-3 * outer(tt, seq_len(Q), function(t, q) cos(q * pi * t))
  } else X <- X.init
  ## 列 RMS=1 に正規化（エントリを O(1) に保ち Θ の巨大化を防ぐ; スケールは Θ へ）
  normX <- function(M) sweep(M, 2, pmax(sqrt(colSums(M^2) / K), .eps), "/")
  X <- normX(X)
  alpha <- tryCatch(coef(coxph(y ~ Z, ties = ties)), error = function(e) rep(0, p))
  names(alpha) <- znames
  Theta <- matrix(0, Q, L)

  ## ---- 平滑化の時間間隔重み（不規則イベント時刻対応; #4）----
  ## 粗さ罰則 Σ_j w_j ||(x_j-x_{j-1})Θ||² を ∫||β'(τ)||²dτ の離散化に合わせる:
  ## w_j ∝ 1/(τ_j-τ_{j-1})。"uniform"=等間隔(=rank; 従来挙動), "gap"=正規化イベント時刻,
  ## "log"=log-time（後期の疎な間隔を相対的に重視）。平均1に正規化して λ スケールを保つ。
  sm.tt <- switch(smooth.time,
    uniform = as.numeric(seq_len(K)),
    gap     = (et - min(et)) / (max(et) - min(et) + .eps),
    log     = { lg <- log(pmax(et, .eps)); (lg - min(lg)) / (max(lg) - min(lg) + .eps) },
    stop("smooth.time must be 'gap', 'log', or 'uniform'"))
  sm.gap <- diff(sm.tt); if (length(sm.gap)) sm.gap <- pmax(sm.gap, max(sm.gap) * 1e-3)  # 微小間隔の暴走回避
  sm.w <- if (length(sm.gap)) { w <- 1 / sm.gap; w / mean(w) } else numeric(0)           # 平均1に正規化

  ## ---- X-block の目的/勾配（負の部分対数尤度 + 平滑化, vec(X)≥0）----
  Sfun <- function(Th) Th %*% A                       # Q×n : s_i
  negll_grad_X <- function(xvec, S, eta0, Msmooth) {
    Xc <- matrix(xvec, K, Q)
    negll <- 0; g <- matrix(0, K, Q)
    for (j in seq_len(K)) {
      rows <- rows_by_j[[j]]
      ids  <- rid[rows]
      Sj   <- S[, ids, drop = FALSE]                  # Q×m
      eta  <- eta0[rows] + as.vector(crossprod(Sj, Xc[j, ]))  # length m
      mx <- max(eta); w <- exp(eta - mx); sw <- sum(w); wn <- w / sw
      dfail <- which(rev[rows] == 1); dj <- length(dfail)       # 同時イベント集合 D_j と数 d_j
      # Breslow 部分尤度寄与: Σ_{i∈D_j} eta_i - d_j·log Σ_{l∈R_j} exp(eta_l)
      # （d_j=1 の無タイ時は従来と一致；#1 査読対応 2026-07-11）
      negll <- negll - (sum(eta[dfail]) - dj * (mx + log(sw)))
      # grad wrt x_j (Q): -(Σ_{i∈D_j} s_i - d_j Σ_l w_jl s_l)
      sbar <- as.vector(Sj %*% wn)
      g[j, ] <- g[j, ] - (rowSums(Sj[, dfail, drop = FALSE]) - dj * sbar)
    }
    # 平滑化: β(t)=XΘ の粗さ (λ/2)Σ_{j≥2}||(x_j-x_{j-1})Θ||² = (λ/2)tr(D M D'), M=ΘΘ'
    # ⇒ 解釈対象 β(t) を直接滑らかにし、Θ の零空間方向（β に無影響）の揺れは罰しない
    if (X.smooth > 0) {
      D  <- diff(Xc)                                  # (K-1)×Q
      DM <- (D %*% Msmooth) * sm.w                     # (K-1)×Q, 各行を時間間隔重み w_j 倍（#4）
      negll <- negll + 0.5 * X.smooth * sum(DM * D)
      gp <- matrix(0, K, Q)
      gp[1:(K-1), ] <- gp[1:(K-1), ] - DM
      gp[2:K, ]     <- gp[2:K, ]     + DM
      g <- g + X.smooth * gp
    }
    list(value = negll, grad = as.vector(g))
  }

  ## ---- 5. 交互反復 ----
  ## 注意: (α,Θ)-block は罰則なし部分尤度, X-block は罰則付き＝厳密には単一目的の座標上昇ではない。
  ##       そこで「罰則付き部分対数尤度の最良反復を保持」し, プラトー(改善停止)で停止, 発散はハルトする。
  alpha.history <- alpha
  converged <- FALSE
  conv.msg  <- sprintf("reached maxit=%d without meeting tol/plateau", maxit)
  opt.conv  <- NA_integer_        # 最終 X-block optim の収束コード（0=収束）
  obj_old   <- NA_real_           # 罰則付き部分対数尤度（前反復）
  beta_old  <- NULL               # β=XΘ（前反復）
  best.obj  <- -Inf               # これまでの最良（最大）罰則付き部分対数尤度
  best.alpha<- alpha; best.Theta <- NULL; best.X <- X; best.it <- 0L
  stall     <- 0L                 # 意味ある改善がない連続反復数
  for (it in seq_len(maxit)) {
    alpha_old <- alpha
    ## (a) (α,Θ)-block: coxph(~ z + F), F_{row,(q,l)} = X[rj,q]·a[l,rid]
    Fmat <- matrix(0, length(rj), Q * L)
    cn <- character(Q * L); col <- 1
    for (q in 1:Q) for (l in 1:L) {
      Fmat[, col] <- X[rj, q] * A[l, rid]
      cn[col] <- paste0("th_", q, "_", l); col <- col + 1
    }
    colnames(Fmat) <- cn
    dat <- data.frame(start = rstart, stop = rstop, ev = rev, id = rid)
    des <- cbind(Zexp, Fmat)
    fit <- tryCatch(
      coxph(Surv(start, stop, ev) ~ des + cluster(id), data = dat, ties = ties),
      error = function(e) NULL)
    if (is.null(fit)) {                        # P1b: (α,Θ)-block 失敗を明示
      if (it == 1)
        stop("nmf.cox2: (alpha,Theta)-block coxph failed at the first iteration; cannot produce a solution.")
      conv.msg <- sprintf("(alpha,Theta)-block coxph failed at iter %d; returning last successful state", it)
      converged <- FALSE; break
    }
    cf <- coef(fit)
    cf[is.na(cf)] <- 0
    alpha_new <- cf[1:p]; names(alpha_new) <- znames
    Theta <- matrix(cf[(p + 1):(p + Q * L)], Q, L, byrow = TRUE)  # rows q, cols l

    ## (b) X-block: L-BFGS-B on vec(X) ≥ 0
    S <- Sfun(Theta)                                   # Q×n
    eta0 <- as.vector(Zexp %*% alpha_new)
    Msm <- tcrossprod(Theta)                           # ΘΘ' (Q×Q) : β(t) 平滑化計量
    opt <- tryCatch(optim(par = as.vector(X), method = "L-BFGS-B",
                 lower = if (nonneg) 0 else -Inf,
                 fn = function(v) negll_grad_X(v, S, eta0, Msm)$value,
                 gr = function(v) negll_grad_X(v, S, eta0, Msm)$grad,
                 control = list(maxit = 100)),
             error = function(e) NULL)
    if (is.null(opt)) {                        # P1b: X-block 失敗を明示
      if (it == 1)
        stop("nmf.cox2: X-block L-BFGS-B failed at the first iteration; cannot produce a solution.")
      conv.msg <- sprintf("X-block L-BFGS-B failed at iter %d; returning last successful state", it)
      converged <- FALSE; break
    }
    X <- matrix(opt$par, K, Q)
    opt.conv <- opt$convergence                # 0=収束, 1=maxit, 51/52=line-search 異常
    ## 識別性: X 列 RMS 正規化（スケールを Θ に移す）
    sc <- pmax(sqrt(colSums(X^2) / K), .eps)
    X <- sweep(X, 2, sc, "/"); Theta <- sweep(Theta, 1, sc, "*")

    ## (c) 緩和・収束（P1a: α だけでなく β=XΘ と罰則付き部分対数尤度の相対変化を併用）
    alpha <- (1 - relax) * alpha_old + relax * alpha_new
    alpha.history <- rbind(alpha.history, alpha)
    beta_new <- X %*% Theta                              # 再スケール後（β は不変）
    obj_new  <- -opt$value                               # 罰則付き部分対数尤度 = ℓ_partial − (λ/2)rough(β)
    dalpha <- sqrt(sum((alpha - alpha_old)^2)) / max(1, sqrt(sum(alpha_old^2)))
    dbeta  <- if (is.null(beta_old)) Inf else
              sqrt(sum((beta_new - beta_old)^2)) / max(1, sqrt(sum(beta_old^2)))
    dobj   <- if (is.na(obj_old)) Inf else abs(obj_new - obj_old) / max(1, abs(obj_old))
    beta_old <- beta_new; obj_old <- obj_new
    if (verbose) cat(sprintf("iter %2d: |Δα|=%.2e |Δβ|=%.2e |Δpll|=%.2e obj=%.4f\n",
                             it, dalpha, dbeta, dobj, obj_new))
    ## 最良反復（最大の罰則付き部分対数尤度）を保持
    scale.obj <- if (is.finite(best.obj)) max(1, abs(best.obj)) else 1
    if (is.finite(obj_new) && obj_new > best.obj + scale.obj * tol) {
      best.obj <- obj_new; best.alpha <- alpha; best.Theta <- Theta; best.X <- X; best.it <- it; stall <- 0L
    } else stall <- stall + 1L
    ## 発散ハルト（目的関数が最良から大きく低下）：最良反復を採用して停止
    if (is.finite(best.obj) && is.finite(obj_new) && obj_new < best.obj - max(1, abs(best.obj)) * 0.5) {
      conv.msg <- sprintf("halted: penalized loglik diverged at iter %d; returning best iterate (iter %d)", it, best.it)
      converged <- best.it > 0; break
    }
    ## 収束: 3量すべて tol 以下（厳格）, または目的関数がプラトー（patience 反復改善なし）
    if (max(dalpha, dbeta, dobj) <= tol) {
      converged <- TRUE
      conv.msg  <- "converged: max rel. change of (alpha, beta=X*Theta, penalized loglik) <= tol"
      break
    }
    if (stall >= patience) {
      converged <- TRUE
      conv.msg  <- sprintf("converged: penalized loglik plateaued (no improvement > tol for %d iters; best iter %d)",
                           patience, best.it)
      break
    }
  }
  ## 最良反復を採用（発散域や振動の最終点でなく, 最大目的関数の解を返す）。it は停止反復数として保持。
  if (!is.null(best.Theta)) { alpha <- best.alpha; Theta <- best.Theta; X <- best.X }
  best.iter <- best.it
  if (converged && !is.na(opt.conv) && opt.conv != 0)     # 収束したが最終 X-block optim が未収束
    conv.msg <- sprintf("%s; note: X-block optim$convergence=%d at some iter", conv.msg, opt.conv)
  if (!converged) warning("nmf.cox2 did NOT converge: ", conv.msg)

  ## ---- 6. 最終 offset で確定 coxph（robust SE）----
  o_final <- numeric(length(rj))
  S <- Sfun(Theta)
  for (j in seq_len(K)) {
    rows <- rows_by_j[[j]]; ids <- rid[rows]
    o_final[rows] <- as.vector(crossprod(S[, ids, drop = FALSE], X[j, ]))
  }
  o_final <- pmax(pmin(o_final, 30), -30)             # 数値安全（exp 発散回避）
  dat2 <- data.frame(start = rstart, stop = rstop, ev = rev, id = rid, off = o_final)
  Zdf <- as.data.frame(Zexp); colnames(Zdf) <- znames
  dat2 <- cbind(dat2, Zdf)
  fml <- as.formula(paste0("Surv(start,stop,ev) ~ ",
                           paste(znames, collapse = " + "),
                           " + offset(off) + cluster(id)"))
  cox.fit <- coxph(fml, data = dat2, ties = ties)

  ## ---- 6b. β(t) の条件付き推論（X̂ 固定; Satoh 2023 NMF-GCM と同型）----
  ## (α,Θ) ブロック coxph の vcov（クラスタ頑健）からデルタ法: β(t_j)_l=Σ_q X[j,q]Θ_ql が Θ に線形
  ## inference=FALSE（CV 内など）ではこの高コスト計算を丸ごとスキップ（beta.t のみ返す）。
  se.beta.t <- matrix(NA_real_, K, L); beta.vcov <- vector("list", L)
  wald <- NULL; wald.global <- NULL; inf.fit <- NULL
  if (inference) {
  Finf <- matrix(0, length(rj), Q * L); ci <- 1
  for (q in 1:Q) for (l in 1:L) { Finf[, ci] <- X[rj, q] * A[l, rid]; ci <- ci + 1 }  # 列順 (q,l)=(q-1)L+l
  ## 頑健な二次形式: W=θ'Σ^+θ（擬似逆），df=Σ の数値ランク（特異でも破綻しない）
  quad.wald <- function(theta, Vmat) {
    Vs <- (Vmat + t(Vmat)) / 2
    ee <- eigen(Vs, symmetric = TRUE)
    tol <- max(ee$values, 0) * 1e-8 * length(ee$values)
    keep <- ee$values > tol
    if (!any(keep)) return(c(chisq = NA_real_, df = 0))
    U <- ee$vectors[, keep, drop = FALSE]
    W <- sum((crossprod(U, theta)^2) / ee$values[keep])
    c(chisq = as.numeric(W), df = sum(keep))
  }
  ## inf.fit を NMF 解 (α,Θ) から初期化 → 収束を安定化（Finf 列順 (q,l)=(q-1)L+l と一致）
  av <- setNames(rep(0, ncol(Zexp)), colnames(Zexp)); cc <- coef(cox.fit)
  av[names(cc)] <- ifelse(is.na(cc), 0, cc)
  init.inf <- c(as.numeric(av), as.vector(t(Theta)))
  inf.fit <- tryCatch(coxph(Surv(rstart, rstop, rev) ~ Zexp + Finf + cluster(rid),
                            ties = ties, init = init.inf,
                            control = coxph.control(iter.max = 200, eps = 1e-9)),
                      error = function(e) NULL)
  if (!is.null(inf.fit)) {
    V <- vcov(inf.fit); cf <- coef(inf.fit); nz <- ncol(Zexp); idxTh <- (nz + 1):(nz + Q * L)
    ## ---- β(t)=0 の Wald 検定（H0: θ_l=Θ[,l]=0 ⇔ β_l(t)≡0；比例ハザード）----
    ## X̂ を既知としたブロック (α,Θ) の同時推定に基づく条件付き Wald（信頼帯と同じ条件付け）
    wald <- data.frame(covariate = anames, df = Q, chisq = NA_real_, p.value = NA_real_,
                       stringsAsFactors = FALSE)
    for (l in 1:L) {
      posl <- idxTh[(0:(Q - 1)) * L + l]                # 共変量 l の Q 係数の位置 Θ[1:Q,l]
      Vl <- V[posl, posl, drop = FALSE]; beta.vcov[[l]] <- Vl
      qw <- quad.wald(cf[posl], Vl)
      wald$chisq[l] <- qw["chisq"]; wald$df[l] <- qw["df"]
      wald$p.value[l] <- if (qw["df"] > 0) stats::pchisq(qw["chisq"], qw["df"], lower.tail = FALSE) else NA_real_
      for (j in 1:K) se.beta.t[j, l] <- sqrt(max(0, as.numeric(crossprod(X[j, ], Vl %*% X[j, ]))))
    }
    ## 全共変量まとめた大域 Wald（H0: β(t)≡0 for all covariates；全体の比例性）
    qwg <- quad.wald(cf[idxTh], V[idxTh, idxTh, drop = FALSE])
    wald.global <- list(chisq = as.numeric(qwg["chisq"]), df = as.integer(qwg["df"]),
                        p.value = if (qwg["df"] > 0) stats::pchisq(qwg["chisq"], qwg["df"], lower.tail = FALSE) else NA_real_)
  }
  }  # end if(inference)

  beta.t <- X %*% Theta                                # K×L : 時間変動係数 β(t)
  colnames(beta.t) <- anames; colnames(se.beta.t) <- anames
  rownames(beta.t) <- rownames(se.beta.t) <- format(et)
  rownames(X) <- format(et); colnames(X) <- paste0("basis", 1:Q)
  rownames(Theta) <- paste0("basis", 1:Q); colnames(Theta) <- anames

  structure(list(
    alpha = coef(cox.fit)[znames], alpha.raw = alpha,
    X = X, Theta = Theta, beta.t = beta.t, se.beta.t = se.beta.t, beta.vcov = beta.vcov,
    wald = wald, wald.global = wald.global,
    inf.fit = inf.fit, event.times = et,
    cox.fit = cox.fit, offset = o_final, alpha.history = alpha.history,
    rank = Q, X.smooth = X.smooth, iters = it, best.iter = best.iter, best.obj = best.obj,
    converged = converged, conv.message = conv.msg, opt.convergence = opt.conv,
    Z = Z, A = A, A.center = A.center, A.scale = A.scale, call = match.call()
  ), class = "nmf.cox2")
}

# =====================================================================
# nmf.cox2.cf  —  二段階クロスフィット版（定理1が対象とする推定量）
#   各 fold I_k を除いた標本で nuisance（X,Θ）を推定 → ô^{(-k)} を out-of-fold で構成
#   → 全標本の offset付き coxph で α を推定（DML: 標的αは全標本, nuisanceのみ cross-fit）。
# =====================================================================
nmf.cox2.cf <- function(formula, data, A, rank = 2, X.smooth = 0, nfolds = 5,
                        ties = "breslow", seed = NULL, verbose = FALSE, ...) {
  if (!identical(ties, "breslow")) stop("nmf.cox2.cf: only ties='breslow' is supported.")
  mf <- model.frame(formula, data); y <- model.response(mf)
  time <- as.numeric(y[, 1]); status <- as.integer(y[, 2]); n <- length(time)
  Z <- model.matrix(attr(mf, "terms"), mf)
  Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  znames <- colnames(Z)
  if (inherits(A, "formula")) {
    Am <- model.matrix(A, model.frame(A, data))
    Am <- Am[, colnames(Am) != "(Intercept)", drop = FALSE]; Araw <- t(Am)
  } else { Araw <- as.matrix(A); if (ncol(Araw) != n && nrow(Araw) == n) Araw <- t(Araw) }

  et <- sort(unique(time[status == 1])); K <- length(et); s_prev <- c(0, et[-K])
  rid <- rj <- rstart <- rstop <- rev <- integer(0)
  for (j in seq_len(K)) { risk <- which(time >= et[j]); m <- length(risk)
    rid <- c(rid, risk); rj <- c(rj, rep(j, m)); rstart <- c(rstart, rep(s_prev[j], m))
    rstop <- c(rstop, rep(et[j], m)); rev <- c(rev, as.integer(status[risk] == 1 & time[risk] == et[j])) }

  if (!is.null(seed)) set.seed(seed)
  fold <- sample(rep_len(1:nfolds, n))
  Ocross <- matrix(0, n, K)
  for (k in 1:nfolds) {
    tr <- which(fold != k); ho <- which(fold == k)
    fit <- tryCatch(nmf.cox2(formula, data = data[tr, , drop = FALSE],
                  A = Araw[, tr, drop = FALSE], rank = rank, X.smooth = X.smooth,
                  ties = ties, verbose = FALSE, inference = FALSE, ...), error = function(e) NULL)
    if (is.null(fit)) stop(sprintf("nmf.cox2.cf: nuisance fit (trained without fold %d) failed; refusing to continue with a zero offset for that fold.", k))
    Xf <- sapply(seq_len(ncol(fit$X)),
                 function(q) approx(fit$event.times, fit$X[, q], xout = et, rule = 2)$y)  # K×Q
    Aho <- (Araw[, ho, drop = FALSE] - fit$A.center) / fit$A.scale
    Sho <- fit$Theta %*% Aho                          # Q×|ho|
    Ocross[ho, ] <- t(Xf %*% Sho)
  }
  o_row <- pmax(pmin(Ocross[cbind(rid, rj)], 30), -30)
  dat <- data.frame(start = rstart, stop = rstop, ev = rev, id = rid, off = o_row,
                    fold = fold[rid])
  Zdf <- as.data.frame(Z[rid, , drop = FALSE]); colnames(Zdf) <- znames
  dat <- cbind(dat, Zdf)
  ## 定理対応（Lemma hess / (S3)）: 評価 fold 内でリスク集合を作る（strata(fold)）。
  ## → fold k のスコアは fold 内メンバーのみ・全員が同一の out-of-fold nuisance ô^{(-k)} を持ち，
  ##   分母に別 fold の nuisance が混ざらない（旧: 全標本1本の coxph はリスク集合内で混在＝定理と不一致）。
  fml <- as.formula(paste0("Surv(start,stop,ev) ~ ", paste(znames, collapse = " + "),
                           " + offset(off) + strata(fold) + cluster(id)"))
  cox.fit <- coxph(fml, data = dat, ties = ties)
  structure(list(alpha = coef(cox.fit)[znames], cox.fit = cox.fit,
                 offset = Ocross, nfolds = nfolds, rank = rank, X.smooth = X.smooth),
            class = "nmf.cox2.cf")
}

# =====================================================================
# nmf.cox2.cv  —  cross-validated partial likelihood (CVPL; Verweij & van
#   Houwelingen 1993) で平滑化 X.smooth（と任意で rank）を選択する。
#   K-fold; 各 fold で学習側に v2 を当て，held-out の寄与を
#   CVPL_k = ℓ_full(θ^{-k}) − ℓ_train(θ^{-k}) で評価，総和を最大化。
#   ℓ は時間変動オフセット η_i(t)=z_i'α + a_i(std)'β(t) の部分対数尤度（Breslow）。
# =====================================================================
nmf.cox2.cv <- function(formula, data, A, rank = 2,
                        smooth.grid = c(0, 30, 100, 300, 1000, 3000),
                        nfolds = 5, ties = "breslow", seed = NULL,
                        mc.cores = 1, verbose = TRUE, ...) {
  if (!identical(ties, "breslow")) stop("nmf.cox2.cv: only ties='breslow' is supported.")
  mf <- model.frame(formula, data); y <- model.response(mf)
  time <- as.numeric(y[, 1]); status <- as.integer(y[, 2]); n <- length(time)
  Z <- model.matrix(attr(mf, "terms"), mf); Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  if (inherits(A, "formula")) {
    Am <- model.matrix(A, model.frame(A, data)); Am <- Am[, colnames(Am) != "(Intercept)", drop = FALSE]
    Araw <- t(Am)
  } else { Araw <- as.matrix(A); if (ncol(Araw) != n && nrow(Araw) == n) Araw <- t(Araw) }
  if (!is.null(seed)) set.seed(seed)
  fold <- sample(rep_len(1:nfolds, n))

  one <- function(si, k) {                               # (smooth, fold) 1タスク
    sm <- smooth.grid[si]; tr <- which(fold != k)
    fit <- tryCatch(nmf.cox2(formula, data = data[tr, , drop = FALSE],
                  A = Araw[, tr, drop = FALSE], rank = rank, X.smooth = sm,
                  ties = ties, verbose = FALSE, inference = FALSE, ...), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    Astd <- (Araw - fit$A.center) / fit$A.scale          # 学習側統計で標準化
    etk <- fit$event.times; B <- fit$beta.t
    betaFun <- function(tj) vapply(seq_len(ncol(B)),
                  function(l) approx(etk, B[, l], xout = tj, rule = 2)$y, numeric(1))
    al <- setNames(rep(0, ncol(Z)), colnames(Z)); al[names(fit$alpha)] <- fit$alpha
    zA <- as.vector(Z %*% al)
    clip <- function(x) pmax(pmin(x, 30), -30)           # held-out eta のクリップ（発散回避）
    plik <- function(sub) {
      et <- sort(unique(time[sub][status[sub] == 1])); ll <- 0
      for (tj in et) {
        risk <- sub[time[sub] >= tj]; fail <- sub[status[sub] == 1 & time[sub] == tj]
        if (!length(fail) || !length(risk)) next
        bt <- betaFun(tj)
        eta <- clip(zA[risk] + as.vector(crossprod(Astd[, risk, drop = FALSE], bt)))
        m <- max(eta); lse <- m + log(sum(exp(eta - m)))
        ef <- clip(zA[fail] + as.vector(crossprod(Astd[, fail, drop = FALSE], bt)))
        ll <- ll + sum(ef - lse)
      }
      ll
    }
    plik(1:n) - plik(tr)
  }
  tasks <- expand.grid(si = seq_along(smooth.grid), k = 1:nfolds)
  vals <- if (mc.cores > 1)
            unlist(parallel::mclapply(seq_len(nrow(tasks)),
                     function(r) one(tasks$si[r], tasks$k[r]), mc.cores = mc.cores))
          else vapply(seq_len(nrow(tasks)),
                     function(r) one(tasks$si[r], tasks$k[r]), numeric(1))
  nS <- length(smooth.grid)
  vmat <- matrix(vals, nrow = nS)                       # si × fold（CVPL の fold 別寄与; tasks は si 最速）
  cvpl <- rowSums(vmat, na.rm = TRUE)
  ## fold 間 SD から CVPL 総和の SE（1-SE 則用）
  se <- apply(vmat, 1, function(v) { v <- v[is.finite(v)]; if (length(v) < 2) NA_real_ else stats::sd(v) * sqrt(length(v)) })
  best.i <- which.max(cvpl); best <- smooth.grid[best.i]
  ## 1-SE 則: 最良 CVPL から 1SE 以内で最も強い平滑化（大きい X.smooth）を選ぶ（過小平滑化の振動を抑制; #4）
  thr <- cvpl[best.i] - ifelse(is.na(se[best.i]), 0, se[best.i])
  elig <- which(cvpl >= thr); best.1se <- smooth.grid[elig[which.max(smooth.grid[elig])]]
  if (verbose) { for (si in seq_along(smooth.grid))
      cat(sprintf("  X.smooth=%-7g  CVPL=%.3f  SE=%.3f\n", smooth.grid[si], cvpl[si], se[si]))
    cat(sprintf("  best(CVPL argmax)=%g ; best(1-SE, 最も滑らか)=%g\n", best, best.1se)) }
  structure(list(smooth.grid = smooth.grid, cvpl = cvpl, se = se, best.smooth = best,
                 best.smooth.1se = best.1se, rank = rank, nfolds = nfolds), class = "nmf.cox2.cv")
}

# =====================================================================
# nmf.cox2.phtest  —  β(t)=0（＝比例ハザード；a の効果は時間変動しない or ゼロ）の
#   クロスフィット Wald 検定。素朴な in-fit Wald は X̂ をデータから推定した後に
#   同じデータで検定するため後選択バイアスで壊滅的に過剰棄却する（帰無で ~80%）。
#   → 基底 X̂^{(-k)} を fold 外で推定し，各被験者の交互作用特徴量を「その人を含まない
#   基底」で構成（DML/クロスフィット）→ 循環を断ち，χ² 分布が回復する。
#   H0: θ_l=Θ[,l]=0 ⇔ β_l(t)≡0。返り値: 共変量別 Wald（df=数値ランク）＋大域 Wald。
# =====================================================================
nmf.cox2.phtest <- function(formula, data, A, rank = 2, X.smooth = 0,
                            nfolds = 5, ties = "breslow", seed = NULL, verbose = FALSE, ...) {
  if (!identical(ties, "breslow")) stop("nmf.cox2.phtest: only ties='breslow' is supported.")
  Q <- rank
  mf <- model.frame(formula, data); y <- model.response(mf)
  time <- as.numeric(y[, 1]); status <- as.integer(y[, 2]); n <- length(time)
  Z <- model.matrix(attr(mf, "terms"), mf); Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  znames <- colnames(Z)
  if (inherits(A, "formula")) {
    Am <- model.matrix(A, model.frame(A, data)); Am <- Am[, colnames(Am) != "(Intercept)", drop = FALSE]
    Araw <- t(Am)
  } else { Araw <- as.matrix(A); if (ncol(Araw) != n && nrow(Araw) == n) Araw <- t(Araw) }
  L <- nrow(Araw); anames <- rownames(Araw); if (is.null(anames)) anames <- paste0("A", 1:L)
  Ac <- rowMeans(Araw); As <- apply(Araw, 1, sd); As[As < 1e-10] <- 1
  Astd <- (Araw - Ac) / As                                  # 応答非依存の標準化（漏洩なし）
  et <- sort(unique(time[status == 1])); K <- length(et); s_prev <- c(0, et[-K])
  if (!is.null(seed)) set.seed(seed)
  fold <- sample(rep_len(1:nfolds, n))

  ## fold 外で基底 X̂^{(-k)} を推定し，全イベント時刻 et に内挿（K×Q）
  Xcf <- vector("list", nfolds)
  for (k in 1:nfolds) {
    tr <- which(fold != k)
    fit <- tryCatch(nmf.cox2(formula, data = data[tr, , drop = FALSE], A = Araw[, tr, drop = FALSE],
                    rank = Q, X.smooth = X.smooth, ties = ties, verbose = FALSE,
                    inference = FALSE, ...), error = function(e) NULL)
    if (is.null(fit)) stop(sprintf("nmf.cox2.phtest: basis fit (trained without fold %d) failed; refusing to continue with a zero basis for that fold.", k))
    Xcf[[k]] <- sapply(1:ncol(fit$X), function(q) approx(fit$event.times, fit$X[, q], xout = et, rule = 2)$y)
  }

  ## 全標本の counting-process 展開
  rid <- rj <- rstart <- rstop <- rev <- integer(0)
  for (j in 1:K) { risk <- which(time >= et[j]); m <- length(risk)
    rid <- c(rid, risk); rj <- c(rj, rep(j, m)); rstart <- c(rstart, rep(s_prev[j], m))
    rstop <- c(rstop, rep(et[j], m)); rev <- c(rev, as.integer(status[risk] == 1 & time[risk] == et[j])) }
  Zexp <- Z[rid, , drop = FALSE]

  ## クロスフィット Finf: 各行は「その被験者を含まない fold の基底」を使用
  frow <- fold[rid]; basisByRow <- matrix(0, length(rj), Q)
  for (k in 1:nfolds) { rk <- which(frow == k); if (length(rk)) basisByRow[rk, ] <- Xcf[[k]][rj[rk], , drop = FALSE] }
  Finf <- matrix(0, length(rj), Q * L); ci <- 1
  for (q in 1:Q) for (l in 1:L) { Finf[, ci] <- basisByRow[, q] * Astd[l, rid]; ci <- ci + 1 }

  inf.fit <- tryCatch(coxph(Surv(rstart, rstop, rev) ~ Zexp + Finf + cluster(rid), ties = ties,
                            control = coxph.control(iter.max = 200, eps = 1e-9)), error = function(e) NULL)
  if (is.null(inf.fit)) return(NULL)
  V <- vcov(inf.fit); cf <- coef(inf.fit); nz <- ncol(Zexp); idxTh <- (nz + 1):(nz + Q * L)
  quad.wald <- function(theta, Vmat) {
    Vs <- (Vmat + t(Vmat)) / 2; ee <- eigen(Vs, symmetric = TRUE)
    tol <- max(ee$values, 0) * 1e-8 * length(ee$values); keep <- ee$values > tol
    if (!any(keep)) return(c(chisq = NA_real_, df = 0))
    U <- ee$vectors[, keep, drop = FALSE]
    c(chisq = as.numeric(sum((crossprod(U, theta)^2) / ee$values[keep])), df = sum(keep))
  }
  wald <- data.frame(covariate = anames, df = NA_integer_, chisq = NA_real_, p.value = NA_real_,
                     stringsAsFactors = FALSE)
  beta.vcov <- vector("list", L)                          # cross-fit vcov ブロック（信頼帯用）
  for (l in 1:L) {
    posl <- idxTh[(0:(Q - 1)) * L + l]
    Vl <- V[posl, posl, drop = FALSE]; beta.vcov[[l]] <- Vl
    qw <- quad.wald(cf[posl], Vl)
    wald$chisq[l] <- qw["chisq"]; wald$df[l] <- qw["df"]
    wald$p.value[l] <- if (qw["df"] > 0) stats::pchisq(qw["chisq"], qw["df"], lower.tail = FALSE) else NA_real_
  }
  qwg <- quad.wald(cf[idxTh], V[idxTh, idxTh, drop = FALSE])
  wald.global <- list(chisq = as.numeric(qwg["chisq"]), df = as.integer(qwg["df"]),
                      p.value = if (qwg["df"] > 0) stats::pchisq(qwg["chisq"], qwg["df"], lower.tail = FALSE) else NA_real_)
  ## ---- 信頼帯用の cross-fit 版 β(t)（代表基底＝fold 平均，係数＝cross-fit Θ̂，分散＝V^cf）----
  ## 返り値は nmf.cox2.band がそのまま使えるフィールド名（X, beta.t, se.beta.t, beta.vcov, event.times）
  Xbar <- Reduce(`+`, Xcf) / nfolds                        # K×Q 代表基底
  Theta.cf <- matrix(cf[idxTh], nrow = Q, ncol = L, byrow = TRUE)  # 列順 (q,l)=(q-1)L+l と整合
  beta.t <- Xbar %*% Theta.cf                              # K×L
  se.beta.t <- matrix(NA_real_, K, L)
  for (l in 1:L) for (j in 1:K)
    se.beta.t[j, l] <- sqrt(max(0, as.numeric(crossprod(Xbar[j, ], beta.vcov[[l]] %*% Xbar[j, ]))))
  colnames(beta.t) <- colnames(se.beta.t) <- anames
  if (verbose) { cat("Cross-fitted Wald test  H0: beta_l(t)=0\n"); print(wald, row.names = FALSE) }
  structure(list(wald = wald, wald.global = wald.global, nfolds = nfolds,
                 rank = Q, X.smooth = X.smooth,
                 X = Xbar, Theta = Theta.cf, beta.t = beta.t, se.beta.t = se.beta.t,
                 beta.vcov = beta.vcov, event.times = et),
            class = "nmf.cox2.phtest")
}

# =====================================================================
# nmf.cox2.band  —  β(t) の信頼帯（点ごと + 同時）
#   点ごと: β̂(t)±z·se(t)。同時: Θ~N(Θ̂,V_l) をシミュレートし sup_t|Δβ|/se の臨界値
#   （Satoh 2016 の同時信頼区間と同型；X̂ 固定に条件づけた推論）。
# =====================================================================
nmf.cox2.band <- function(fit, covariate = 1, level = 0.95, nsim = 2000, seed = 1) {
  l  <- if (is.character(covariate)) match(covariate, colnames(fit$beta.t)) else covariate
  bt <- fit$beta.t[, l]; se <- fit$se.beta.t[, l]; X <- fit$X; Vl <- fit$beta.vcov[[l]]
  z  <- stats::qnorm(1 - (1 - level) / 2)
  pw <- cbind(lo = bt - z * se, hi = bt + z * se)
  crit <- NA_real_; sim <- pw * NA
  if (!any(is.na(se)) && !is.null(Vl)) {
    if (!is.null(seed)) set.seed(seed)
    ch <- chol(Vl + diag(1e-10, nrow(Vl)))
    supstat <- vapply(seq_len(nsim), function(b) {
      d <- as.vector(X %*% crossprod(ch, stats::rnorm(nrow(Vl))))   # β の摂動 = X (Θ_sim-Θ̂)
      max(abs(d) / pmax(se, 1e-8))
    }, numeric(1))
    crit <- as.numeric(stats::quantile(supstat, level))
    sim <- cbind(lo = bt - crit * se, hi = bt + crit * se)
  }
  list(t = fit$event.times, beta = bt, se = se, pointwise = pw, simultaneous = sim,
       crit = crit, level = level, covariate = colnames(fit$beta.t)[l])
}

print.nmf.cox2 <- function(x, ...) {
  cat("NMF-COX v2 (log-link, X>=0 / Theta signed; partial-likelihood)\n")
  cat(sprintf("rank=%d, X.smooth=%g, iters=%d, converged=%s\n",
              x$rank, x$X.smooth, x$iters, ifelse(isTRUE(x$converged), "TRUE", "FALSE")))
  if (!isTRUE(x$converged)) cat("  ** ", x$conv.message, " **\n", sep = "")
  cat("\nCox coefficients (alpha, net of a's time-varying effect):\n")
  print(round(coef(summary(x$cox.fit)), 4))
  if (!is.null(x$wald)) {
    cat("\nWald test for time-varying effect  H0: beta_l(t)=0  (given basis X):\n")
    w <- x$wald; w$chisq <- round(w$chisq, 3)
    w$p.value <- format.pval(w$p.value, digits = 3, eps = 1e-4)
    print(w, row.names = FALSE)
    if (!is.null(x$wald.global))
      cat(sprintf("Global (all covariates):  chisq=%.3f, df=%d, p=%s\n",
                  x$wald.global$chisq, x$wald.global$df,
                  format.pval(x$wald.global$p.value, digits = 3, eps = 1e-4)))
  }
  invisible(x)
}

# =====================================================================
# phtest.spline  —  非PH検定 H0: delta_l(t)=0（主効果 gamma_l を残す；FIXED 自然スプライン基底）
#   データ駆動 NMF 基底の非識別性・後選択バイアス（in-fit Wald の過剰棄却，
#   cross-fit Xbar の基底混線）を回避するため，時間の事前固定基底
#   B(t)=[1, ns(g(t), df-1)] で beta_l(t)=sum_m theta_lm B_m(t) を張り，
#   coxph の頑健 Wald で共変量ごとに H0: theta_l=0（= beta_l(t)≡0）を検定する。
#   基底がデータに依存しない → クロスフィット不要で素朴 Wald が公称水準に較正。
#   g(t): 時間変換 "identity"/"log"/"rank"（cox.zph と同様）。df>=1（df=1 は定数のみ＝
#   通常の PH 検定，df>=2 で時間変動を許す）。返り値は print.nmf.cox2 と同じ
#   $wald / $wald.global を持つので既存表示・sim と互換。
# =====================================================================
phtest.spline <- function(formula, data, A, df = 3, transform = "identity",
                          ties = "breslow", verbose = FALSE) {
  if (!requireNamespace("splines", quietly = TRUE)) stop("package 'splines' required.")
  mf <- model.frame(formula, data); y <- model.response(mf)
  time <- as.numeric(y[, 1]); status <- as.integer(y[, 2]); n <- length(time)
  Z <- model.matrix(attr(mf, "terms"), mf); Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  znames <- colnames(Z)
  if (inherits(A, "formula")) {
    Am <- model.matrix(A, model.frame(A, data)); Am <- Am[, colnames(Am) != "(Intercept)", drop = FALSE]
    Araw <- t(Am)
  } else { Araw <- as.matrix(A); if (ncol(Araw) != n && nrow(Araw) == n) Araw <- t(Araw) }
  L <- nrow(Araw); anames <- rownames(Araw); if (is.null(anames)) anames <- paste0("A", 1:L)
  Ac <- rowMeans(Araw); As <- apply(Araw, 1, sd); As[As < 1e-10] <- 1
  Astd <- (Araw - Ac) / As                                   # 応答非依存の標準化

  et <- sort(unique(time[status == 1])); K <- length(et); s_prev <- c(0, et[-K])
  gt <- switch(transform, identity = et, log = log(et), rank = rank(et),
               stop("transform must be identity/log/rank"))
  ## 非PH検定: beta_l(t)=gamma_l+delta_l(t)．主効果 gamma_l（定数）は常にモデルへ入れ非検定，
  ## 時間変動 delta_l(t)=ns(g(t),df-1)（中心化，定数を含まない）だけを検定する（H0: delta_l(t)=0）．
  ## 「効果なし」ではなく「比例ハザード」の検定であり，gamma を検定から外すぶん自由度が df-1 に減る．
  if (df < 2) stop("phtest.spline: df>=2 required (tests df-1 time-varying basis functions).")
  Btv <- splines::ns(gt, df = df - 1)                       # K×(df-1) 時間変動基底（定数なし）
  Btv <- sweep(Btv, 2, colMeans(Btv))                       # 中心化→主効果 gamma と分離
  Mtv <- ncol(Btv)                                          # = df-1

  rid <- rj <- rstart <- rstop <- rev <- integer(0)
  for (j in 1:K) { risk <- which(time >= et[j]); m <- length(risk)
    rid <- c(rid, risk); rj <- c(rj, rep(j, m)); rstart <- c(rstart, rep(s_prev[j], m))
    rstop <- c(rstop, rep(et[j], m)); rev <- c(rev, as.integer(status[risk] == 1 & time[risk] == et[j])) }
  Zexp <- Z[rid, , drop = FALSE]
  Amain <- t(Astd[, rid, drop = FALSE]); colnames(Amain) <- paste0("main", seq_len(L))  # 定数主効果 gamma_l（非検定）
  ## 検定対象の交互作用: a_l(std) * Btv_m（列順 (l,m)）
  Ftv <- matrix(0, length(rj), L * Mtv); ci <- 1; colidx <- vector("list", L)
  for (l in 1:L) { idx <- integer(Mtv)
    for (mm in 1:Mtv) { Ftv[, ci] <- Astd[l, rid] * Btv[rj, mm]; idx[mm] <- ci; ci <- ci + 1 }
    colidx[[l]] <- idx }
  fit <- tryCatch(coxph(Surv(rstart, rstop, rev) ~ Zexp + Amain + Ftv + cluster(rid), ties = ties,
                        control = coxph.control(iter.max = 200, eps = 1e-9)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  V <- vcov(fit); cf <- coef(fit); nz <- ncol(Zexp); off <- nz + L    # Zexp(nz)+Amain(L) の後に Ftv
  idxF <- (off + 1):(off + L * Mtv)
  quad.wald <- function(theta, Vmat) {
    Vs <- (Vmat + t(Vmat)) / 2; ee <- eigen(Vs, symmetric = TRUE)
    tol <- max(ee$values, 0) * 1e-8 * length(ee$values); keep <- ee$values > tol
    if (!any(keep)) return(c(chisq = NA_real_, df = 0))
    U <- ee$vectors[, keep, drop = FALSE]
    c(chisq = as.numeric(sum((crossprod(U, theta)^2) / ee$values[keep])), df = sum(keep)) }
  wald <- data.frame(covariate = anames, df = NA_integer_, chisq = NA_real_, p.value = NA_real_,
                     stringsAsFactors = FALSE)
  for (l in 1:L) { pos <- off + colidx[[l]]
    qw <- quad.wald(cf[pos], V[pos, pos, drop = FALSE])
    wald$chisq[l] <- qw["chisq"]; wald$df[l] <- qw["df"]
    wald$p.value[l] <- if (qw["df"] > 0) stats::pchisq(qw["chisq"], qw["df"], lower.tail = FALSE) else NA_real_ }
  qwg <- quad.wald(cf[idxF], V[idxF, idxF, drop = FALSE])
  wald.global <- list(chisq = as.numeric(qwg["chisq"]), df = as.integer(qwg["df"]),
                      p.value = if (qwg["df"] > 0) stats::pchisq(qwg["chisq"], qwg["df"], lower.tail = FALSE) else NA_real_)
  if (verbose) { cat(sprintf("Fixed-spline non-PH test (df=%d, tv=%d, transform=%s)  H0: delta_l(t)=0\n", df, Mtv, transform))
                 print(wald, row.names = FALSE) }
  structure(list(wald = wald, wald.global = wald.global, df = df, transform = transform,
                 basis = Btv, event.times = et, alpha = coef(fit)[znames], cox.fit = fit),
            class = "phtest.spline")
}

# =====================================================================
# nmf.cox2.boot.band  —  β(t) の信頼帯（個体単位 full-refit ブートストラップ; #5）
#   X̂ を固定した条件付き帯は基底推定の不確実性を伝播せず，真の時間変動効果の同時被覆が
#   崩れる（band_coverage_sim 参照）。ここでは各ブートストラップ標本で NMF-COX 全体を
#   再推定（任意で λ_X も CVPL 1-SE 再選択）し，β(t)=XΘ を直接比較して帯を作る。
#   積 β(t) は基底の列置換・符号に不変なので NMF 基底の整列は不要。
#   select=TRUE で各標本の λ_X を CVPL 1-SE で再選択（選択の不確実性も帯に反映）。
#   返り値は nmf.cox2.band と同じ $t/$beta/$pointwise/$simultaneous を持つ。
# =====================================================================
nmf.cox2.boot.band <- function(formula, data, A, rank = 2, X.smooth = 0, smooth.time = "gap",
                               covariate = 1, grid = NULL, B = 400, level = 0.95, ties = "breslow",
                               select = FALSE, smooth.grid = c(100, 300, 1000, 3000),
                               nfolds = 3, seed = 1, mc.cores = 1, maxit = 20, ...) {
  if (!inherits(A, "formula")) stop("nmf.cox2.boot.band: A は formula で渡してください（各標本で再構成するため）。")
  f0 <- nmf.cox2(formula, data = data, A = A, rank = rank, X.smooth = X.smooth,
                 smooth.time = smooth.time, ties = ties, verbose = FALSE, maxit = maxit, ...)
  cn <- colnames(f0$beta.t)
  l  <- if (is.character(covariate)) match(covariate, cn) else covariate
  if (is.null(grid)) grid <- f0$event.times
  b0 <- approx(f0$event.times, f0$beta.t[, l], xout = grid, rule = 2)$y
  n  <- nrow(data); cvar <- cn[l]
  one <- function(b) {
    set.seed(seed + b); idx <- sample.int(n, n, replace = TRUE)
    d <- data[idx, , drop = FALSE]
    sm <- X.smooth
    if (select) { cv <- tryCatch(nmf.cox2.cv(formula, data = d, A = A, rank = rank,
                    smooth.grid = smooth.grid, nfolds = nfolds, seed = 1, verbose = FALSE,
                    smooth.time = smooth.time, ...), error = function(e) NULL)
      if (!is.null(cv)) sm <- cv$best.smooth.1se }
    fb <- tryCatch(nmf.cox2(formula, data = d, A = A, rank = rank, X.smooth = sm,
                    smooth.time = smooth.time, ties = ties, verbose = FALSE, maxit = maxit, ...),
                   error = function(e) NULL)
    if (is.null(fb)) return(rep(NA_real_, length(grid)))
    lb <- match(cvar, colnames(fb$beta.t))                    # 同一共変量（積 β は置換不変）
    approx(fb$event.times, fb$beta.t[, lb], xout = grid, rule = 2)$y
  }
  BT <- if (mc.cores > 1) do.call(rbind, parallel::mclapply(1:B, one, mc.cores = mc.cores))
        else t(vapply(1:B, one, numeric(length(grid))))
  ok <- rowSums(is.na(BT)) == 0; BT <- BT[ok, , drop = FALSE]; Bok <- nrow(BT)
  ## 点ごと: パーセンタイル帯
  pw <- t(apply(BT, 2, function(v) stats::quantile(v, c((1 - level)/2, 1 - (1 - level)/2), na.rm = TRUE)))
  ## 同時 (sup-t): ブートストラップ中心・SD で標準化した sup 統計量の分位点
  bmean <- colMeans(BT); bsd <- apply(BT, 2, stats::sd); bsd <- pmax(bsd, 1e-8)
  supstat <- apply(BT, 1, function(v) max(abs(v - bmean) / bsd))
  crit <- as.numeric(stats::quantile(supstat, level))
  sim <- cbind(lo = b0 - crit * bsd, hi = b0 + crit * bsd)
  list(t = grid, beta = b0, se = bsd, pointwise = cbind(lo = pw[, 1], hi = pw[, 2]),
       simultaneous = sim, crit = crit, level = level, B = Bok, covariate = cvar)
}
