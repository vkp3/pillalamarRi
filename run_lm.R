#' Runs linear regression to estimate effect of polygenic risk SCORE (PRS) on gene expression
#' @description formula: Expr ~ PRS + Cov, where Cov = covariates
#' 
#' @param expr Expression vector (numeric) with length = #samples
#' @param cov Regression covariates in form [cov x samples]
#' @param SCORE The PRS, not included in `cov`
#' @param method Choose between 'default' or 'two-stage' for lm() method (see desc. in support functions below)
#' @return A [1 x 8] vector output from an lm() like below:
#'         ['intercept', 'beta', 'SE', 't_value', 'pval', 'beta.conf.low', 'beta.conf.high', 'corr.rho']
#'
#' @example using `pbapply::pblapply` to parallelize run_lm() over all genes
#' num.cores = 10
#' lm.res <-
#' pblapply(tx_expr,            # Expression vector list for `pbapply::pblapply`
#'         run_lm,              # This function
#'         cov = cov,           # Covariate matrix, as desribed above
#'         SCORE = SCORE,       # PRS
#'         method = 'two-stage',# Choose between 'default' or 'two-stage'
#'         cl = num.cores)      # Number of cores to parallelize over
#'         
#' lm.res <- simplify2array(lm.res, higher=F)
#' rownames(lm.res) <-
#'  c('intercept',
#'    'beta',
#'    'SE',
#'    't_value',
#'    'pval',
#'    'conf.low',
#'    'conf.high',
#'    'corr.rho')
#' colnames(lm.res) <- get.gene_id_stable(gene.ids)
#'
#' # Sort results by p-value
#' lm.res <- as.data.frame(t(lm.res))
#' lm_res.sort <- lm.res[order(lm.res$pval), ]
#'
#' 
#' Author: Vamsee Pillalamarri
run_lm <- function(expr, cov, SCORE, method='default') {
  res <- switch (method,
    "default" = run_lm_default(expr, cov, SCORE),
    "two-stage" = run_lm_two_stage(expr, cov, SCORE)
  )
  
  return(res)
}

# run_lm() support functions
run_lm_default <- function(expr, cov, SCORE) {
  expr <- as.numeric(expr)
  expr_cov <- cbind(SCORE, expr, cov)
  
  # # Run lm() normal procedure
  lm.fit <- lm(expr ~ ., data = expr_cov)
  lm.fit.summary <- summary(lm.fit)
  # print(lm.fit.summary)

  # Get corr
  cor_expr_score <-
    cor(expr, SCORE)

  # Capture p-val, etc.
  lm.res_ <-
    as.data.frame(t(coef(lm.fit.summary)['SCORE',]))
  intercept <- coef(lm.fit)[1]
  names(intercept) <- NULL
  lm.res_ <-
    cbind(data.frame(intercept),
          lm.res_,
          t(confint(lm.fit)['SCORE', ]),
          cor_expr_score)
  
  return(as.matrix(lm.res_))
}

run_lm_two_stage <- function(expr, cov, SCORE) {
  expr <- as.numeric(expr)
  expr_cov <- cbind(SCORE, expr, cov)
  
  # Run lm with fully-adjusted 2 stage regression residual / adjustment procedure
  # https://aeolister.wordpress.com/2016/07/21/regressing-out-a-covariate-is-problematic/
  # https://stat.ethz.ch/pipermail/r-help/2005-April/068856.html
  # Partly adjusted model issues: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3201714/
  # model 1: lm(Expr ~ Covariates)
  # model 2: lm(PRS ~ Covariates)
  # model 3: either:
  #                 lm(rank normalized residuals of model 1 ~ residuals of model 2) OR
  #                 lm(residuals of model 1 ~ residuals of model 2)
  # This captures the conditional effect of PRS on adjusted gene expression
  lm.fit1 <- lm(expr ~ . - SCORE, data = expr_cov) # regress expr onto cov
  lm.fit2 <- lm(SCORE ~ . - expr, data = expr_cov) # regress prs onto cov
  resid.dat <- data.frame(expr=resid(lm.fit1),
                          SCORE=resid(lm.fit2)) # don't rank normalize adj. gene expression residuals
  # resid.dat <- data.frame(expr=rankNorm(resid(lm.fit1)), 
  #                       SCORE=resid(lm.fit2)) # do rank normalize adj. gene expression residuals

  lm.fit3 <- lm(expr ~ SCORE, data = resid.dat) # regress adjusted variables onto each other
  lm.fit.3.summary <- summary(lm.fit3)

  cor_expr_score <-
    with(resid.dat, cor(expr, SCORE))

  # Capture p-val, etc.
  lm.res_ <-
    as.data.frame(t(coef(lm.fit.3.summary)['SCORE',]))
  intercept <- coef(lm.fit3)[1]
  names(intercept) <- NULL
  lm.res_ <-
    cbind(data.frame(intercept),
          lm.res_,
          t(confint(lm.fit3)['SCORE', ]),
          cor_expr_score)

  return(as.matrix(lm.res_))
}

# OLD Code ----
# # Run lm with partly-adjusted 2-stage "regressing out" procedure (potentially incorrect)
# # Partly-adjusted model issues: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3201714/
# # lm.fit <- lm(expr ~ . - SCORE, data = expr_cov)
# # #lm.fit <- lm(expr ~ ., data = expr_cov)
# # #lm.fit.summary <- summary(lm.fit)
# # lm.fit2 <- lm(residuals(lm.fit) ~ SCORE)
# # lm.fit.summary <- summary(lm.fit2)
#
# # Capture p-val, etc.
# lm.res_ <-
#   as.data.frame(t(coef(lm.fit.summary)['SCORE',]))
# intercept <- coef(lm.fit2)[1]
# names(intercept) <- NULL
# lm.res_ <-
#   cbind(data.frame(intercept),
#         lm.res_,
#         t(confint(lm.fit2)['SCORE', ]),
#         cor_expr_score)
# lm_res[k, ] <- lm.res_
