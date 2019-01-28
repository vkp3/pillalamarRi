#' Permutation analysis for tissue
#' 
#' @param permutations Number of permutations
#' @param tx_expr Expression matrix in form: [genes x samples]. Will be converted to a list(!) of gene-expr vectors (if not input as list)
#' @param gene.ids Character vector of gene IDs, corresponding to rows in `tx_expr` [genes x samples]
#' @param cov Regression covariates [cov x samples]
#' @param SCORE Main covariate to be permuted (not included in `cov`)
#' 
#' @export
#' 
#' @return A 2-element list: [[1]]=95% Cutoff, [[2]]=[pvals x permutations] (pvals = all pvals for each gene / permutation)
#' 
#' @examples 
#' res <- perm.for.cutoff2(permutations=2, tx_expr, gene.ids, cov, SCORE, num.cores=10)
#' 
#' Author: Vamsee Pillalamarri
#' Adapted from Stephanie Yang's code: https://github.com/syyang93/analyzeR/blob/master/R/perm.for.cutoff.R

perm.for.cutoff2 <- function(permutations = 100, tx_expr, gene.ids, cov, SCORE, num.cores = 10) {
  suppressPackageStartupMessages(library(parallel))
  suppressPackageStartupMessages(library(stringr))
  
  # Support functions
  get.gene_id_stable <- 
    function(x) {
      return(str_extract(string = x, pattern = '(ENSG)[0-9]+'))
    }
  # source('/dcs01/arking/arkinglab/active/projects/npd.PRS.GTEx/code/run_lm.R') # run_LM()
  
  # If `tx_expr` is not a list,
  # then transform it (assumes tx_expr is a [genes x smpls] expression matrix) to list
  if(!any(class(tx_expr)=="list")) {
    tx_expr <-
      as.list(as.data.frame(t(as.matrix(tx_expr))))
  }
  
  print(paste0('Permutations running on: ', length(tx_expr),' genes.'))
  
  perm.res <- rep(NA, permutations) # min p-values by permutation
  all.pvals <- matrix(NA, nrow=length(tx_expr), 
                      ncol=permutations, dimnames = list(get.gene_id_stable(gene.ids)))
  
  for(i in 1:permutations) {
    print(paste0('Permutation #: ', i, sep=''))
    SCORE.perm <- sample(SCORE)
    
    # Use `parallel::mclapply` to parallelize run_lm over gene expr vectors, using `num.cores` threads
    # To use `mclapply``, make sure to convert expr vectors to a list, not matrix (requirement of `lapply`)
    # Each gene's expr vector is one element of the list `tx_expr`.
    # Note: `run_lm` must return [1 x 8] vector output from an lm() like below:
    # ['intercept', 'beta', 'SE', 't_value', 'pval', 'beta.conf.low', 'beta.conf.high', 'corr.rho']
    lm.res <-
      mclapply(tx_expr,
               run_lm,
               cov = cov,
               SCORE = SCORE.perm,
               mc.cores = num.cores)
    lm.res <- simplify2array(lm.res, higher=F)
    rownames(lm.res) <-
      c('intercept',
        'beta',
        'SE',
        't_value',
        'pval',
        'conf.low',
        'conf.high',
        'corr.rho')
    colnames(lm.res) <- get.gene_id_stable(gene.ids)
    
    # Harvest results
    all.pvals[,i] = lm.res['pval',]
    perm.res[i] <- min(lm.res['pval',])
    print(paste0('min p-val: ', perm.res[i], sep=''))
  }
  
  # Gather 95% cutoff
  cutoff_95pct <- sort(perm.res)[ceiling(0.05 * permutations)]
  print(paste0('95% Cutoff: ', formatC(cutoff_95pct, format = 'e', digits = 3), sep=''))
  
  return(list(cutoff_95pct, all.pvals))
}
