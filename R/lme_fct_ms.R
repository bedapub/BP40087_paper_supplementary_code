# helper function to call lmer and glht on all parameters
# code modified from Nicolas Staedler
wrap_lme_ms <- function(x, Kstar) {
  t.contr <- trimws(gsub("\\(.*", "", rownames(Kstar)))
  t.contr <- factor(t.contr, levels = unique(t.contr))

  lme_fit <- lme4::lmer(cfb ~ VISITDY + (1 | USUBJID), data = x)
  Kstar.l <- split.data.frame(Kstar, t.contr)

  t.res <- lapply(Kstar.l, function(x) {
    lhl <- multcomp::glht(lme_fit, linfct = x)

    res.contr <- data.frame(
      cbind(
        confint(lhl, calpha = multcomp::univariate_calpha())$confint,
        summary(lhl, test = multcomp::adjusted("none"))$test$pvalues,
        rep(summary(lhl, test = multcomp::Chisqtest())$test$pvalue, nrow(x))
      )
    )
    colnames(res.contr) <- c("estimate", "lower.limit", "upper.limit", "pval", "global.pval")
    return(res.contr)
  })
  res <- do.call(rbind, t.res)
  rownames(res) <- gsub(".*\\.", "", rownames(res))

  res <- data.frame(res) %>%
    mutate(
      CONTR = rownames(.),
      VISITDY = as.numeric(rownames(.))
    )
  return(res)
}

lme_fct_ms <- function(mydf) {

  # linear mixed effects modeling
  # set-up contrast
  tmp <- expand.grid(VISITDY = levels(mydf$VISITDY))
  lh <- model.matrix(~ VISITDY, data = tmp)
  K <- diag(rep(1,nlevels(mydf$VISITDY)))
  Kstar <- K %*% lh
  rownames(Kstar) <- levels(mydf$VISITDY)

  # run mixed effects analysis and collect results
  res.contr.lme <- mydf %>%
    group_by(LABORESU, LBTSTDTL, LBMTDTL) %>%
    do(wrap_lme_ms(., Kstar = Kstar)) %>%
    ungroup() %>%
    select(LABORESU, LBTSTDTL, LBMTDTL, CONTR, VISITDY, estimate, lower.limit, upper.limit, pval, global.pval) %>%
    arrange(LABORESU, LBMTDTL, LBTSTDTL, CONTR, VISITDY)

  # adjustment for multiple comparison
  # effect estimate on odds-ratio and fold change scale
  res.contr.lme.adj <- res.contr.lme %>%
    #group_by(CONTR) %>%
    mutate(pval.adj = p.adjust(pval, method = "fdr"), g.pval.adj = p.adjust(global.pval, method = "fdr")) %>%
    ungroup() %>%
    mutate(
      Effect = ifelse(LABORESU == "%", exp(estimate),
        ifelse(LABORESU == "cells/uL", 2^{
          estimate
        }, estimate)
      ),
      LO = ifelse(LABORESU == "%", exp(lower.limit),
        ifelse(LABORESU == "cells/uL", 2^{
          lower.limit
        }, lower.limit)
      ),
      UP = ifelse(LABORESU == "%", exp(upper.limit),
        ifelse(LABORESU == "cells/uL", 2^{
          upper.limit
        }, upper.limit)
      )
    ) %>%
    arrange(LABORESU, CONTR, pval.adj, desc(abs(Effect))) %>%
    mutate(VISITDY = factor(VISITDY))

  return(res.contr.lme.adj)
}
