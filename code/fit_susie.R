library(susieR)
source("code/utils.R")

# thin wrapper for susie
# puts model in a tibble with alpha, cs, and pip columns
fit.susie <- function(X, y, ...){
  s <- susieR::susie(X, y, ...)
  return(tibble(
    model=list(s),
    alpha=list(s$alpha),
    cs=list(alpha2cs(s$alpha)),
    pip=list(susie_get_pip(s)),
    coef=list(coef.susie(s)[-1])))
}

# SuSiE
fit.susie.auto <- function(X, y, ...){
  s <- susieR::susie_auto(X, y, ...)
  return(tibble(model=list(s), alpha=list(s$alpha), cs=list(alpha2cs(s$alpha)), pip=list(susie_get_pip(s))))
}
