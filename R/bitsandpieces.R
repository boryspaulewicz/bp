#' @export
af = as.factor
#' @export
an = as.numeric
#' @export
ac = as.character
#' @export
st = function(x)ct(x) / sd(x, na.rm = TRUE)
#' @export
ct = function(x)x - mean(x, na.rm = TRUE)

#' Poziom istotności zamieniony na gwiazdki
#' @param x Wektor poziomów istotności
#' @return Wektor tekstowy zawierający 0 lub więcej gwiazdek
#' @export
p2star = function(x)ifelse(x < .001, '***', ifelse(x < .01, '** ', ifelse(x < .05, '*  ', '   ')))
