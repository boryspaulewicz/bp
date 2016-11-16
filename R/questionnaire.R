#' @export
quest.sum.compute = function(x, n, rev = NULL, rev.val = NULL,
  means = TRUE){
  if(ncol(x) != n)stop(paste("Argument nie ma", n, 'kolumn'))
  if(is.null(rev)){
    res = x
  }else{
    res = cbind(x[,c(1:n)[-rev]], rev.val-x[,rev])
  }
  neg.cors = sum(cor(res, use = 'pair') < 0) / 2
  if(neg.cors > 0){
    warning(paste(neg.cors, "korelacji miedzy itemami jest ujemnych"))
    ## print(csm(res))
  }
  ## W ten sposób radzimy sobie z brakującymi danymi
  if(means){
    apply(res, 1, function(x)mean(x, na.rm = TRUE) * ncol(res))
  }else{
    rowSums(res, na.rm = TRUE)
  }
}

#' Oblicza wynik surowy dla kwestionariusza CES-D
#' @export
cesd.compute = function(x, ...)quest.sum.compute(x, 20, c(8, 12, 16), 3, ...)

#' Oblicza wynik surowy dla kwestionariusza STAI-X1
#' @export
stai1.compute = function(x, ...)quest.sum.compute(x, 20, c(21, 22, 25, 28, 30, 31, 35, 36, 39, 40) - 20, 5, ...)

#' Oblicza wynik surowy dla kwestionariusza STAI-X2
#' @export
stai2.compute = function(x, ...)quest.sum.compute(x, 20, c(21, 26, 27, 30, 33, 36, 39) - 20, 5, ...)

#' Oblicza wynik surowy dla kwestionariusza KAS
#' @export
kas.compute = function(x, ...)quest.sum.compute(x, 29, c(2, 3, 5, 7, 9, 12, 13, 18, 22, 23, 25, 27, 28), 1, ...)

#' Oblicza wynik surowy dla kwestionariusza Liebovitz'a - skala lęk
#' @export
ll.compute = function(x, ...)quest.sum.compute(x, 24, ...)

#' Oblicza wynik surowy dla kwestionariusza Liebovitz'a - skala unikanie
#' @export
lu.compute = function(x, ...)quest.sum.compute(x, 24, ...)

#' Oblicza wynik surowy dla kwestionariusza SKU
#' @export
sku.compute = function(x, ...)quest.sum.compute(x, 20, c(1, 2, 3, 6, 7, 8, 11, 12, 15, 16, 20), 5, ...)

#' Oblicza wynik surowy dla kwestionariusza SPIN
#' @export
spin.compute = function(x, means = TRUE){
  if(ncol(x) != 17){
    stop('SPIN has 17 items')
  }else{
    if(means){
      apply(x, 1, function(y)mean(y, na.rm = TRUE) * ncol(x))
    }else{ apply(x, 1, sum) }
  }
}

#' Oblicza wynik surowy dla kwestionariusza INTE
#' @export
inte.compute = function(x, ...)quest.sum.compute(x, 33, c(5, 28, 33), 6, ...)

#' Oblicza wynik surowy dla kwestionariusza TAS-26
#' @export
tas26.compute = function(x, ...){
  ddf.pos = c(3, 17, 21)
  ddf.neg = c(6, 9)
  dif.pos = c(4, 8, 10, 14, 20, 22, 23, 25, 26)
  dif.neg = c(1, 12)
  eot.pos = c(7, 19)
  eot.neg = c(11, 13, 24)
  rdd.pos = c(2, 18)
  rdd.neg = c(5, 15, 16)
  ## tas.scales.i <<- list(ddf = c(ddf.pos,ddf.neg),
  ##   dif = c(dif.pos,dif.neg),
  ##   eot = c(eot.pos,eot.neg),
  ##   rdd = c(rdd.pos,rdd.neg))
  ## res = c(ddf.neg, ddf.pos, dif.neg, dif.pos, eot.neg, eot.pos, rdd.neg, rdd.pos)
  ## Ok
  if(ncol(x) != 26){
    stop(paste("TAS ma 26 pozycji, a dostalem", ncol(x)))
  }else{
    res = data.frame(ddf = apply(cbind(x[,ddf.pos], 6-x[,ddf.neg]) ,1, sum),
      dif = apply(cbind(x[,dif.pos], 6-x[,dif.neg]) ,1, sum),
      eot = apply(cbind(x[,eot.pos], 6-x[,eot.neg]) ,1, sum),
      rdd = apply(cbind(x[,rdd.pos], 6-x[,rdd.neg]) ,1, sum))
    res$tas = rowSums(res)
  }
  res
}

#' Oblicza wynik surowy dla kwestionariusza POMS
#' @export
poms.compute = function(x, ...){
  anxiety.pos = c(2, 10, 16, 20, 26, 27, 34, 41)
  anxiety.neg = c(22)
  anger = c(3, 12, 17, 24, 31, 33, 39, 42, 47, 52, 53, 57)
  fatigue = c(4, 11, 29, 40, 46, 49, 65)
  depression = c(5, 9, 14, 18, 21, 23, 32, 35, 36, 44, 45, 48, 58, 61, 62)
  vigor = c(7, 15, 19, 38, 51, 56, 60, 63)
  confusion.pos = c(8, 28, 37, 50, 59, 64)
  confusion.neg = c(54)
  friendliness = c(1, 6, 13, 25, 30, 43, 55)
  ## res = c(anxiety.pos, anxiety.neg, anger, fatigue, depression, vigor, confusion.pos, confusion.neg, friendliness)
  ## plot(jitter(sort(res)))
  ##
  ## Wszystkie itemy są sklasyfikowane dokladnie raz
  if(ncol(x) != 65){
    stop(paste("POMS ma 65 pozycji, a dostalem", ncol(x)))
  }else{
    data.frame(anxiety = apply(cbind(x[,anxiety.pos], 5-x[,anxiety.neg]) ,1, sum),
               anger = apply(x[,anger] ,1, sum),
               fatigue = apply(x[,fatigue] ,1, sum),
               depression = apply(x[,depression] ,1, sum),
               vigor = apply(x[,vigor] ,1, sum),
               confusion = apply(cbind(x[,confusion.pos], 5-x[,confusion.neg]) ,1, sum),
               friendliness = apply(x[,friendliness] ,1, sum))
  }
}

#' Oblicza wynik surowy dla kwestionariusza PANAS-C 20
#' @export
panas20.compute = function(x, means = TRUE, ...){
  panas.pos = c(1, 3, 5, 6, 8, 9, 14, 16, 19, 20)
  if(ncol(x) != 20){
    stop(paste("PANAS C20 ma 20 pozycji, a dostalem", ncol(x)))
  }else{
    if(means){
      list(pa = apply(x[,panas.pos], 1, function(y)mean(y, na.rm = TRUE) * ncol(x)),
           na = apply(x[,-panas.pos], 1, function(y)mean(y, na.rm = TRUE) * ncol(x)))
    }else{
      list(pa = rowSums(x[,panas.pos], ...), na = rowSums(x[,-panas.pos], ...))
    }
  }
}

#' Oblicza wynik surowy dla kwestionariusza PANAS-S 20
#' @export
panas20s.compute = function(x, means = TRUE, ...){
  panas.pos = c(1, 3, 5, 6, 8, 9, 10, 16, 19, 20)
  if(ncol(x) != 20){
    stop(paste("PANAS C20 ma 20 pozycji, a dostalem", ncol(x)))
  }else{
    if(means){
      list(pa = apply(x[,panas.pos], 1, function(y)mean(y, na.rm = TRUE) * ncol(x)),
           na = apply(x[,-panas.pos], 1, function(y)mean(y, na.rm = TRUE) * ncol(x)))
    }else{
      list(pa = rowSums(x[,panas.pos], ...), na = rowSums(x[,-panas.pos], ...))
    }
  }
}

#' Oblicza wynik surowy dla kwestionariusza KRE
#' @export
kre.compute = function(x, ...){
  rea.pos = c(1, 3, 5, 7, 8, 10)
  if(ncol(x) != 10){
    stop(paste("KRE ma 10 pozycji, a dostalem", ncol(x)))
  }else{
    list(rea = rowSums(x[,rea.pos], ...), sup = rowSums(x[,-rea.pos], ...))
  }
}

#' Oblicza wynik surowy dla kwestionariusza FNE
#' @export
fne.compute = function(x, ...)quest.sum.compute(x, 30, c(1, 4, 6, 8, 10, 12, 15, 16, 18, 21, 23, 26, 27), 1, ...)
