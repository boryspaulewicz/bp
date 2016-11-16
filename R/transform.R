#' Oblicz wartość funkcji dla podzbiorów i rozwiń do oryginalnej
#' długości.
#' @param dv Zmienna agregowana.
#' @param f Czynnik grupujący.
#' @param fun Funkcja agregująca (ma zwracać skalar).
#' @return Wektor zagregowanych wartości rozwinięty do długości length(dv).
#' @export
aggregate_expand = function(dv, f, fun = mean){
  res = aggregate(dv ~ f, FUN = fun)
  rownames(res) = as.character(res[,1])
  return(res[as.character(f), 2])
}
