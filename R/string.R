#' Zwraca wektor logiczny odpowiadający literałom pasującym do wzorca
#' @export
rmatch = function(pattern, vector){
  res = TRUE
  for(i in 1:length(vector)){
      if(length(grep(pattern, vector[i])) > 0){
          res[i] = TRUE
      }else{
          res[i] = FALSE
      }
  }
  res
}
