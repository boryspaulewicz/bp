#' Test katalogowości pliku (string)
#' @export
is.dir = function(f){
  if(!is.na(file.info(f)$isdir) && file.info(f)$isdir){
    res = TRUE;
  }else{ res = FALSE; }
  res;
}

#' Przeszukaj (pod)katalogi i sklej ramki wczytane z plików pasujących
#' do podanego wzorca.
#' @export
df.from.files = function(test.fun, verbose = TRUE, header = FALSE, sep = '\t',
  fun = function(fname){
    res = read.table(fname, header = header, sep = sep, ...)
    res
  }, ...){
  do.files(test.fun, function(x, result){
    if(verbose)print(getwd())
    df = fun(x)
    ## Dodajemy nazwe pliku, czasem tam jest zakodowany warunek albo inne informacje
    df = cbind(df, paste(getwd(), x, sep = '/'))
    if(is.logical(result)){
      result = df
    }else{
      result = rbind(result, df)
    }
    result
  }, result = TRUE)
}

#' Rekurencyjne operowanie na plikach/scierzkach.
#'
#' Dla kazdego pliku, jesli spelnia zadany test to wykonaj na nim
#' zadana operacje. Jesli jest katalogiem to wejdz i rekurencyjnie rób
#' dalej swoje. dir to katalog początkowy.
#' @export
do.files = function(test = function(x) TRUE, fun = function(x, result) message(x) , dir = ".", result = NULL){
  if(dir != ".")setwd(dir);
  res = list.files();
  if(length(res) > 0){
    for(f in res){
      ## Jesli to zwykly plik to test i ewentualnie funkcja
      if(!is.dir(f)){
        if(test(f)){
          if(!is.null(result)){
            result = fun(f, result = result)
          }else{
            fun(f)
          }
        }
      }else{
      ##Jesli to jest katalog to do opracuj
        if(!is.null(result)){
          result = do.files(test, fun, f, result = result)
        }else{
          result = do.files(test, fun, f)
        }
        setwd("..")
      }
    }
    result
  }
}

#' Test funkcji \code{do.files} - wypisanie wszystkich plikóww w
#' podkatalogach.
#' @export
do.files.test = function(dir){
  test = function(f){ res = !isdir(f); res; }
  fun = function(f){ message(f); }
  do.files(test = test, fun = fun, dir);
}
