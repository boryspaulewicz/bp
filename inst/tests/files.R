do_files_test = function(dir){
  test = function(f){ res = !is_dir(f); res; }
  fun = function(f){ message(f); }
  do_files(test = test, fun = fun, dir);
}

## Powinien pokazać wszystkie pliki w podkatalogach, bez nazw
## katalogów
do_files_test('../')
