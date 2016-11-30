do_files_test = function(dir){
  test = function(f){ res = !isdir(f); res; }
  fun = function(f){ message(f); }
  do_files(test = test, fun = fun, dir);
}
