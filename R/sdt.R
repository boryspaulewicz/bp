## TODO: NA tam, gdzie nie było odpowiednich ratingów

#' Oblicza proporcje trafień i fałszywych alarmów dla poszczególnych ocen
#' @param stim Wektor binarny reprezentujący bodziec.
#' @param dec Wektor binarny reprezentujący decyzję.
#' @param rating Wektor reprezentujący oceny na skali 1-max.rating.
#' @param max.rating Maksymalna wartość oceny. Domyślnie max(rating).
#' @param f Czynnik, który ma oddzielać podzbiory danych do obliczania punktów ROC.
#' @return Ramka z kolumnami f, resp (przekodowane oceny), hi i fa.
#' @export
roc_points = function(stim, dec, rating, max.rating = max(rating), f = 1){
    resp = NA
    resp[dec == 0] = max.rating - rating[dec == 0] + 1
    resp[dec == 1] = rating[dec == 1] + max.rating
    res = expand.grid(f = unique(f), resp = 1:(max.rating * 2 - 1), hi = NA, fa = NA)
    for(g in unique(f)){
        hi = fa = NA
        for(i in 1:(max.rating * 2 - 1)){
            ss = (stim == 1) & (g == f)
            hi[i] = mean(resp[ss] > i)
            ss = (stim == 0) & (g == f)
            fa[i] = mean(resp[ss] > i)
        }
        res[res$f == g, 'hi'] = hi
        res[res$f == g, 'fa'] = fa
    }
    res
}
