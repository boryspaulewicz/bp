% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat.R
\name{contrast_matrix}
\alias{contrast_matrix}
\title{Oblicza macierz kontrastów dla tabelki (summary) modelu}
\usage{
contrast_matrix(m, rnd = 2, draw_stars = TRUE, fun = I)
}
\description{
Zakładamy, że druga kolumna ramki (m) zawiera błędy standardowe i
dla każdego porównania współczynników stosujemy średnią
geometryczną ich błędów standardowych jako błąd standardowy
różnicy.
}

