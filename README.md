# sentimia

Solving nonlinear equations and application to hydraulics.

R package created **Jemima Nkama**, as part of the introduction to R course, under the supervision of **Mr Yvan Ngneunmeu**.

## Contents

- `bissection()`, `point_fixe()`, `newton()`, `secante()` : four methods to solve f(x) = 0.
- `frottement_colebrook()` : friction factor of a pipe
  (équation de Colebrook-White).
- `perte_de_charge()` : head loss in the pipe (Darcy-Weisbach).
- `reynolds()`, `colebrook_g()`, `colebrook_dg()` : Utility functionc.

## Installation

```r
# from the provided .tar.gz fourni
install.packages("sentimia_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Exemple

```r
library(sentimia)

newton(function(x) x^2 - 2, function(x) 2 * x, 1)

# steel pipe de 100 m, diameter 10 cm, eau à 1 m/s
perte_de_charge(D = 0.10, L = 100, V = 1)
```

## Documentation

- help for each : `?newton`, `?frottement_colebrook`, ...
- User guide : `vignette("introduction", package = "sentimia")`
- Presentation scripts : dossier `inst/scripts/`
# sentimia
