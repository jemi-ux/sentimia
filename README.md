# sentimia

Résolution d'équations non linéaires et application à l'hydraulique.

Package R réalisé par **Jemima Nkama**, dans le cadre du cours d'introduction
à R, sous la supervision de **Mr Yvan Ngneunmeu**.

## Contenu

- `bissection()`, `point_fixe()`, `newton()`, `secante()` : quatre méthodes
  pour résoudre f(x) = 0.
- `frottement_colebrook()` : coefficient de frottement d'une conduite
  (équation de Colebrook-White).
- `perte_de_charge()` : perte de charge de la conduite (Darcy-Weisbach).
- `reynolds()`, `colebrook_g()`, `colebrook_dg()` : fonctions utilitaires.

## Installation

```r
# depuis le fichier .tar.gz fourni
install.packages("sentimia_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Exemple

```r
library(sentimia)

newton(function(x) x^2 - 2, function(x) 2 * x, 1)

# Conduite d'acier de 100 m, diamètre 10 cm, eau à 1 m/s
perte_de_charge(D = 0.10, L = 100, V = 1)
```

## Documentation

- Aide de chaque fonction : `?newton`, `?frottement_colebrook`, ...
- Guide d'utilisation : `vignette("introduction", package = "sentimia")`
- Scripts de l'exposé : dossier `inst/scripts/`
# sentimia
