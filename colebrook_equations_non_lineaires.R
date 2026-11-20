# =====================================================================
# EXPOSÉ : ÉQUATIONS NON LINÉAIRES  (cours d'introduction à R)
# Application : génie civil / hydraulique
#
# Problème : calculer le coefficient de frottement f d'une conduite
# d'eau avec l'équation de Colebrook-White :
#
#   1/sqrt(f) = -2 * log10( eps/(3.7*D) + 2.51/(Re*sqrt(f)) )
#
# f est des deux côtés : on ne peut pas l'isoler, il faut une méthode
# numérique. On pose x = 1/sqrt(f) et on cherche x tel que g(x) = 0.
# =====================================================================


# ---------------------------------------------------------------------
# ÉTAPE 1 : les données (variables)
# ---------------------------------------------------------------------
D   <- 0.10      # diamètre de la conduite (m)
eps <- 0.000045  # rugosité de l'acier (m)
Re  <- 100000    # nombre de Reynolds (écoulement turbulent)
tol <- 1e-10     # précision demandée


# ---------------------------------------------------------------------
# ÉTAPE 2 : la fonction g(x) et sa dérivée dg(x)
# ---------------------------------------------------------------------
# g(x) = 0  <=>  équation de Colebrook avec x = 1/sqrt(f)
g <- function(x) {
  x + 2 * log10(eps / (3.7 * D) + 2.51 * x / Re)
}

# dg(x) : dérivée de g, nécessaire pour la méthode de Newton
dg <- function(x) {
  1 + (2 / log(10)) * (2.51 / Re) / (eps / (3.7 * D) + 2.51 * x / Re)
}

# phi(x) : on réécrit g(x) = 0 sous la forme x = phi(x) (point fixe)
phi <- function(x) {
  -2 * log10(eps / (3.7 * D) + 2.51 * x / Re)
}


# ---------------------------------------------------------------------
# ÉTAPE 3 : les quatre méthodes, écrites comme des fonctions R
# Chaque fonction renvoie la racine et le nombre d'itérations.
# ---------------------------------------------------------------------

# --- Méthode 1 : bissection -------------------------------------------
# On coupe l'intervalle [a, b] en deux jusqu'à atteindre la précision.
bissection <- function(f, a, b, tol) {
  k <- 0                              # compteur d'itérations
  repeat {
    k <- k + 1
    m <- (a + b) / 2                  # milieu de l'intervalle
    if (abs(f(m)) < tol || (b - a) / 2 < tol) break   # précision atteinte
    if (f(a) * f(m) < 0) {            # la racine est dans [a, m]
      b <- m
    } else {                          # la racine est dans [m, b]
      a <- m
    }
  }
  c(racine = m, iterations = k)
}

# --- Méthode 2 : point fixe -------------------------------------------
# On répète x <- phi(x) jusqu'à ce que x ne change presque plus.
point_fixe <- function(phi, x0, tol) {
  x <- x0
  k <- 0
  repeat {
    k <- k + 1
    x_new <- phi(x)
    if (abs(x_new - x) < tol) break
    x <- x_new
  }
  c(racine = x_new, iterations = k)
}

# --- Méthode 3 : Newton-Raphson ---------------------------------------
# x_new = x - f(x) / f'(x)   (on suit la tangente de la courbe)
newton <- function(f, df, x0, tol) {
  x <- x0
  k <- 0
  repeat {
    k <- k + 1
    x_new <- x - f(x) / df(x)
    if (abs(x_new - x) < tol) break
    x <- x_new
  }
  c(racine = x_new, iterations = k)
}

# --- Méthode 4 : sécante ----------------------------------------------
# Comme Newton, mais la dérivée est remplacée par une pente entre
# deux points : pas besoin de calculer f'.
secante <- function(f, x0, x1, tol) {
  k <- 0
  repeat {
    k <- k + 1
    x2 <- x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    if (abs(x2 - x1) < tol) break
    x0 <- x1
    x1 <- x2
  }
  c(racine = x2, iterations = k)
}


# ---------------------------------------------------------------------
# ÉTAPE 4 : on lance les quatre méthodes
# ---------------------------------------------------------------------
r1 <- bissection(g, 3, 12, tol)       # on cherche entre x = 3 et x = 12
r2 <- point_fixe(phi, 8, tol)         # point de départ x0 = 8
r3 <- newton(g, dg, 8, tol)           # point de départ x0 = 8
r4 <- secante(g, 6, 10, tol)          # deux points de départ

# Solution de référence donnée directement par R
x_ref <- uniroot(g, c(3, 12), tol = 1e-12)$root


# ---------------------------------------------------------------------
# ÉTAPE 5 : tableau de comparaison (data.frame)
# ---------------------------------------------------------------------
x_sol   <- unname(c(r1["racine"], r2["racine"], r3["racine"], r4["racine"]))
nb_iter <- unname(c(r1["iterations"], r2["iterations"],
                    r3["iterations"], r4["iterations"]))

resultats <- data.frame(
  Methode    = c("Bissection", "Point fixe", "Newton", "Secante"),
  x          = x_sol,
  f          = 1 / x_sol^2,             # on revient à f = 1/x^2
  Iterations = nb_iter
)

print(resultats, digits = 8)
cat("\nSolution de référence (uniroot) : x =", x_ref,
    "  f =", 1 / x_ref^2, "\n")


# ---------------------------------------------------------------------
# ÉTAPE 6 : vérification de l'équation de Colebrook
# ---------------------------------------------------------------------
f_sol  <- 1 / x_sol[3]^2                # coefficient de frottement (Newton)
gauche <- 1 / sqrt(f_sol)
droite <- -2 * log10(eps / (3.7 * D) + 2.51 / (Re * sqrt(f_sol)))
cat("\nVérification : 1/sqrt(f) =", gauche,
    "  second membre =", droite, "\n")
cat("Écart =", abs(gauche - droite), "(proche de 0 : l'équation est satisfaite)\n")


# ---------------------------------------------------------------------
# ÉTAPE 7 : application, la perte de charge (Darcy-Weisbach)
# ---------------------------------------------------------------------
L   <- 100                  # longueur de la conduite (m)
V   <- 1                    # vitesse de l'eau (m/s)
grav <- 9.81                # pesanteur (m/s2)

hf <- f_sol * (L / D) * V^2 / (2 * grav)
cat("\nPerte de charge sur", L, "m : hf =", round(hf, 2), "m d'eau\n")


# ---------------------------------------------------------------------
# ÉTAPE 8 : graphiques
# ---------------------------------------------------------------------
par(mfrow = c(1, 2))        # deux graphiques côte à côte

# Graphique 1 : la courbe g(x). La racine est là où elle coupe l'axe x.
curve(g(x), from = 3, to = 12, col = "blue", lwd = 2,
      xlab = "x = 1/sqrt(f)", ylab = "g(x)",
      main = "Équation de Colebrook : g(x) = 0")
abline(h = 0, lty = 2)                       # axe horizontal en pointillés
points(x_ref, 0, pch = 19, col = "red", cex = 1.5)

# Graphique 2 : nombre d'itérations de chaque méthode
barplot(nb_iter, names.arg = resultats$Methode,
        col = c("orange", "darkgreen", "red", "purple"),
        ylab = "Nombre d'itérations",
        main = "Quelle méthode converge le plus vite ?")
