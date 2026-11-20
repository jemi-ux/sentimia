# =====================================================================
# APPLICATION : HYDRAULIQUE (génie civil)
# Reproduit les deux graphiques de la présentation :
#   1) le coefficient de frottement f selon le nombre de Reynolds
#   2) la perte de charge selon le diamètre de la conduite
#
# Chaque point demande de résoudre l'équation de Colebrook-White,
# qui est non linéaire : on utilise uniroot() de R.
# =====================================================================


# ---------------------------------------------------------------------
# ÉTAPE 1 : fonction qui calcule f (équation de Colebrook)
# ---------------------------------------------------------------------
# Re  : nombre de Reynolds
# rel : rugosité relative = eps / D
# On pose x = 1/sqrt(f) et on cherche x tel que g(x) = 0.
colebrook <- function(Re, rel) {
  g <- function(x) x + 2 * log10(rel / 3.7 + 2.51 * x / Re)
  x <- uniroot(g, c(3, 30), tol = 1e-12)$root   # racine entre 3 et 30
  1 / x^2                                        # on revient à f
}

# Test : acier (eps = 0.045 mm), D = 0.10 m, Re = 100000
colebrook(1e5, 0.000045 / 0.10)                  # donne environ 0.02012


# ---------------------------------------------------------------------
# ÉTAPE 2 : f en fonction du nombre de Reynolds
# ---------------------------------------------------------------------
Re_vals <- c(4000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)

f_lisse   <- numeric(length(Re_vals))   # tuyau lisse
f_acier   <- numeric(length(Re_vals))   # acier neuf  (eps/D = 0.00045)
f_rugueux <- numeric(length(Re_vals))   # très rugueux (eps/D = 0.005)

for (i in 1:length(Re_vals)) {
  f_lisse[i]   <- colebrook(Re_vals[i], 0)
  f_acier[i]   <- colebrook(Re_vals[i], 0.000045 / 0.10)
  f_rugueux[i] <- colebrook(Re_vals[i], 0.005)
}

tab_f <- data.frame(Re = Re_vals, lisse = f_lisse,
                    acier = f_acier, rugueux = f_rugueux)
print(tab_f, digits = 4)


# ---------------------------------------------------------------------
# ÉTAPE 3 : perte de charge en fonction du diamètre
# Débit fixé : celui de l'exemple (V = 1 m/s dans la conduite de 0,10 m)
# ---------------------------------------------------------------------
L    <- 100                       # longueur (m)
eps  <- 0.000045                  # rugosité de l'acier (m)
nu   <- 1e-6                      # viscosité cinématique de l'eau (m2/s)
grav <- 9.81                      # pesanteur (m/s2)
Q    <- pi / 4 * 0.10^2 * 1       # débit (m3/s) = 7,85 L/s

diam <- c(0.05, 0.075, 0.10, 0.125, 0.15, 0.20)   # diamètres testés (m)

vitesse <- numeric(length(diam))
hf      <- numeric(length(diam))

for (i in 1:length(diam)) {
  d <- diam[i]
  vitesse[i] <- 4 * Q / (pi * d^2)               # V = Q / section
  Re <- vitesse[i] * d / nu                      # nombre de Reynolds
  f  <- colebrook(Re, eps / d)                   # coefficient de frottement
  hf[i] <- f * (L / d) * vitesse[i]^2 / (2 * grav)   # Darcy-Weisbach
}

tab_hf <- data.frame(Diametre_cm = diam * 100,
                     Vitesse_m_s = round(vitesse, 2),
                     Perte_m     = round(hf, 2))
print(tab_hf)

cat("\nDiamètre divisé par 2 (10 cm -> 5 cm) : perte multipliée par",
    round(hf[1] / hf[3]), "\n")


# ---------------------------------------------------------------------
# ÉTAPE 4 : graphiques
# ---------------------------------------------------------------------
par(mfrow = c(1, 2))

# Graphique 1 : f en fonction de Re (échelle logarithmique en x)
plot(Re_vals, f_rugueux, type = "b", log = "x", col = "orange", pch = 16,
     ylim = c(0, 0.05), xlab = "Nombre de Reynolds Re",
     ylab = "Coefficient de frottement f",
     main = "Le frottement selon l'écoulement")
lines(Re_vals, f_acier, type = "b", col = "darkcyan", pch = 17)
lines(Re_vals, f_lisse, type = "b", col = "navy",     pch = 15)
legend("topright", bty = "n", pch = c(15, 17, 16),
       col = c("navy", "darkcyan", "orange"),
       legend = c("Tuyau lisse", "Acier neuf", "Très rugueux"))

# Graphique 2 : perte de charge selon le diamètre
barplot(hf, names.arg = paste(diam * 100, "cm"), col = "darkcyan",
        xlab = "Diamètre de la conduite",
        ylab = "Perte de charge sur 100 m (m d'eau)",
        main = "Choisir le diamètre")
