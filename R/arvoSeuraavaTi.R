arvoSeuraavaTila <- function(muutokset, logml) {
  # Suorittaa yksil�n seuraavan tilan arvonnan

  y <- logml + muutokset # siirron j�lkeiset logml:t
  y <- y - base::max(y)
  y <- exp(y)
  summa <- sum(y)
  y <- y / summa
  y <- cumsum(y)

  i2 <- rand_disc(y) # uusi kori
  suurin <- muutokset(i2)
  return(list(suurin = suurin, i2 = i2))
}
