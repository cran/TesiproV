# ======================================================================
# File: tests/testthat/test_rng_reproducibility.R
# ======================================================================
# library(testthat)
# library(TesiproV)   # das zu testende Paket

test_that("gleiche RNG-Seeds -> identische MC-Ergebnisse (reproduzierbar)", {

  ## ---------------------------------------------------------------
  ## 1) Plattform‑abhängige Thread‑Anzahl (Windows → 1, sonst 2)
  ## ---------------------------------------------------------------
  if (.Platform$OS.type == "windows") {
    future::plan(future::multisession)
    use_threads <- 2L
  } else {
    use_threads <- 2L # klein, damit der Test schnell läuft
  }

  ## ---------------------------------------------------------------
  ## 2) Basis‑Variablen (ohne DistributionPackage)
  ## ---------------------------------------------------------------
  var1 <- PROB_BASEVAR(
    Name             = "E",
    Description      = "Effect",
    DistributionType = "norm",
    Mean             = 10,
    Sd               = 0.2
  )
  var2 <- PROB_BASEVAR(
    Name             = "R",
    Description      = "Resistance",
    DistributionType = "norm",
    Mean             = 35,
    Sd               = 0.5
  )

  ## ---------------------------------------------------------------
  ## 3) Limit‑State‑Funktion
  ## ---------------------------------------------------------------
  lsf <- SYS_LSF(vars = list(var1, var2), name = "LSF XY")
  lsf$func <- function(E, R) {
    R - E^1.5
  }
  lsf$check() # stellt sicher, dass alles konsistent ist

  ## ---------------------------------------------------------------
  ## 4) Hilfsfunktion: erstellt ein SYS_PROB‑Objekt, führt die
  ##    Maschinen aus und gibt die Ergebnis‑Liste zurück.
  ## ---------------------------------------------------------------
  run_one <- function(seed_val) {
    ## – FORM – (nur für den Vergleich, kein Einfluss auf RNG)
    machine_form <- PROB_MACHINE(
      name = "FORM Rack.-Fieß.",
      fCall = "FORM",
      options = list(
        n_optim = 20,
        loctol = 0.001,
        optim_type = "rackfies"
      )
    )

    ## – MC‑CRUDE – das ist die Maschine, deren RNG wir testen
    machine_mc <- PROB_MACHINE(
      name = "MC CoV 0.05",
      fCall = "MC_CRUDE",
      options = list(
        n_max = 5e5, # genug, damit die Simulation fertig wird
        cov_user = 0.05,
        use_threads = use_threads,
        seed = seed_val # <-- entscheidender Seed
      )
    )

    ## – SORM – (nur für Vollständigkeit)
    machine_sorm <- PROB_MACHINE(
      name   = "SORM",
      fCall  = "SORM"
    )

    ## Gesamtsystem zusammenbauen
    prob_system <- SYS_PROB(
      sys_input    = list(lsf),
      probMachines = list(machine_form, machine_mc, machine_sorm),
      debug.level  = 0
    )

    ## Maschinen ausführen → Ergebnis‑Liste (benannt nach den Maschinen)
    prob_system$runMachines()
  }

  ## ---------------------------------------------------------------
  ## 5) Zwei unabhängige Läufe mit exakt demselben Seed
  ## ---------------------------------------------------------------
  set.seed(99) # Master‑Seed – hat keinen Einfluss,
  # weil jede Maschine ihren eigenen Seed nutzt
  res_a <- run_one(seed_val = 1234) # erster Durchlauf
  res_b <- run_one(seed_val = 1234) # zweiter Durchlauf (gleicher Seed)

  ## ---------------------------------------------------------------
  ## 6) Ergebnis‑Auswertung
  ## ---------------------------------------------------------------
  # Wir prüfen nur die MC‑Maschine (der Rest ist nicht von RNG‑Repro-
  # duzierbarkeit betroffen).  Das Ergebnis‑Objekt hat das Element $pf.
  pf_a <- res_a[["MC CoV 0.05"]][["pf"]]
  pf_b <- res_b[["MC CoV 0.05"]][["pf"]]

  # Beide Schätzungen müssen exakt übereinstimmen (bis auf Rundungs‑Fehler
  # des R‑Doubles).  Ein sehr kleine Toleranz von 1e‑14 ist sicher.
  expect_equal(pf_a, pf_b,
    tolerance = 1e-14,
    info = "Mit identischem Seed und gleicher Thread‑Zahl muss das MC‑Ergebnis reproduzierbar sein"
  )
})


# ------------------------------------------------------------------
# 0) Hilfsfunktion – Minimal‑Beispiel (ohne DistributionPackage)
# ------------------------------------------------------------------
make_demo_objects <- function() {
  var1 <- PROB_BASEVAR(
    Name             = "E",
    Description      = "Effect",
    DistributionType = "norm",
    Mean             = 10,
    Sd               = 0.2
  )
  var2 <- PROB_BASEVAR(
    Name             = "R",
    Description      = "Resistance",
    DistributionType = "norm",
    Mean             = 35,
    Sd               = 0.5
  )

  lsf_obj <- SYS_LSF(vars = list(var1, var2), name = "LSF XY")
  lsf_obj$func <- function(E, R) {
    R - E^1.5
  }
  lsf_obj$check()

  lsf_fun <- function(x) lsf_obj$func(x[1], x[2])
  lDistr_list <- list(list(var1), list(var2))
  dp <- c(var1$Mean, var2$Mean)

  list(
    lsf_fun = lsf_fun,
    lDistr = lDistr_list,
    dp = dp,
    lsf_obj = lsf_obj
  )
}

# ------------------------------------------------------------------
# 2) Ungültiger densityType → erklärender Fehler
# ------------------------------------------------------------------
test_that("invalid densityType erzeugt einen erklärenden Fehler", {
  demo <- make_demo_objects()
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

  expect_error(
    MC_IS(
      lsf = demo$lsf_fun,
      lDistr = demo$lDistr,
      dps = demo$dp,
      cov_user = 0.10,
      n_batch = 1000,
      n_max = 1e5,
      use_threads = 2,
      backend = backend,
      densityType = "not_a_distribution", # <‑ falscher Typ
      seed = 111,
      dataRecord = FALSE
    ),
    regexp = "unsupported|unknown|invalid|densityType",
    info = "Ein unbekannter densityType muss mit einer klaren Fehlermeldung abbrechen"
  )
})


# ------------------------------------------------------------------
# 3) Edge‑Case n_batch = 0 -> sinnvolle Fehlermeldung
# ------------------------------------------------------------------
test_that("n_batch = 0 führt zu einem Fehler", {
  demo <- make_demo_objects()
  if (.Platform$OS.type == "windows") {
    backend <- "future"
  } else {
    backend <- "parallel"
  }

  expect_error(
    MC_IS(
      lsf = demo$lsf_fun,
      lDistr = demo$lDistr,
      dps = demo$dp,
      cov_user = 0.10,
      n_batch = 0,
      n_max = 1e5,
      use_threads = 2,
      backend = backend,
      seed = 111,
      dataRecord = FALSE
    ),
    regexp = "n_batch.*positive|must be > 0"
  )
})
