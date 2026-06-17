local({
  required <- "4.6.0"
  current  <- paste(R.version$major, R.version$minor, sep = ".")
  if (current != required) {
    warning(
      "Este proyecto requiere R ", required,
      " (activo: R ", current, "). Cambia con: rig default ", required
    )
  }
})

source("renv/activate.R")
