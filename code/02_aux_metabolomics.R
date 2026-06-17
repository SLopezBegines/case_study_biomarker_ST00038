# =============================================================================
# I/O helpers for Metabolomics Workbench REST "MSdata" tables.
# Format: row 1 = header (Metabolite_name, RefMet_name, <sample IDs...>);
#         one row labelled "Factors" holding "Key:Value | Key:Value | ..."
#         per sample; remaining rows = one identified metabolite each.
# All comments in English.
# =============================================================================

# Parse a single Factors cell like
#   "Organ:Plasma | Cancer status:Adenocarcinoma | Smoker:Current | Gender:F"
# into a named character vector with canonical column names.
.parse_factor_cell <- function(s) {
  out <- c(Organ = NA_character_, Cancer_status = NA_character_,
           Smoker = NA_character_, Gender = NA_character_)
  if (is.na(s)) return(out)
  for (part in strsplit(s, "\\|")[[1]]) {
    part <- trimws(part)
    idx <- regexpr(":", part, fixed = TRUE)
    if (idx > 0) {
      k <- trimws(substr(part, 1, idx - 1))
      v <- trimws(substr(part, idx + 1, nchar(part)))
      key <- switch(k,
                    "Organ"         = "Organ",
                    "Cancer status" = "Cancer_status",
                    "Smoker"        = "Smoker",
                    "Gender"        = "Gender",
                    NA_character_)
      if (!is.na(key)) out[key] <- v
    }
  }
  out
}

# Read one MSdata file. Returns a list with:
#   $mat     : numeric matrix [metabolites x samples], rownames = Metabolite_name
#   $refmet  : data.frame mapping Metabolite_name -> RefMet_name
#   $factors : data.frame [samples x {Organ, Cancer_status, Smoker, Gender}]
read_msdata <- function(path) {
  df <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "")
  names(df) <- trimws(names(df))
  metcol <- names(df)[1]
  refcol <- names(df)[2]
  df[[metcol]] <- trimws(df[[metcol]])
  sample_ids <- names(df)[-(1:2)]

  fac_row   <- df[df[[metcol]] == "Factors", , drop = FALSE]
  data_rows <- df[df[[metcol]] != "Factors", , drop = FALSE]
  # Drop empty/blank spacer rows (no metabolite name -> cannot be keyed)
  data_rows <- data_rows[!is.na(data_rows[[metcol]]) &
                           nzchar(data_rows[[metcol]]), , drop = FALSE]

  # Factors -> data.frame
  fac_strings <- as.character(unlist(fac_row[1, sample_ids]))
  fac <- t(vapply(fac_strings, .parse_factor_cell,
                  FUN.VALUE = c(Organ = "", Cancer_status = "",
                                Smoker = "", Gender = "")))
  rownames(fac) <- sample_ids
  factors_df <- as.data.frame(fac, stringsAsFactors = FALSE)

  # Intensity matrix (metabolites x samples), coerced to numeric
  mat <- as.matrix(data_rows[, sample_ids, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- data_rows[[metcol]]

  refmet <- data.frame(Metabolite_name = data_rows[[metcol]],
                       RefMet_name     = data_rows[[refcol]],
                       stringsAsFactors = FALSE)

  list(mat = mat, refmet = refmet, factors = factors_df)
}
