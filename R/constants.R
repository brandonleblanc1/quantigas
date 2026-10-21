#' Miscellaneous Constants
#'
#' A locked environment holding fixed constants for use across the package.
#'
#' @export
misc <- new.env(parent = emptyenv())

misc$SUMMARIZED <- "SUMMARIZED"
misc$TIME_SERIES <- "TIME_SERIES"
misc$DO_NOT_INCLUDE <- "DNI"
misc$PRECONDITION <- "PRECONDITION"
misc$NU_MIN <- 3
misc$NU_MAX <- 100

lockEnvironment(misc, bindings = TRUE)

#' Default Values
#'
#' A locked environment holding fixed constants for use across the package.
#'
#' @export
defaults <- new.env(parent = emptyenv())

defaults$N <- 60
defaults$K <- 10000
defaults$CF_SAMPLES <- "none"
defaults$CF_OPERATION <- "*"
defaults$METHOD <- misc$SUMMARIZED
defaults$ZERO_GAS <- "N2"
defaults$CF <- 1
defaults$UCF <- 0
defaults$SUMM_COLS <- "sample_id"
defaults$MODEL_NUM <- 1
defaults$UNITS <- "mu*mol ~~ mol^-1"
defaults$FONTFAMILY <- "sans"
defaults$FONTSIZE <- 11

lockEnvironment(defaults, bindings = TRUE)

#' Compound Options
#'
#' A locked environment holding fixed constants for use across the package.
#'
#' @export
compounds <- new.env(parent = emptyenv())

compounds$ARGON <- "7440-37-1"
compounds$CARBON_DIOXIDE <- "124-38-9"
compounds$CARBON_MONOXIDE <- "630-08-0"
compounds$HELIUM <- "7440-59-7"
compounds$HYDROGEN <- "1333-74-0"
compounds$HYDROGEN_SULFIDE <- "7783-06-4"
compounds$METHANE <- "74-82-8"
compounds$NITROGEN_MONOXIDE <- "10102-43-9"
compounds$NITROGEN_DIOXIDE <- "10102-44-0"
compounds$NITROGEN_OXIDES <- "nox"
compounds$OXYGEN <- "7782-44-7"
compounds$PROPANE <- "74-98-6"
compounds$SULFUR_DIOXIDE <- "7446-09-5"

lockEnvironment(compounds, bindings = TRUE)

#' Instrument Options
#'
#' A locked environment holding fixed constants for use across the package.
#'
#' @export
instruments <- new.env(parent = emptyenv())

instruments$HORIBA_NDIR <- "HORIBA_NDIR"
instruments$XSTREAM_NDIR <- "XSTREAM_NDIR"
instruments$AMETEK_NDUV <- "AMETEK_NDUV"
instruments$LOW_CHEMI <- "LOW_CHEMI"
instruments$HIGH_CHEMI <- "HIGH_CHEMI"
instruments$GAS_CHROM <- "GAS_CHROM"

lockEnvironment(instruments, bindings = TRUE)
