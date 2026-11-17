#' [R6::R6Class] for Processing COGAS Data
#'
#' Base class for processing data from COGAS. Essentially just a container for
#' the process of reading, cleaning, and analyzing raw data from GSMG
#' measurements. This object is typically constructed via its derived classes,
#' [CalData] and [SampData] (and variants).
#'
#' @export
COGASData <- R6::R6Class(
  classname = "COGASData",
  lock_class = TRUE,
  public = list(
    #' @description
    #' Create a new `COGASData` object. Automatically calls
    #' [`clean_data()`](#method-clean_data).
    #'
    #' @param params (`list`)
    #'   Parameters pertaining to the current analysis. Must contain the
    #'   following at minimum:
    #'
    #'   * instrument (`character(1)`): instrument used
    #'   * compound (`character(1)`): compound of interest
    #'   * control (`character(1)`): name of the control sample
    #'   * filepath (`character(1)`): passed to `private$read_data()`
    #'
    #'   Additional parameters for other private attributes or arguments to
    #'   helper functions can also be given and will be used in the proper
    #'   location. Any additional parameters left unspecified will be given
    #'   default values.
    #'
    #' @return A new `COGASData` object.
    initialize = function(params) {
      # Validating inputs.
      stopifnot(
        length(params) > 0,
        c("instrument", "compound", "control", "filepath") %in% names(params)
      )

      private$instrument <- params$instrument
      private$compound <- params$compound
      private$control <- stringr::str_to_upper(params$control)
      private$method <- if ("method" %in% names(params)) {
        params$method
      } else {
        defaults$METHOD
      }
      private$zero_gas <- if ("zero_gas" %in% names(params)) {
        stringr::str_to_upper(params$zero_gas)
      } else {
        defaults$ZERO_GAS
      }
      private$cf <- if ("cf" %in% names(params)) {
        params$cf
      } else {
        defaults$CF
      }
      private$ucf <- if ("ucf" %in% names(params)) {
        params$ucf
      } else {
        defaults$UCF
      }

      rlang::inject(self$clean_data(
        !!!params[!names(params) %in% names(private)]
      ))
    },

    #' @description
    #' Print useful information about the state of the object.
    #'
    #' @return The calling object invisibly.
    print = function() {
      cat(
        "<COGASData>\n",
        "  Private Attributes:\n",
        "    instrument: ",
        private$instrument,
        "\n",
        "    compound: ",
        private$compound,
        "\n",
        "    collection method: ",
        private$method,
        "\n",
        "    control: ",
        private$control,
        "\n",
        "    zero gas: ",
        private$zero_gas,
        "\n",
        "    correction factor: ",
        private$cf,
        "\n",
        "    correction factor u: ",
        private$ucf,
        "\n",
        "    response mc samples: matrix [",
        paste0(dim(private$response_mc), collapse = " x "),
        "]\n",
        "    control mc samples: matrix [",
        paste0(dim(private$control_mc), collapse = " x "),
        "]\n",
        "    ratio mc samples: matrix [",
        paste0(dim(private$ratio_mc), collapse = " x "),
        "]\n",
        "    response covariance: matrix [",
        paste0(dim(private$response_sigma), collapse = " x "),
        "]\n",
        "    ratio covariance: matrix [",
        paste0(dim(private$ratio_sigma), collapse = " x "),
        "]\n",
        "    raw data: tbl_df [",
        paste0(dim(private$raw_data), collapse = " x "),
        "]\n",
        "    cleaned data: tbl_df [",
        paste0(dim(private$cleaned_data), collapse = " x "),
        "]\n",
        "    summary stats: tbl_df [",
        paste0(dim(private$summ_stats), collapse = " x "),
        "]\n",
        "    drift plot: ",
        ifelse(is.null(private$drift_plot), "NULL\n", "ggplot\n"),
        "    replicates plot: ",
        ifelse(is.null(private$replicate_plot), "NULL\n", "ggplot\n"),
        sep = ""
      )

      invisible(self)
    },

    #' @description
    #' Read in and clean the data from the current analysis. Currently only has
    #' support for reading data from a CSV file containing all of the raw data
    #' scraped from a COGAS Excel workbook. Based on the instrument used,
    #' compound of interest, and data collection method, different data
    #' operations may be performed.
    #'
    #' One thing to note is that unless an entire break set went awry, any data
    #' being discarded from the analysis should be left intact with the rest of
    #' the data. This "bad" data should be discriminated against by renaming the
    #' samples to "DNI" (Do Not Include). This ensures that any metadata
    #' associated with the "bad" data is still accurately reflected in the
    #' "good" data, and these data points can be easily filtered out.
    #'
    #' @param ... Additional named arguments to pass to helper functions.
    #'
    #' @return The calling object invisibly.
    clean_data = function(...) {
      private$raw_data <- private$read_data(...) %>%
        private$method_summ(...) %>%
        private$instrument_calc() %>%
        private$calc_ratios(...) %>%
        dplyr::mutate(dplyr::across(port:day_set, forcats::as_factor))

      # Generating the cleaned data set. Any data labeled "DNI" gets filtered
      # out, but these entries are initially necessary to ensure metadata is
      # correct.
      private$cleaned_data <- private$raw_data %>%
        dplyr::select(!tidyselect::ends_with("control")) %>%
        dplyr::filter_out(dplyr::when_any(
          sample_id %in%
            c(
              private$control,
              private$zero_gas,
              misc$DO_NOT_INCLUDE,
              misc$PRECONDITION
            ),
          control_set == 0
        )) %>%
        dplyr::mutate(dplyr::across(
          tidyselect::where(is.factor),
          forcats::fct_drop
        )) %>%
        dplyr::relocate(tidyselect::ends_with("ratio"), .after = day_set) %>%
        dplyr::arrange(sample_id, sample_set)

      # Finding subset of raw data that corresponds to samples kept in the
      # cleaned dataset.
      cleaned_subset <- as.numeric(as.character(
        private$cleaned_data$sample_set
      ))

      private$response_sigma <- robustbase::covMcd(t(private$response_mc[
        cleaned_subset,
      ]))$cov
      private$ratio_sigma <- robustbase::covMcd(t(private$ratio_mc[
        cleaned_subset,
      ]))$cov

      invisible(self)
    },

    #' @description
    #' Compute summary statistics on the cleaned data according to the supplied
    #' grouping variables.
    #'
    #' @param summ_cols (`character`)
    #'   Variables to group by when performing calculation. Defaults to
    #'   "sample_id".
    #'
    #' @return The calling object invisibly.
    summarize_data = function(summ_cols = defaults$SUMM_COLS) {
      # Gathering covariance matrices and centering MC samples based on cleaned
      # data subset.
      sigma <- list(
        response = private$response_sigma,
        ratio = private$ratio_sigma
      )
      centered <- purrr::imap(
        list(response = private$response_mc, ratio = private$ratio_mc),
        \(mat, name) {
          sweep(
            mat[as.numeric(as.character(private$cleaned_data$sample_set)), ],
            1,
            private$cleaned_data[[name]]
          )
        }
      )

      # Finding subsets of cleaned data dictated by the summary columns.
      var_subsets <- split(
        seq_len(nrow(private$cleaned_data)),
        purrr::map(summ_cols, \(var) private$cleaned_data[[var]]),
        drop = TRUE,
        lex.order = TRUE
      )

      private$summ_stats <- purrr::map(var_subsets, \(sub) {
        # Calculating weighted average and the corresponding variance for each
        # subset.
        n <- length(sub)
        ones <- rep(1, n)
        summary <- purrr::imap(sigma, \(mat, name) {
          sigma_inv <- MASS::ginv(mat[sub, sub])
          denom <- as.numeric(t(ones) %*% sigma_inv %*% ones)
          w <- sigma_inv %*% ones / denom

          list(
            mu = as.numeric(t(w) %*% private$cleaned_data[[name]][sub]),
            var = 1 / denom
          )
        })

        tibble::tibble(
          avg_response = summary$response$mu,
          sd_w = robustbase::Qn(centered$response[sub, ]),
          sd_b = sqrt(sum(sigma$response[sub, sub])),
          sd_total = sqrt(summary$response$var),
          rsd_total = 100 * sd_total / avg_response,
          avg_ratio = summary$ratio$mu,
          u_w = robustbase::Qn(centered$ratio[sub, ]),
          u_b = sqrt(sum(sigma$ratio[sub, sub])),
          u_total = sqrt(summary$ratio$var),
          ru_total = 100 * u_total / avg_ratio,
          n_repl = n
        )
      }) %>%
        dplyr::bind_rows(.id = "subset") %>%
        tidyr::separate_wider_delim(subset, ".", names = summ_cols)

      invisible(self)
    },

    #' @description
    #' Generate useful diagnostic plots for the current analysis.
    #'
    #' @return The calling object invisibly.
    plot_data = function() {
      # Plotting ratios and responses, normalized to their respective averages,
      # for each `break_set` to understand how well the drift correction
      # procedure is working.
      private$drift_plot <- private$raw_data %>%
        dplyr::filter_out(dplyr::when_any(
          sample_id %in% c(private$zero_gas, misc$PRECONDITION),
          control_set == 0
        )) %>%
        dplyr::mutate(
          norm_response = response / mean(response),
          norm_ratio = ratio / mean(ratio),
          .by = c(break_set, sample_id)
        ) %>%
        dplyr::mutate(dplyr::across(
          tidyselect::where(is.factor),
          forcats::fct_drop
        )) %>%
        tidyr::pivot_longer_spec(tibble::tribble(
          ~.name          , ~.value , ~type      ,
          "norm_response" , "value" , "Response" ,
          "norm_ratio"    , "value" , "Ratio"
        )) %>%
        ggplot2::ggplot(ggplot2::aes(
          x = as.numeric(sample_set),
          y = value,
          color = sample_id
        )) +
        ggplot2::geom_line(ggplot2::aes(
          linetype = type,
          group = interaction(break_set, sample_id, type)
        )) +
        ggplot2::geom_point(ggplot2::aes(shape = type), size = 2) +
        ggplot2::scale_color_viridis_d(name = "Sample ID", option = "turbo") +
        ggplot2::scale_linetype_manual(name = NULL, values = c(1, 3)) +
        ggplot2::scale_shape_manual(name = NULL, values = c(19, 1)) +
        ggplot2::guides(
          shape = ggplot2::guide_legend(order = 1),
          linetype = ggplot2::guide_legend(order = 1),
          color = ggplot2::guide_legend(override.aes = list(linetype = 0))
        ) +
        ggplot2::labs(x = NULL, y = "Normalized Value") +
        ggplot2::theme_bw()

      # Getting unnested ratios for violin plots.
      unnested <- private$ratio_mc[
        as.numeric(as.character(private$cleaned_data$sample_set)),
      ] %>%
        tibble::as_tibble(.name_repair = "unique_quiet") %>%
        dplyr::bind_cols(private$cleaned_data) %>%
        tidyr::pivot_longer(tidyselect::starts_with("..."))

      # Plotting distribution of each replicate ratio by sample id, overlaid
      # with their location and uncertainty estimates.
      private$replicate_plot <- ggplot2::ggplot(
        mapping = ggplot2::aes(x = sample_id)
      ) +
        ggplot2::geom_violin(
          ggplot2::aes(y = value, fill = break_set, group = sample_set),
          data = unnested,
          position = "identity",
          quantile.linetype = 3,
          quantile.color = "black",
          alpha = 0.6
        ) +
        ggplot2::geom_pointrange(
          ggplot2::aes(y = ratio, ymin = ratio - uratio, ymax = ratio + uratio),
          data = private$cleaned_data,
          size = 0.25
        ) +
        ggplot2::scale_fill_viridis_d(name = "Break Set") +
        ggplot2::labs(x = NULL, y = "Ratio") +
        ggplot2::theme_bw()

      invisible(self)
    },

    #' @description
    #' Get all private attributes that contain general information on the
    #' analysis. If all attributes are not needed at the same time, it can be
    #' convenient to select a specific one from the function call, i.e.
    #' `get_info()$instrument`.
    #'
    #' @return A `list` of attributes.
    get_info = function() {
      list(
        instrument = private$instrument,
        compound = private$compound,
        method = private$method,
        control = private$control,
        zero_gas = private$zero_gas,
        cf = private$cf,
        ucf = private$ucf
      )
    },

    #' @description
    #' Get all private attributes containing MC samples. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_mc()$response_mc`.
    #'
    #' @return A `list` of attributes.
    get_mc = function() {
      list(
        response_mc = private$response_mc,
        control_mc = private$control_mc,
        ratio_mc = private$ratio_mc,
        response_sigma = private$response_sigma,
        ratio_sigma = private$ratio_sigma
      )
    },

    #' @description
    #' Get all private attributes that are data frames. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_data()$cleaned_data`.
    #'
    #' @return A `list` of attributes.
    get_data = function() {
      list(
        raw_data = private$raw_data,
        cleaned_data = private$cleaned_data,
        summ_stats = private$summ_stats
      )
    },

    #' @description
    #' Get all private attributes that are plot objects. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_plots()$drift_plot`.
    #'
    #' @return A `list` of attributes.
    get_plots = function() {
      list(
        drift_plot = private$drift_plot,
        replicate_plot = private$replicate_plot
      )
    },

    #' @description
    #' Save the object to an RDS file.
    #'
    #' @param filepath (`character(1)`)
    #'   Path to save the object to.
    #'
    #' @return The calling object invisibly.
    save_object = function(filepath) {
      saveRDS(self, file = filepath)

      invisible(self)
    }
  ),
  private = list(
    # `character(1)`: Instrument used for the current analysis.
    instrument = NULL,

    # `character(1)`: Compound of interest for the current analysis.
    compound = NULL,

    # `character(1)`: Method used to collect the data.
    method = NULL,

    # `character(1)`: Name of the control sample used throughout the data file.
    control = NULL,

    # `character(1)`: Name of the zero gas sample used throughout the data file.
    zero_gas = NULL,

    # `numeric(1)`: Correction factor to apply to sample responses.
    cf = NULL,

    # `numeric(1)`: Uncertainty associated with the correction factor.
    ucf = NULL,

    # `matrix`: Monte Carlo samples for responses.
    response_mc = NULL,

    # `matrix`: Monte Carlo samples for interpolated control responses.
    control_mc = NULL,

    # `matrix`: Monte Carlo samples for ratios.
    ratio_mc = NULL,

    # `matrix`: Covariance matrix for responses.
    response_sigma = NULL,

    # `matrix`: Covariance matrix for ratios
    ratio_sigma = NULL,

    # [tibble::tibble]: Raw data from the current analysis.
    raw_data = NULL,

    # [tibble::tibble]: Cleaned data from the current analysis.
    cleaned_data = NULL,

    # [tibble::tibble]: Summary statistics on the cleaned data.
    summ_stats = NULL,

    # [ggplot2::ggplot]: Plot of the normalized responses and ratios for each
    # break set.
    drift_plot = NULL,

    # [ggplot2::ggplot]: Plot of replicate ratios with uncertainty bars for
    # each sample.
    replicate_plot = NULL,

    make_col_names = function() {
      default_names <- c("date_time", "port", "sample_id")
      data_names <- if (private$method == misc$TIME_SERIES) {
        "response_ts"
      } else {
        c("response", "uresponse")
      }

      if (private$instrument == instruments$GAS_CHROM) {
        data_names <- c("file_name", "rt", data_names, "height")
      }

      c(default_names, data_names)
    },

    # This function reads in the raw data from a CSV file and adds the
    # following metadata variables:
    #
    # * `sample_set`: represents the individual sample measurements
    # * `control_set`: represents the different measurements of the control
    # * `break_set`: represents the separate sub-analyses performed
    # * `day_set`: represents the distinct analytical days of an analysis
    #
    # If analyzing data from multiple analyses, it is expected that the data
    # from each is contained in the same file, separated by a header row. This
    # is how different break sets are determined.
    read_data = function(filepath, ...) {
      readr::read_file(filepath) %>%
        stringr::str_to_upper() %>%
        # Searching for lines that contain the header row, which should contain
        # "DATE".
        stringr::str_split_1("DATE.+\\r?\\n{1}") %>%
        purrr::discard(\(txt) txt == "") %>%
        # Reading in each chunk of data found.
        purrr::map(\(txt) {
          readr::read_csv(
            I(txt),
            col_names = private$make_col_names(),
            show_col_types = FALSE
          )
        }) %>%
        dplyr::bind_rows(.id = "break_set") %>%
        dplyr::mutate(
          # Removing any newline characters still present from some of the GC
          # processing software.
          dplyr::across(tidyselect::where(is.character), \(txt) {
            stringr::str_remove_all(txt, "\\n")
          }),
          break_set = as.numeric(break_set),
          date = as.Date(date_time),
          is_diff_break = break_set != dplyr::lag(break_set),
          # `sample_set` is determined from differences between consecutive
          # `sample_id` values.
          is_diff_samp = sample_id != dplyr::lag(sample_id),
          # `control_set` is determined from measurements of the control sample.
          is_control = sample_id == private$control,
          # Replacing resulting `NA` values from `lag`.
          dplyr::across(tidyselect::where(is.logical), \(col) {
            dplyr::replace_values(col, NA ~ TRUE)
          }),
          sample_set = cumsum(is_diff_break | is_diff_samp),
          control_set = cumsum(is_control & (is_diff_break | is_diff_samp))
        ) %>%
        # `day_set` is just the different `date` groupings.
        dplyr::mutate(day_set = dplyr::cur_group_id(), .by = date) %>%
        dplyr::relocate(
          c(sample_set, control_set, break_set, day_set),
          .after = sample_id
        ) %>%
        dplyr::select(!date:is_control)
    },

    # This function takes the raw data from the current analysis and produces
    # summary statistics for the associated time series data. Implementation of
    # a robust time series framework will be provided in the future for
    # modelling instrument response dynamics.
    method_summ = function(raw, n = defaults$N, ...) {
      if (private$method == misc$TIME_SERIES) {
        raw %>%
          dplyr::slice_tail(n = n, by = sample_set) %>%
          dplyr::summarise(
            dplyr::across(date_time:day_set, dplyr::last),
            response = mean(response_ts),
            uresponse = sd(response_ts) / sqrt(n),
            .by = sample_set
          ) %>%
          dplyr::relocate(sample_set, .after = sample_id)
      } else {
        raw %>% dplyr::mutate(uresponse = uresponse / sqrt(n))
      }
    },

    # This function takes the raw data from the current analysis and performs
    # operations depending on the selected instrument.
    instrument_calc = function(raw) {
      if (private$instrument == instruments$HORIBA_NDIR) {
        # Adding quantization error to uncertainties before propagating.
        # These values are based on digital resolution of instrument
        # and assuming a uniform distribution for the uncertainty about any
        # particular response.
        resolution <- if (private$compound == compounds$CARBON_DIOXIDE) {
          0.01
        } else {
          0.1
        }
        raw %>% dplyr::mutate(uresponse = sqrt(uresponse^2 + resolution^2 / 12))
      } else if (
        private$instrument %in% c(instruments$LOW_CHEMI, instruments$HIGH_CHEMI)
      ) {
        raw %>% dplyr::mutate(uresponse = sqrt(uresponse^2 + (1e-5)^2 / 12))
      } else {
        raw
      }
    },

    # This function takes the raw data from the current analysis and performs
    # the calculation of response ratios and propagation of their respective
    # uncertainties for the sampling sequence to implement a basic drift
    # correction procedure. This procedure assumes the control sample acts as
    # the conditional mean for the entire sampling sequence. Response ratios
    # are calculated by dividing sample response measurements by the
    # linearly-interpolated response measurements of the control. Uncertainty
    # propagation is performed via a simple Monte Carlo simulation in which
    # the responses for the analyte of interest are assumed to be normally
    # distributed.
    calc_ratios = function(
      raw,
      K = defaults$K,
      cf_samples = defaults$CF_SAMPLES,
      cf_operation = defaults$CF_OPERATION,
      ...
    ) {
      # Getting necessary indices for interpolation and indices to apply
      # correction factor to, if specified.
      x_out <- seq_len(nrow(raw))
      x_control <- which(raw$sample_id == private$control)
      x_cf <- switch(
        cf_samples,
        none = NULL,
        all = x_out,
        control = x_control,
        not_control = which(raw$sample_id != private$control),
        which(raw$sample_id %in% cf_samples)
      )

      # Monte Carlo simulation for uncertainty propagation; assuming responses
      # and correction factors come from a normal distribution.
      private$response_mc <- stats::rnorm(
        nrow(raw) * K,
        mean = raw$response,
        sd = raw$uresponse
      ) %>%
        matrix(ncol = K, dimnames = list(raw$sample_id))

      if (!is.null(x_cf)) {
        cf_mc <- stats::rnorm(K, mean = private$cf, sd = private$ucf)

        # Applying correction factor to specified MC samples.
        private$response_mc[x_cf, ] <- sweep(
          private$response_mc[x_cf, ],
          2,
          cf_mc,
          FUN = cf_operation
        )

        # Applying correction factor to specified samples in original data.
        raw$response[x_cf] <- eval(call(
          cf_operation,
          raw$response[x_cf],
          private$cf
        ))
      }

      # Performing interpolation on each column of `response_mc`. Each column
      # corresponds to a MC sample of the entire analytical sequence so that the
      # rows match up with the raw data.
      private$control_mc <- apply(private$response_mc, 2, \(col) {
        stats::approx(
          x = x_control,
          y = col[x_control],
          xout = x_out,
          na.rm = FALSE
        )$y
      })

      # Calculating response ratios.
      private$ratio_mc <- private$response_mc / private$control_mc

      # Interpolating control vals and calculating ratios in original data,
      # as well as their respective uncertainty estimates from the simulation.
      raw %>%
        dplyr::mutate(
          control = stats::approx(
            x = x_control,
            y = response[x_control],
            xout = x_out,
            na.rm = FALSE
          )$y,
          ucontrol = apply(private$control_mc, 1, robustbase::Qn),
          ratio = response / control,
          uratio = apply(private$ratio_mc, 1, robustbase::Qn)
        )
    }
  )
)
