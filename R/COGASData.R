#' @title
#' COGASData Class
#'
#' @description
#' Abstract base [R6::R6Class] for data from COGAS. This class acts as a
#' container for the process of reading, cleaning, and analyzing raw data from
#' GSMG measurements. There are general methods for cleaning and summarizing
#' depending on the instrument used and the compound of interest, propagating
#' measurement uncertainties, and visualization capabilities.
#'
#' @examples
COGASData <- R6::R6Class(
  classname = "COGASData",
  lock_class = TRUE,
  private = list(
    #' @field instrument (`character(1)`)
    #' Instrument used for the current analysis.
    instrument = NULL,

    #' @field compound (`character(1)`)
    #' Compound of interest for the current analysis.
    compound = NULL,

    #' @field method (`character(1)`)
    #' Method used to collect the data.
    method = NULL,

    #' @field control (`character(1)`)
    #' Name of the control sample used throughout the data file.
    control = NULL,

    #' @field zero_gas (`character(1)`)
    #' Name of the zero gas sample used throughout the data file.
    zero_gas = NULL,

    #' @field cf (`numeric(1)`)
    #' Correction factor to apply to sample responses.
    cf = NULL,

    #' @field ucf (`numeric(1)`)
    #' Uncertainty associated with the correction factor.
    ucf = NULL,

    #' @field response_mc (`matrix()`)
    #' Monte Carlo samples for responses.
    response_mc = NULL,

    #' @field control_mc (`matrix()`)
    #' Monte Carlo samples for interpolated control responses.
    control_mc = NULL,

    #' @field ratio_mc (`matrix()`)
    #' Monte Carlo samples for ratios.
    ratio_mc = NULL,

    #' @field raw_data (`data.frame()`)
    #' Raw data from the current analysis.
    raw_data = NULL,

    #' @field cleaned_data (`data.frame()`)
    #' Cleaned data from the current analysis.
    cleaned_data = NULL,

    #' @field summ_stats (`data.frame()`)
    #' Summary statistics associated with $cleaned_data.
    summ_stats = NULL,

    #' @field break_set_plot ([ggplot2::ggplot])
    #' Plot of the normalized responses and ratios for each break set.
    break_set_plot = NULL,

    #' @field replicate_plot ([ggplot2::ggplot])
    #' Plot of replicate ratios with uncertainty bars for each sample.
    replicate_plot = NULL,

    #' @description
    #' Make Column Names
    #'
    #' @details
    #' This function generates the column names for the data being read in
    #' based on the value of `private$instrument` and `private$method`.
    #'
    #' @return `character()`
    make_col_names = function() {
      default_names <- c("date_time", "port", "sample_id")
      data_names <- if (private$method == TIME_SERIES) {
        "response_ts"
      } else {
        c("response", "uresponse")
      }

      if (private$instrument == GAS_CHROM) {
        data_names <- c("file_name", "rt", data_names, "height")
      }

      c(default_names, data_names)
    },

    #' @description
    #' Read Data from CSV
    #'
    #' @details
    #' This function reads in the raw data from a CSV file and adds the
    #' following metadata variables:
    #'
    #' * `sample_set`: represents the individual sample measurements
    #' * `control_set`: represents the different measurements of the control
    #' * `break_set`: represents the separate sub-analyses performed
    #' * `day_set`: represents the distinct analytical days of an analysis
    #'
    #' If analyzing data from multiple analyses, it is expected that the data
    #' from each is contained in the same file, separated by a header row. This
    #' is how different break sets are determined.
    #'
    #' @param filepath (`character(1)`)
    #'   Path to the data file.
    #' @param ... Unused arguments.
    #'
    #' @return `data.frame()`
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

    #' @description
    #' Summarize Time Series Data
    #'
    #' @details
    #' This function takes the raw data from the current analysis and produces
    #' summary statistics for the associated time series data. Its main purpose
    #' is for GC data and the cases where data is collected using older COGAS
    #' versions, which don't automatically produce summaries of the data.
    #' Implementation of a robust time series framework will be provided in
    #' the future for modelling instrument response dynamics.
    #'
    #' @param raw (`data.frame()`)
    #'   Raw data for the current analysis.
    #' @param n (`numeric(1)`)
    #'   Number of data points to use from the end of each `sample_set`.
    #'   Defaults to 60.
    #' @param ... Unused arguments.
    #'
    #' @return `data.frame()`
    summarize_if_ts = function(raw, n = N_DEFAULT, ...) {
      if (private$method == TIME_SERIES) {
        raw %>%
          dplyr::slice_tail(n = n, by = sample_set) %>%
          dplyr::summarise(
            dplyr::across(date_time:day_set, dplyr::last),
            response = mean(response_ts),
            uresponse = sd(response_ts),
            .by = sample_set
          ) %>%
          dplyr::relocate(sample_set, .after = sample_id)
      } else {
        raw
      }
    },

    #' @description
    #' Perform Instrument-Specific Calculations
    #'
    #' @details
    #' This function takes the raw data from the current analysis and performs
    #' operations depending on the selected instrument.
    #'
    #' @param raw (`data.frame()`)
    #'   Raw data for the current analysis.
    #'
    #' @return `data.frame()`
    instrument_calc = function(raw) {
      if (private$instrument == HORIBA_NDIR) {
        # Adding quantization error to uncertainties before propagating.
        # These values are based on digital resolution of instrument
        # and assuming a uniform distribution for the uncertainty about any
        # particular response.
        resolution <- if (private$compound == CARBON_DIOXIDE) 0.01 else 0.1
        raw %>% dplyr::mutate(uresponse = sqrt(uresponse^2 + resolution^2 / 12))
      } else {
        raw
      }
    },

    #' @description
    #' Calculate Response Ratios
    #'
    #' @details
    #' This function takes the raw data from the current analysis and performs
    #' the calculation of response ratios and propagation of their respective
    #' uncertainties for the sampling sequence to implement a basic drift
    #' correction procedure. This procedure assumes the control sample acts as
    #' the conditional mean for the entire sampling sequence. Response ratios
    #' are calculated by dividing sample response measurements by the
    #' linearly-interpolated response measurements of the control. Uncertainty
    #' propagation is performed via a simple Monte Carlo simulation in which
    #' the responses for the analyte of interest are assumed to be normally
    #' distributed.
    #'
    #' @param raw (`data.frame()`)
    #'   Raw data for the current analysis.
    #' @param K (`numeric(1)`)
    #'   Number of Monte Carlo samples to draw. Defaults to 10,000.
    #' @param cf_samples (`character()`)
    #'   Samples to apply the correction factor to. Defaults to "none". Must be
    #'   one of the following:
    #'
    #'   * "none" to skip
    #'   * "all" for all samples
    #'   * "control" for just the control
    #'   * "not_control" for all samples but the control
    #'   * vector of specific sample ids to apply to only those samples
    #'
    #' @param cf_operation (`character(1)`)
    #'   Operation to use when applying the correction factor. Defaults to "*".
    #' @param ... Unused arguments.
    #'
    #' @return `data.frame()`
    calc_ratios = function(
      raw,
      K = K_DEFAULT,
      cf_samples = CF_SAMPLES_DEFAULT,
      cf_operation = CF_OPERATION_DEFAULT,
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
        matrix(ncol = K)

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
  ),
  public = list(
    #' @description
    #' Create a New [COGASData] Object
    #'
    #' @details
    #' This object is typically constructed via its derived classes [CalData]
    #' and [SampData]. Automatically calls `clean_data`.
    #'
    #' @param params (`list()`)
    #'   Parameters pertaining to the current analysis. Must contain the
    #'   following at minimum:
    #'
    #'   * instrument (`character(1)`): instrument used
    #'   * compound (`character(1)`): compound of interest
    #'   * control (`character(1)`): name of the control sample
    #'   * filepath (`character(1)`): passed to `private$read_data`
    #'
    #'   Additional parameters for other private attributes or arguments to
    #'   helper functions can also be given and will be used in the proper
    #'   location. Any additional parameters left unspecified will be given
    #'   default values.
    #'
    #' @return [COGASData]
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
        METHOD_DEFAULT
      }
      private$zero_gas <- if ("zero_gas" %in% names(params)) {
        stringr::str_to_upper(params$zero_gas)
      } else {
        ZERO_GAS_DEFAULT
      }
      private$cf <- if ("cf" %in% names(params)) {
        params$cf
      } else {
        CF_DEFAULT
      }
      private$ucf <- if ("ucf" %in% names(params)) {
        params$ucf
      } else {
        UCF_DEFAULT
      }

      rlang::inject(self$clean_data(
        !!!params[!names(params) %in% names(private)]
      ))
    },

    #' @description
    #' Print [COGASData] Object
    #'
    #' @return [COGASData]
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
        "    raw data: tibble [",
        paste0(dim(private$raw_data), collapse = " x "),
        "]\n",
        "    cleaned data: tibble [",
        paste0(dim(private$cleaned_data), collapse = " x "),
        "]\n",
        "    summary stats: tibble [",
        paste0(dim(private$summ_stats), collapse = " x "),
        "]\n",
        "    break set plot: ",
        ifelse(is.null(private$break_set_plot), "NULL\n", "ggplot\n"),
        "    replicates plot: ",
        ifelse(is.null(private$replicate_plot), "NULL\n", "ggplot\n"),
        sep = ""
      )

      invisible(self)
    },

    #' @description
    #' Clean Data
    #'
    #' @details
    #' This function reads in and cleans the data from the current analysis.
    #' Currently only has support for reading in data from a CSV file containing
    #' all of the raw data scraped from a COGAS Excel workbook. Based on the
    #' instrument used, compound of interest, and data collection method,
    #' different data operations may be performed.
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
    #' @return [COGASData]
    clean_data = function(...) {
      private$raw_data <- private$read_data(...) %>%
        private$summarize_if_ts(...) %>%
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
            c(private$control, private$zero_gas, DO_NOT_INCLUDE, PRECONDITION),
          control_set == 0
        )) %>%
        dplyr::mutate(dplyr::across(
          tidyselect::where(is.factor),
          forcats::fct_drop
        )) %>%
        dplyr::relocate(tidyselect::ends_with("ratio"), .after = day_set) %>%
        dplyr::arrange(sample_id, sample_set)

      invisible(self)
    },

    #' @description
    #' Summarize Data
    #'
    #' @details
    #' This function takes the cleaned data from the current analysis and
    #' computes summary statistics according to the supplied grouping
    #' variables.
    #'
    #' @param summ_cols (`character()`)
    #'   Variables to group by when performing calculation. Defaults to
    #'   "sample_id".
    #'
    #' @return [COGASData]
    summarize_data = function(summ_cols = SUMM_COLS_DEFAULT) {
      private$summ_stats <- private$cleaned_data %>%
        dplyr::summarise(
          avg_response = mean(response),
          sd_w = sqrt(mean(uresponse^2)),
          sd_b = stats::sd(response),
          sd_total = sqrt(sd_w^2 + sd_b^2),
          rsd_total = 100 * sd_total / avg_response,
          avg_ratio = mean(ratio),
          u_w = sqrt(mean(uratio^2)),
          u_b = stats::sd(ratio),
          u_total = sqrt(u_w^2 + u_b^2),
          ru_total = 100 * u_total / avg_ratio,
          n_repl = dplyr::n(),
          .by = tidyselect::all_of(summ_cols)
        ) %>%
        dplyr::relocate(tidyselect::all_of(summ_cols)) %>%
        dplyr::arrange(dplyr::pick(tidyselect::all_of(summ_cols)))

      invisible(self)
    },

    #' @description
    #' Plot Data
    #'
    #' @details
    #' This function generates useful diagnostic plots for the data from the
    #' current analysis.
    #'
    #' @return [COGASData]
    plot_data = function() {
      # Plotting ratios and responses, normalized to their respective averages,
      # for each `break_set` to understand how well the drift correction
      # procedure is working.
      private$break_set_plot <- private$raw_data %>%
        dplyr::filter_out(dplyr::when_any(
          sample_id %in% c(private$zero_gas, PRECONDITION),
          control_set == 0
        )) %>%
        dplyr::summarise(
          avg_response = mean(response),
          avg_ratio = mean(ratio),
          .by = c(break_set, sample_id)
        ) %>%
        dplyr::inner_join(
          private$raw_data,
          .,
          by = c("break_set", "sample_id")
        ) %>%
        dplyr::mutate(
          norm_response = response / avg_response,
          norm_ratio = ratio / avg_ratio,
          dplyr::across(tidyselect::where(is.factor), forcats::fct_drop)
        ) %>%
        ggplot2::ggplot(ggplot2::aes(
          x = as.numeric(sample_set),
          color = sample_id
        )) +
        ggplot2::geom_line(ggplot2::aes(
          y = norm_response,
          linetype = "Response"
        )) +
        ggplot2::geom_point(
          ggplot2::aes(y = norm_response, shape = "Response"),
          size = 2
        ) +
        ggplot2::geom_line(ggplot2::aes(y = norm_ratio, linetype = "Ratio")) +
        ggplot2::geom_point(
          ggplot2::aes(y = norm_ratio, shape = "Ratio"),
          size = 2
        ) +
        ggplot2::scale_color_viridis_d(name = "Sample ID", option = "turbo") +
        ggplot2::scale_linetype_manual(
          name = NULL,
          values = c(1, 3),
          breaks = c("Ratio", "Response")
        ) +
        ggplot2::scale_shape_manual(
          name = NULL,
          values = c(19, 1),
          breaks = c("Ratio", "Response")
        ) +
        ggplot2::guides(
          shape = ggplot2::guide_legend(order = 1),
          linetype = ggplot2::guide_legend(order = 1),
          color = ggplot2::guide_legend(override.aes = list(linetype = 0))
        ) +
        ggplot2::facet_wrap(
          ggplot2::vars(break_set),
          ncol = 2,
          scales = "free",
          dir = "v",
          labeller = "label_both"
        ) +
        ggplot2::labs(x = NULL, y = "Normalized Value") +
        ggplot2::theme_bw()

      # Getting unnested ratios for violin plots.
      unnested <- private$ratio_mc %>%
        tibble::as_tibble(
          .name_repair = "unique_quiet",
          rownames = "sample_set"
        ) %>%
        dplyr::left_join(private$cleaned_data, ., by = "sample_set") %>%
        tidyr::pivot_longer(
          tidyselect::starts_with("..."),
          names_to = "mc_repl",
          values_to = "ratio_mc"
        ) %>%
        dplyr::select(sample_id, sample_set, break_set, ratio_mc)

      # Plotting distribution of each replicate ratio by sample id, overlaid
      # with their location and uncertainty estimates.
      private$replicate_plot <- ggplot2::ggplot(
        mapping = ggplot2::aes(x = sample_set)
      ) +
        ggplot2::geom_violin(
          ggplot2::aes(y = ratio_mc, fill = break_set),
          data = unnested,
          quantile.linetype = 3,
          quantile.color = "black",
          scale = "width",
          alpha = 0.6
        ) +
        ggplot2::geom_pointrange(
          ggplot2::aes(y = ratio, ymin = ratio - uratio, ymax = ratio + uratio),
          data = private$cleaned_data,
          size = 0.25
        ) +
        ggplot2::scale_fill_viridis_d(name = "Break Set") +
        ggplot2::facet_wrap(
          ggplot2::vars(sample_id),
          ncol = 2,
          scales = "free"
        ) +
        ggplot2::labs(x = "Sample Set", y = "Ratio") +
        ggplot2::theme_bw()

      invisible(self)
    },

    #' @description
    #' Get Analysis Info
    #'
    #' @details
    #' This functions returns a list of all private attributes that contain
    #' information related to the analysis. If all attributes are not needed at
    #' the same time, it can be convenient to select a specific one from the
    #' function call, i.e. `get_info()$instrument`.
    #'
    #' @return `list()`
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
    #' Get MC Samples
    #'
    #' @details
    #' This function returns a list of all private attributes containing MC
    #' samples. If all attributes are not needed at the same time, it can be
    #' convenient to select a specific one from the function call, i.e.
    #' `get_mc()$response_mc`.
    #'
    #' @return `list()`
    get_mc = function() {
      list(
        response_mc = private$response_mc,
        control_mc = private$control_mc,
        ratio_mc = private$ratio_mc
      )
    },

    #' @description
    #' Get Data Frames
    #'
    #' @details
    #' This function returns a list of all private attributes that are data
    #' frames. If all attributes are not needed at the same time, it can be
    #' convenient to select a specific one from the function call, i.e.
    #' `get_data()$cleaned_data`.
    #'
    #' @return `list()`
    get_data = function() {
      list(
        raw_data = private$raw_data,
        cleaned_data = private$cleaned_data,
        summ_stats = private$summ_stats
      )
    },

    #' @description
    #' Get Plots
    #'
    #' @details
    #' This function returns a list of all private attributes that are plots.
    #' If all attributes are not needed at the same time, it can be convenient
    #' to select a specific one from the function call, i.e.
    #' `get_plots()$break_set_plot`.
    #'
    #' @return `list()`
    get_plots = function() {
      list(
        break_set_plot = private$break_set_plot,
        replicate_plot = private$replicate_plot
      )
    },

    #' @description
    #' Save Object
    #'
    #' @details
    #' This function saves the object as an RDS file.
    #'
    #' @param filepath (`character(1)`)
    #'   Path to save the object to.
    #'
    #' @return [COGASData]
    save_object = function(filepath) {
      saveRDS(self, file = filepath)

      invisible(self)
    }
  )
)
