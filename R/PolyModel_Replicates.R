#' [R6::R6Class] for Model Building
#'
#' Derived class for modelling GSMG measurement data. Specifically deals with
#' building an error-in-variables model with all replicates data. It inherits
#' from the [PolyModel] class and thus has access to its public and private
#' attributes and methods.
#'
#' @export
PolyModel_Replicates <- R6::R6Class(
  classname = "PolyModel_Replicates",
  lock_class = TRUE,
  inherit = PolyModel,
  public = list(
    #' @description
    #' Create a new `PolyModel_Replicates` object.
    #'
    #' @param cal_obj ([CalData])
    #'   Calibration data object.
    #' @param filt_samples (`character`)
    #'   Sample(s) of interest to exclude from the model. Defaults to an empty
    #'   string.
    #'
    #' @return A new `PolyModel_Replicates` object.
    initialize = function(cal_obj, filt_samples = "") {
      stopifnot(inherits(cal_obj, "CalData"))

      cal_data <- cal_obj$get_data()
      cal_mc <- cal_obj$get_mc()

      # Combining certification data with cleaned data.
      private$cal_data <- cal_data$cleaned_data %>%
        dplyr::right_join(cal_data$cert_data, by = "sample_id") %>%
        dplyr::filter_out(sample_id %in% filt_samples) %>%
        dplyr::select(
          sample_id,
          sample_set,
          r = ratio,
          ur = uratio,
          c = c_cert,
          uc = uc_cert
        ) %>%
        dplyr::arrange(c, sample_set) %>%
        dplyr::mutate(
          dplyr::across(tidyselect::where(is.factor), forcats::fct_drop),
          sample_id = forcats::fct_reorder(sample_id, c),
          sample_set = forcats::fct_inorder(sample_set)
        )

      # Getting initial values for rho and sigma, as well as lower and upper
      # bounds for rho.
      private$init_vals <- private$cal_data %>%
        dplyr::summarise(
          rho0 = mean(r),
          sigma0 = stats::sd(r),
          rho_lb = rho0 - 6 * sigma0,
          rho_ub = rho0 + 6 * sigma0,
          dplyr::across(c(c, uc), dplyr::last),
          .by = sample_id
        )

      private$counts <- private$cal_data %>% dplyr::count(sample_id)
      private$n_std <- nrow(private$counts)
      private$n_repl <- sum(private$counts$n)

      # Subsetting ratio MC samples for relevant data.
      private$r_mc <- cal_mc$ratio_mc[
        as.numeric(levels(private$cal_data$sample_set)),
      ] %>%
        `rownames<-`(private$cal_data$sample_id)
    },

    #' @description
    #' Print useful information about the state of the object.
    #'
    #' @return The calling object invisibly.
    print = function() {
      super$print()

      cat(
        "<PolyModel_Replicates>\n",
        "  Private Attributes:\n",
        "    number standards: ",
        private$n_std,
        "\n",
        "    number replicates: ",
        private$n_repl,
        "\n",
        "    selected model: ",
        private$selected_model,
        "\n",
        "    dark uncertainty mc samples: numeric [",
        length(private$lambda_mc),
        "]\n",
        "    counts data: tbl_df [",
        paste0(dim(private$counts), collapse = " x "),
        "]\n",
        sep = ""
      )

      invisible(self)
    },

    #' @description
    #' Build candidate polynomial models and run the optimization algorithm for
    #' each to assess their fit of the data. It also generates residual plots
    #' for each candidate model. The best model has the smallest signed value
    #' of AIC, BIC, or both; however, the model selected for use is usually the
    #' one in which the values of AIC and BIC show their first substantial drop
    #' as model complexity increases. These statistical criteria should be
    #' taken into account in conjunction with the residual plots, which should
    #' show no discernable patterns.
    #'
    #' @return The calling object invisibly.
    model_selection = function() {
      # Setting up for parallel processing.
      future::plan(future::multisession, workers = future::availableCores() - 2)

      private$models <- private$build_models() %>%
        dplyr::mutate(
          # Running the optimization algorithm for each candidate model.
          opt = furrr::future_map2(forms, funcs, \(form, func) {
            private$run_opt(form, func)
          })
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          # Extracting parameters and calculating performance metrics.
          rho = list(opt$solution[1:private$n_std]),
          sigma = list(exp(opt$solution[
            (private$n_std + 1):(2 * private$n_std)
          ])),
          nu = opt$solution[(2 * private$n_std + 1)],
          tau = exp(opt$solution[(2 * private$n_std + 2)]),
          beta = list(opt$solution[-(1:(2 * private$n_std + 2))]),
          xi = list(vapply(
            rho,
            \(r) eval(str2expression(funcs))(r, beta),
            numeric(1)
          )),
          r_resid = list(private$cal_data$r - rep(rho, private$counts$n)),
          r_zscore = list(
            r_resid /
              sqrt(
                private$cal_data$ur^2 +
                  rep(sigma^2 * nu / (nu - 2), private$counts$n)
              )
          ),
          c_resid = list(private$init_vals$c - xi),
          c_zscore = list(c_resid / sqrt(private$init_vals$uc^2 + tau^2)),
          n_params = length(opt$solution),
          objective = opt$objective,
          aic = 2 * (n_params + opt$objective),
          bic = log(private$n_repl) * n_params + 2 * opt$objective,
          rmse = sqrt(mean(c_resid^2)),
          rse = sqrt(sum(c_resid^2) / (private$n_repl - n_params - 1)),
          r_sq = 1 -
            sum(c_resid^2) /
              sum((private$init_vals$c - mean(private$init_vals$c))^2),
          adj_r_sq = 1 -
            (1 - r_sq) * (private$n_repl - 1) / (private$n_repl - n_params - 1)
        ) %>%
        dplyr::ungroup()

      # Done with parallel processing.
      future::plan(future::sequential)

      # Generating plots of standardized residuals for each candidate model.
      private$r_resid_plot <- private$plot_models_resid("r")
      private$c_resid_plot <- private$plot_models_resid("c")

      invisible(self)
    },

    #' @description
    #' Perform a Monte Carlo simulation on the selected candidate model to
    #' propagate measurement uncertainty through the analysis function using
    #' the replicate MC samples of the calibration object. Assumes the certified
    #' amount fractions are normally distributed.
    #'
    #' @param model_num (`numeric(1)`)
    #'   Number identifying the candidate polynomial model to use in the
    #'   simulation. Defaults to 1.
    #' @param ... Additional options for the optimization routine.
    #'
    #' @return The calling object invisibly.
    simulate = function(model_num = defaults$MODEL_NUM, ...) {
      private$selected_model <- model_num

      # Extracting necessary parameters from the selected model.
      model_info <- private$models %>%
        dplyr::filter(model_num == private$selected_model)

      solution <- model_info$opt[[1]]$solution
      n_coef <- length(model_info$beta[[1]])
      K <- ncol(private$r_mc)

      # Generating random amount fraction samples from a normal distribution.
      private$c_mc <- stats::rnorm(
        private$n_std * K,
        mean = private$init_vals$c,
        sd = private$init_vals$uc
      ) %>%
        matrix(ncol = K, dimnames = list(private$counts$sample_id))

      # Generating random dark uncertainty samples from a normal distribution.
      private$lambda_mc <- stats::rnorm(K, sd = model_info$tau)

      # Holding all params except beta constant at their optimized values.
      lb <- c(solution[1:(2 * private$n_std + 2)], rep(-Inf, n_coef))
      ub <- c(solution[1:(2 * private$n_std + 2)], rep(Inf, n_coef))

      # Setting up for parallel processing.
      n_workers <- future::availableCores() - 2
      chunk_size <- K / n_workers
      col_chunks <- split(seq_len(K), ceiling(seq_len(K) / chunk_size))

      future::plan(future::multisession, workers = n_workers)

      # Running the optimization algorithm for each set of ratio and amount
      # fraction random samples. Each worker will get a certain chunk of columns
      # to optimize in serial.
      private$solution_mc <- furrr::future_map2(
        lapply(col_chunks, \(cols) private$r_mc[, cols]),
        lapply(col_chunks, \(cols) private$c_mc[, cols]),
        \(r_mc, c_mc) {
          lapply(seq_len(ncol(r_mc)), \(i) {
            private$optimizer(
              solution,
              lb,
              ub,
              eval(str2expression(model_info$funcs)),
              r_mc[, i],
              c_mc[, i],
              ...
            )$solution[-(1:(2 * private$n_std + 2))]
          })
        }
      ) %>%
        unlist() %>%
        matrix(ncol = K)

      # Done with parallel processing.
      future::plan(future::sequential)

      # Estimating residuals and uncertainties from simulation.
      r_resid <- sweep(
        private$r_mc,
        1,
        rep(model_info$rho[[1]], private$counts$n)
      )
      c_resid <- private$c_mc - self$predict(model_info$rho[[1]])

      # Summarizing residuals for both ratios and amounts fractions.
      private$mc_summ <- list(r = r_resid, c = c_resid) %>%
        purrr::map(\(mat) {
          mat %>%
            tibble::as_tibble(
              .name_repair = "unique_quiet",
              rownames = "sample_id"
            ) %>%
            tidyr::pivot_longer(tidyselect::starts_with("...")) %>%
            dplyr::summarise(
              resid = robustbase::huberM(value)$mu,
              uresid = robustbase::Qn(value),
              .by = sample_id
            )
        }) %>%
        dplyr::bind_rows(.id = "source") %>%
        tidyr::pivot_wider(
          names_from = source,
          values_from = c(resid, uresid),
          names_glue = "{source}_{.value}"
        ) %>%
        dplyr::mutate(
          r_zscore = r_resid / r_uresid,
          c_zscore = c_resid / c_uresid
        ) %>%
        dplyr::relocate(
          tidyselect::starts_with("c"),
          .after = tidyselect::last_col()
        )

      # Generating summary plots of residuals.
      private$r_mc_plot <- private$plot_mc_resid("r", r_resid)
      private$c_mc_plot <- private$plot_mc_resid("c", c_resid)

      # Getting covariance of coefficients.
      private$beta_cov <- private$calc_beta_cov(
        model_info$beta[[1]],
        private$solution_mc
      )

      invisible(self)
    },

    #' @description
    #' Predict unknowns utilizing the results obtained from
    #' [`simulate()`](#method-simulate). MC samples are returned and must be
    #' summarized by the caller.
    #'
    #' @param data (NULL | `numeric` | `matrix`)
    #'   Behavior will change based on the value provided:
    #'
    #'   * NULL: will calculate the value for a control sample by taking the
    #'   sums of each set of coefficients. Returns a vector.
    #'   * vector: will calculate the outer product of the vector with the
    #'   coefficients for each power of r. Returns a matrix.
    #'   * matrix: will calculate the product of each element of the matrix
    #'   with the corresponding coefficients for each power of r. Returns a
    #'   matrix.
    #'
    #' @return A `vector` or `matrix` with the calculated MC samples.
    predict = function(data = NULL) {
      if (is.null(private$solution_mc)) {
        self$simulate()
      }

      if (is.null(data)) {
        # Control sample is a ratio of 1, so just taking the sums of each pair
        # of coefficients.
        colSums(private$solution_mc) + private$lambda_mc
      } else {
        # Extracting the powers of r from the analysis function.
        powers <- eval(str2expression(stringr::str_extract(
          private$models$funcs[
            private$models$model_num == private$selected_model
          ],
          "c\\(.+\\)(?=\\))"
        )))

        seq_along(powers) %>%
          purrr::map(\(i) {
            if (is.matrix(data)) {
              sweep(data^powers[i], 2, private$solution_mc[i, ], FUN = "*")
            } else {
              outer(data^powers[i], private$solution_mc[i, ])
            }
          }) %>%
          purrr::reduce(`+`) %>%
          sweep(2, private$lambda_mc, FUN = "+")
      }
    },

    #' @description
    #' Get all private attributes that contain general information on the
    #' analysis. If all attributes are not needed at the same time, it can be
    #' convenient to select a specific one from the function call, i.e.
    #' `get_info()$n_std`.
    #'
    #' @return A `list` of attributes.
    get_info = function() {
      list(
        n_std = private$n_std,
        n_repl = private$n_repl,
        selected_model = private$selected_model
      )
    },

    #' @description
    #' Get all private attributes containing MC samples. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_mc()$r_mc`.
    #'
    #' @return A `list` of attributes.
    get_mc = function() {
      rlang::list2(!!!super$get_mc(), lambda_mc = private$lambda_mc)
    },

    #' @description
    #' Get all private attributes that are data frames. If all attributes are
    #' not needed at the same time, it can be convenient to select a specific
    #' one from the function call, i.e. `get_data()$cal_data`.
    #'
    #' @return A `list` of attributes.
    get_data = function() {
      rlang::list2(
        !!!super$get_data(),
        init_vals = private$init_vals,
        counts = private$counts
      )
    }
  ),
  private = list(
    # `numeric(1)`: Number of calibration standards.
    n_std = NULL,

    # `numeric(1)`: Total number of replicate response ratios.
    n_repl = NULL,

    # `numeric(1)`: Number identifying the model selected in the uncertainty
    # simulation.
    selected_model = NULL,

    # `matrix`: Monte Carlo samples for the estimated dark uncertainty.
    lambda_mc = NULL,

    # [tibble::tibble]: Initial values to use for model setup.
    init_vals = NULL,

    # [tibble::tibble]: Count of the number of replicates for each sample.
    counts = NULL,

    # This function serves as the objective during optimization and represents
    # the likelihood of observing the data given the measurement model, the
    # calibration relationship, and the values of the model's adjustable
    # parameters.
    objective = function(theta, func, x, y) {
      # Extracting parameters from `theta`. `rho` is the true ratio means,
      # `sigma` is the true ratio uncertainty not accounted for by the replicate
      # uncertainties, `nu` is the degrees of freedom for the measurement
      # process, `tau` is the dark uncertainty lurking in the standards, and
      # `beta` is the polynomial coefficients.
      rho <- theta[1:private$n_std]
      sigma <- exp(theta[(private$n_std + 1):(2 * private$n_std)])
      nu <- theta[(2 * private$n_std + 1)]
      tau <- exp(theta[(2 * private$n_std + 2)])
      beta <- theta[-(1:(2 * private$n_std + 2))]

      xi <- vapply(rho, \(r) func(r, beta), numeric(1))
      ur_total <- sqrt(
        private$cal_data$ur^2 + rep(sigma^2 * nu / (nu - 2), private$counts$n)
      )
      uc_total <- sqrt(private$init_vals$uc^2 + tau^2)

      # Likelihood of observing the replicate ratios.
      loglik_r <- stats::dnorm(
        x,
        mean = rep(rho, private$counts$n),
        sd = ur_total,
        log = TRUE
      )

      # Likelihood of observing the certified amount fractions.
      loglik_c <- stats::dnorm(y, mean = xi, sd = uc_total, log = TRUE)

      -sum(loglik_r, loglik_c)
    },

    # This function runs the optimization algorithm for the specified model.
    run_opt = function(form, func) {
      lm0 <- stats::lm(
        stats::formula(form),
        tibble::tibble(r = private$init_vals$rho0, c = private$init_vals$c)
      )
      beta0 <- unname(lm0$coefficients)
      tau0 <- stats::sigma(lm0)
      nu0 <- min(private$n_repl - private$n_std, 100)
      theta0 <- c(
        private$init_vals$rho0,
        log(private$init_vals$sigma0),
        nu0,
        log(tau0),
        beta0
      )

      lb <- c(
        private$init_vals$rho_lb,
        rep(-Inf, private$n_std),
        3,
        rep(-Inf, 1 + length(beta0))
      )
      ub <- c(
        private$init_vals$rho_ub,
        rep(Inf, private$n_std),
        100,
        rep(Inf, 1 + length(beta0))
      )

      private$optimizer(
        theta0,
        lb,
        ub,
        eval(str2expression(func)),
        private$cal_data$r,
        private$init_vals$c,
        maxtime = 300
      )
    }
  )
)
