estimate_R_adj_inc <- function(incidence_data, P, si_shape, si_rate, 
                               si_t_max, start_day, window_length) {
  ## load EpiEstim and supress the error message for having negative incidence
  library(EpiEstim)

  process_I <- function(incid) {
    # If the input is an incidence object, we want to convert it to a data frame
    # that EpiEstim understands, which contains a single column for the I counts.
    if (inherits(incid, "incidence")) {
      I_inc <- incid
      incid <- as.data.frame(I_inc)
      incid$I <- rowSums(incidence::get_counts(I_inc))
    }
    vector_I <- FALSE
    single_col_df_I <- FALSE
    if (is.vector(incid)) {
      vector_I <- TRUE
    } else if (is.data.frame(incid)) {
      if (ncol(incid) == 1) {
        single_col_df_I <- TRUE
      }
    }
    if (vector_I | single_col_df_I) {
      if (single_col_df_I) {
        I_tmp <- incid[[1]]
      } else {
        I_tmp <- incid
      }
      incid <- data.frame(local = I_tmp, imported = rep(0, length(I_tmp)))
      I_init <- sum(incid[1, ])
      incid[1, ] <- c(0, I_init)
    } else {
      if (!is.data.frame(incid) |
        (!("I" %in% names(incid)) &
          !all(c("local", "imported") %in% names(incid)))) {
        stop("incid must be a vector or a dataframe with either i) a column 
           called 'I', or ii) 2 columns called 'local' and 'imported'.")
      }
      if (("I" %in% names(incid)) &
        !all(c("local", "imported") %in% names(incid))) {
        incid$local <- incid$I
        incid$local[1] <- 0
        incid$imported <- c(incid$I[1], rep(0, nrow(incid) - 1))
      }
      if (incid$local[1] > 0) {
        warning("incid$local[1] is >0 but must be 0, as all cases on the first 
              time step are assumed imported. This is corrected automatically 
              by cases being transferred to incid$imported.")
        I_init <- sum(incid[1, c("local", "imported")])
        incid[1, c("local", "imported")] <- c(0, I_init)
      }
    }

    incid[which(is.na(incid))] <- 0
    date_col <- names(incid) == "dates"
    # if (any(date_col)) {
    #   if (any(incid[, !date_col] < 0)) {
    #     stop("incid must contain only non negative integer values.")
    #   }
    # } else {
    #   if (any(incid < 0)) {
    #     stop("incid must contain only non negative integer values.")
    #   }
    # }

    return(incid)
  }

  ## assignInNamespace("process_I", process_I, ns = "EpiEstim") 
  ## overwrite the process_I function in EpiEstim pacakage to supress 
  ## the error while there are negative values for the incidences, 
  ## need to restart R to restore this action.

  ##
  m <- dim(incidence_data)[2]

  ## create discrete gamma distribution with specified lenght of t with non-zero probability
  w <- sapply(1:si_t_max, function(x) {
    pgamma(x, si_shape, si_rate) - pgamma(x - 1, si_shape, si_rate)
  })

  w <- w / sum(w)

  ## adjust the incidence with P matrix
  if (is.list(P) == T) {
    incidence_data_nomix <- do.call("rbind", lapply(1:nrow(case_matrix), function(x) {
      unlist(case_matrix[x, ]) %*% solve(P[[x]])
    }))
  } else if (is.matrix(P) == T) {
    incidence_data_nomix <- do.call("rbind", lapply(1:nrow(incidence_data), function(x) {
      unlist(incidence_data[x, ]) %*% solve(P)
    }))
  }

  ## set negative incidence to zero for computation considerations
  incidence_data_nomix[incidence_data_nomix < 0] <- exp(incidence_data_nomix[incidence_data_nomix < 0])

  ## obtain imported incidence that is equivalant to the incidence adjustment by P matrix so that the original EpiEstim can run
  import_matrix <- incidence_data[, 1:m] - incidence_data_nomix

  region_names <- names(incidence_data)

  ## estimating R(t) for each region with adjusted incidence
  result_list <- lapply(1:m, function(j) {
    try_data <- data.frame(local = incidence_data_nomix[, j], imported = import_matrix[, j])

    try_data[1, ] <- c(0, try_data[1, 1])


    res_non_parametric_si <- estimate_R(try_data,
      method = "non_parametric_si",
      config = make_config(list(
        t_start = start_day:(dim(try_data)[1] - window_length),
        t_end = (start_day + window_length):(dim(try_data)[1]),
        si_distr = c(0, w)
      ))
    )

    res_non_parametric_si
  })

  names(result_list) <- region_names

  return(result_list)
}