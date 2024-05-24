library(rstan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(gridExtra)
library(parallel)
library(data.table)

MA_cases_count <- readRDS("MA_all_counties_cases_count.rds")

avg_mobility_matrix_list <- readRDS("avg_mobility_matrix_ma_all_counties_list.rds")

p <- sqrt(length(avg_mobility_matrix_list[[1]]))


geoid_link <- fread("all-geocodes-v2017.csv")
geoid_link_mass <- geoid_link[geoid_link$`State Code (FIPS)` == 25, ]
ma_fips <- seq(25001, 25027, 2)

geoid_link_mass_unique <- geoid_link_mass[!duplicated(geoid_link_mass$`County Code (FIPS)`), ] %>% as.data.frame()
county_names <- geoid_link_mass_unique$`Area Name (including legal/statistical area description)`[-1]
county_names <- gsub(" County", "", county_names)



geoids <- ma_fips

case_count <- data.frame(date = rep(seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day"), each = p), geoid = rep(geoids, length(seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day"))))

MA_cases_count$cdc_report_dt <- as.Date(MA_cases_count$cdc_report_dt)

case_count <- dplyr::left_join(case_count, MA_cases_count, by = c("date" = "cdc_report_dt", "geoid" = "county_fips_code"))

case_count$n[is.na(case_count$n)] <- 0

case_count$county <- factor(case_count$geoid, levels = geoids, labels = county_names)

library(tidyr)

MA_cases_count_wide <- case_count %>%
  select(., date, geoid, n) %>%
  spread(., key = geoid, value = n)

MA_cases_count_wide_sub <- MA_cases_count_wide[MA_cases_count_wide$date >= "2020-07-01", ]

case_matrix <- as.matrix(MA_cases_count_wide_sub[, -1])

P_list <- lapply(avg_mobility_matrix_list[1:length(seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day"))], function(x) {
  matrix(x, p)
})

stan_data <- list()

stan_data$cases <- case_matrix

stan_data$K <- dim(stan_data$cases)[2]

stan_data$N <- dim(stan_data$cases)[1]

w <- sapply(1:nrow(stan_data$cases), function(x) {
  pgamma(x, 3.45, 0.66) - pgamma(x - 1, 3.45, 0.66)
})

w <- w / sum(w)
w <- c(w, rep(0, nrow(stan_data$cases) - length(w)))

stan_data$SI <- w


## Load R files in the EpiEstim package
for (i in dir("EpiEstim/R/", full.names = T)) {
  source(i)
}

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

P_list <- lapply(avg_mobility_matrix_list[1:length(seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day"))], function(x) {
  matrix(x, p)
})


case_matrix_nomix <- do.call("rbind", lapply(1:nrow(case_matrix), function(x) {
  unlist(case_matrix[x, ]) %*% solve(P_list[[x]])
}))


shrink_scale <- 1
round_num <- 0
while (sum(sapply(order(apply(case_matrix, 2, mean)), function(x) {
  sum(round(case_matrix_nomix[, x]) < 0)
}) > 0)) {
  print(round_num <- round_num + 1)

  for (k in order(apply(case_matrix, 2, mean))) {
    for (d in 1:212) {
      print(c(k, d))
      shrink_scale <- 1

      while (case_matrix_nomix[d, k] < 0) {
        shrink_scale <- max(c(0, shrink_scale - 0.01))

        for (i in c(1:14)[-k]) {
          P_list[[d]][i, i] <- P_list[[d]][i, i] + P_list[[d]][k, i] * (1 - shrink_scale)

          P_list[[d]][k, i] <- P_list[[d]][k, i] * shrink_scale

          P_list[[d]][k, k] <- P_list[[d]][k, k] + sum(P_list[[d]][-k, k] * (1 - shrink_scale))

          P_list[[d]][-k, k] <- P_list[[d]][-k, k] * shrink_scale
        }

        case_matrix_nomix <- do.call("rbind", lapply(1:nrow(case_matrix), function(x) {
          unlist(case_matrix[x, ]) %*% solve(P_list[[x]])
        }))
      }
    }
  }
}


case_matrix_nomix <- do.call("rbind", lapply(1:nrow(case_matrix), function(x) {
  unlist(case_matrix[x, ]) %*% solve(P_list[[x]])
}))


import_matrix <- case_matrix - case_matrix_nomix

r_scale_result <- lapply(1:ncol(case_matrix), function(j) {
  try_data <- data.frame(local = case_matrix_nomix[, j], imported = import_matrix[, j])

  try_data[1, 2] <- sum(try_data[1, ])
  try_data[1, 1] <- 0
  res_non_parametric_si <- estimate_R(try_data,
    method = "non_parametric_si",
    config = make_config(list(
      si_distr = c(0, w)
    ))
  )

  res_non_parametric_si
})

geoid_link <- fread("all-geocodes-v2017.csv")
geoid_link_mass <- geoid_link[geoid_link$`State Code (FIPS)` == 25, ]
ma_fips <- seq(25001, 25027, 2)

geoid_link_mass_unique <- geoid_link_mass[!duplicated(geoid_link_mass$`County Code (FIPS)`), ] %>% as.data.frame()
county_names <- geoid_link_mass_unique$`Area Name (including legal/statistical area description)`[-1]
county_names <- gsub(" County", "", county_names)




dat_regions <- do.call("rbind", lapply(1:ncol(case_matrix), function(x) {
  dat <- as.data.frame(r_scale_result[[x]]$R[, c(3, 5, 11)])
  dat$region <- county_names[x]
  dat$day <- 1:nrow(dat)
  dat
}))

dat_regions_R_est <- spread(dat_regions[, c(1, 4, 5)], region, `Mean(R)`)[, -1]
dat_regions_R_lb_est <- spread(dat_regions[, c(2, 4, 5)], region, `Quantile.0.025(R)`)[, -1]
dat_regions_R_ub_est <- spread(dat_regions[, c(3, 4, 5)], region, `Quantile.0.975(R)`)[, -1]
dat_regions_inc <- as.matrix(case_matrix)



estimatated_incidence <- sapply((1:nrow(dat_regions_R_est)) + 7, function(t) {
  rev(w[1:min(length(w), (t - 1))]) %*% dat_regions_inc[1:(t - 1), , drop = F] %*% diag(dat_regions_R_est[t - 7, ]) %*% P_list[[t]]
})




estimatated_incidence_lb <- sapply((1:nrow(dat_regions_R_est)) + 7, function(t) {
  rev(w[1:min(length(w), (t - 1))]) %*% dat_regions_inc[1:(t - 1), , drop = F] %*% diag(dat_regions_R_lb_est[t - 7, ]) %*% P_list[[t]]
})



estimatated_incidence_ub <- sapply((1:nrow(dat_regions_R_est)) + 7, function(t) {
  rev(w[1:min(length(w), (t - 1))]) %*% dat_regions_inc[1:(t - 1), , drop = F] %*% diag(dat_regions_R_ub_est[t - 7, ]) %*% P_list[[t]]
})


included_days <- seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day")

estimatated_incidence <- as.data.frame(t(estimatated_incidence))





names(estimatated_incidence) <- county_names

estimatated_incidence_long <- gather(estimatated_incidence)
estimatated_incidence_long$day <- rep(included_days[8:212], p)

estimatated_incidence_long$yl <- gather(as.data.frame(t(estimatated_incidence_lb)))$value

estimatated_incidence_long$yh <- gather(as.data.frame(t(estimatated_incidence_ub)))$value


### scaled_nomob


r_scale_result_nomob <- lapply(1:ncol(case_matrix), function(j) {
  try_data <- data.frame(local = case_matrix[, j], imported = 0)

  try_data[1, 2] <- sum(try_data[1, ])
  try_data[1, 1] <- 0
  res_non_parametric_si <- estimate_R(try_data,
    method = "non_parametric_si",
    config = make_config(list(
      si_distr = c(0, w)
    ))
  )

  res_non_parametric_si
})


dat_regions_nomob <- do.call("rbind", lapply(1:ncol(case_matrix), function(x) {
  dat <- as.data.frame(r_scale_result_nomob[[x]]$R[, c(3, 5, 11)])
  dat$region <- county_names[x]
  dat$day <- 1:nrow(dat)
  dat
}))

dat_regions_R_est_nomob <- spread(dat_regions_nomob[, c(1, 4, 5)], region, `Mean(R)`)[, -1]
dat_regions_R_lb_est_nomob <- spread(dat_regions_nomob[, c(2, 4, 5)], region, `Quantile.0.025(R)`)[, -1]
dat_regions_R_ub_est_nomob <- spread(dat_regions_nomob[, c(3, 4, 5)], region, `Quantile.0.975(R)`)[, -1]
dat_regions_inc <- as.matrix(case_matrix)


estimatated_incidence_nomob <- sapply((1:nrow(dat_regions_R_est_nomob)) + 7, function(t) {
  rev(w[1:min(length(w), (t - 1))]) %*% dat_regions_inc[1:(t - 1), , drop = F] %*% diag(dat_regions_R_est_nomob[t - 7, ])
})




estimatated_incidence_lb_nomob <- sapply((1:nrow(dat_regions_R_est_nomob)) + 7, function(t) {
  rev(w[1:min(length(w), (t - 1))]) %*% dat_regions_inc[1:(t - 1), , drop = F] %*% diag(dat_regions_R_lb_est_nomob[t - 7, ])
})



estimatated_incidence_ub_nomob <- sapply((1:nrow(dat_regions_R_est_nomob)) + 7, function(t) {
  rev(w[1:min(length(w), (t - 1))]) %*% dat_regions_inc[1:(t - 1), , drop = F] %*% diag(dat_regions_R_ub_est_nomob[t - 7, ])
})



included_days <- seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day")

estimatated_incidence_nomob <- as.data.frame(t(estimatated_incidence_nomob))





names(estimatated_incidence_nomob) <- county_names

estimatated_incidence_long_nomob <- gather(estimatated_incidence_nomob)
estimatated_incidence_long_nomob$day <- rep(included_days[8:212], p)

estimatated_incidence_long_nomob$yl <- gather(as.data.frame(t(estimatated_incidence_lb_nomob)))$value

estimatated_incidence_long_nomob$yh <- gather(as.data.frame(t(estimatated_incidence_ub_nomob)))$value


### inc_renewal_equ

data <- readRDS(paste0("ma_all_counties_adjusted_P_out_extracted.rds"))

data$County <- factor(data$County, levels = geoids, labels = county_names)

date_range <- seq(as.Date("2020-07-08"), as.Date("2021-01-28"), by = "day")


names(estimatated_incidence_long)[c(1:2)] <- c("region", "y")

inc_scale_plot_data <- estimatated_incidence_long

inc_scale_plot_data <- left_join(inc_scale_plot_data, dat_plot, by = c("region", "day"))

inc_scale_plot_data <- na.omit(inc_scale_plot_data)

inc_scale_plot_data <- inc_scale_plot_data[inc_scale_plot_data$day %in% date_range, ]

inc_renewal_equ_plot_data <- data[data$x %in% date_range, ]

names(inc_renewal_equ_plot_data) <- c("day", "y_real", "y", "yl", "yh", "r", "rl", "ru", "region")

real_inc_data <- inc_renewal_equ_plot_data[, c(1:2, 9)]

inc_renewal_equ_plot_data <- inc_renewal_equ_plot_data[inc_renewal_equ_plot_data$day %in% inc_scale_plot_data$day, ]

inc_renewal_equ_plot_data <- inc_renewal_equ_plot_data[, -2]


#### inc_renewal_equ_nomob


data <- readRDS(paste0("ma_all_counties_adjusted_P_no_mobility_out_extracted.rds"))


data$County <- factor(data$County, levels = geoids, labels = county_names)


date_range <- seq(as.Date("2020-07-08"), as.Date("2021-01-28"), by = "day")


names(estimatated_incidence_long_nomob)[c(1:2)] <- c("region", "y")

inc_scale_plot_data_nomob <- estimatated_incidence_long_nomob

inc_scale_plot_data_nomob <- left_join(inc_scale_plot_data_nomob, dat_plot, by = c("region", "day"))

inc_scale_plot_data_nomob <- na.omit(inc_scale_plot_data_nomob)

inc_scale_plot_data_nomob <- inc_scale_plot_data_nomob[inc_scale_plot_data_nomob$day %in% date_range, ]

inc_renewal_equ_plot_data_seperate <- data[data$x %in% date_range, ]

names(inc_renewal_equ_plot_data_seperate) <- c("day", "y_real", "y", "yl", "yh", "r", "rl", "ru", "region")

real_inc_data <- inc_renewal_equ_plot_data_seperate[, c(1:2, 9)]

inc_renewal_equ_plot_data_seperate <- inc_renewal_equ_plot_data_seperate[inc_renewal_equ_plot_data_seperate$day %in% inc_scale_plot_data$day, ]

inc_renewal_equ_plot_data_seperate <- inc_renewal_equ_plot_data_seperate[, -2]



####


inc_renewal_equ_plot_data$Method <- "II"
inc_scale_plot_data$Method <- "I"
inc_scale_plot_data_nomob$Method <- "Without Mobility"

inc_both_method_plot_data <- rbind(rbind(inc_scale_plot_data, inc_renewal_equ_plot_data), inc_scale_plot_data_nomob)

date_range <- seq(as.Date("2020-07-15"), as.Date("2021-01-28"), by = "day")
inc_plot_both <- ggplot(inc_both_method_plot_data[inc_both_method_plot_data$day %in% date_range, ], aes(x = day)) +
  geom_line(inc_both_method_plot_data[inc_both_method_plot_data$day %in% date_range, ], mapping = aes(y = y, x = day, color = Method), alpha = 0.7, size = 0.6) +
  geom_ribbon(inc_both_method_plot_data[inc_both_method_plot_data$day %in% date_range, ], mapping = aes(x = day, ymin = yl, ymax = yh, color = Method, fill = Method), alpha = 0.3, colour = NA) +
  geom_col(real_inc_data[real_inc_data$day %in% date_range, ], mapping = aes(y = y_real, x = day), alpha = 0.7, size = 1.2) +
  facet_wrap(~region, nrow = 2, scales = "free_y") +
  theme_classic() +
  xlab("Days") +
  ylab("Incidence") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.45, hjust = 1.3))


rt_plot_both <- ggplot(inc_both_method_plot_data[inc_both_method_plot_data$day %in% date_range, ], aes(x = day, color = Method, fill = Method)) +
  geom_line(inc_both_method_plot_data[inc_both_method_plot_data$day %in% date_range, ], mapping = aes(y = r, x = day, color = Method), alpha = 0.7, size = 0.6) +
  geom_ribbon(inc_both_method_plot_data[inc_both_method_plot_data$day %in% date_range, ], mapping = aes(x = day, ymin = rl, ymax = ru, color = Method), alpha = 0.3, colour = NA) +
  facet_wrap(~region, nrow = 2, scales = "free_y") +
  theme_classic() +
  xlab("Days") +
  ylab("Reproductive Number") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.45, hjust = 1.3))


############# Add pop change


library(data.table)


library(dplyr)
ma_fips <- seq(25001, 25027, 2)
geoid_link <- fread("all-geocodes-v2017.csv")
geoid_link_mass <- geoid_link[geoid_link$`State Code (FIPS)` == 25, ]


geoid_link_mass_unique <- geoid_link_mass[!duplicated(geoid_link_mass$`County Code (FIPS)`), ] %>% as.data.frame()
county_names <- geoid_link_mass_unique$`Area Name (including legal/statistical area description)`[-1]
county_names <- gsub(" County", "", county_names)

county_ids <- data.frame(
  geoid = ma_fips,
  name = county_names
)


mobility_ma_list <- readRDS("mobility_ma_list.rds")


p <- nrow(county_ids)




library(dplyr)
mobility_ma_list_nonames <- lapply(mobility_ma_list, function(x) {
  mobility_meta <- data.frame(
    geoid_o = rep(county_ids$geoid, each = p),
    geoid_d = rep(county_ids$geoid, p)
  )

  mobility_meta <- left_join(mobility_meta, x, by = c("geoid_o", "geoid_d"))
  names(mobility_meta) <- names(mobility_ma_list[[1]])

  mobility_meta$pop_flows[is.na(mobility_meta$pop_flows)] <- 0
  mobility_meta$visitor_flows[is.na(mobility_meta$visitor_flows)] <- 0
  mobility_meta$date_range[is.na(mobility_meta$date_range)] <- na.omit(unique(mobility_meta$date_range))
  return(mobility_meta)
})



s_p_vec <- do.call(
  "rbind",
  lapply(2:length(mobility_ma_list_nonames), function(x) {
    print(x)
    matrix_p_i <- matrix(mobility_ma_list_nonames[[x]]$pop_flows, p)
    matrix_p_i_1day_earlier <- matrix(mobility_ma_list_nonames[[x - 1]]$pop_flows, p)
    sapply(1:p, function(j) {
      matrix_p_i_a <- matrix_p_i_1day_earlier[j, j] + sum(matrix_p_i_1day_earlier[j, -j])

      matrix_p_i_w <- sum(matrix_p_i[j, -j])

      matrix_p_i_z <- sum(matrix_p_i[-j, j])

      s_p_i <- (-matrix_p_i_z + matrix_p_i_w) / (matrix_p_i_a)

      s_p_i
    })
  })
)
library(tidyr)
s_p_vec <- as.data.frame(s_p_vec)
names(s_p_vec) <- county_ids$name
s_p_vec_long <- gather(s_p_vec, key = "region", value = "pop")

date_range_all <- seq(as.Date("2020-02-02"), as.Date("2021-03-31"), by = "day")

date_range <- seq(as.Date("2020-07-15"), as.Date("2021-01-28"), by = "day")

s_p_vec_long$day <- rep(seq(as.Date("2020-02-02"), as.Date("2021-03-31"), by = "day"), ncol(s_p_vec))

s_p_vec_long_sub <- s_p_vec_long[s_p_vec_long$day %in% date_range, ]

s_p_vec_long_sub_change <- s_p_vec_long_sub

g_scale_r <- ggplot(s_p_vec_long_sub, aes(y = pop, x = day)) +
  geom_line() +
  facet_wrap(~region, nrow = 3) +
  geom_hline(yintercept = 0, color = "black", size = 1, alpha = 0.5) +
  facet_wrap(~region, nrow = 2) +
  theme_classic() +
  xlab("Days") +
  ylab("Population Change") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.45, hjust = 1.3))


library(gridExtra)
library(cowplot)

jpeg("real_data_inc_rt_both_with_pop_change.jpg", res = 300, w = 3670, h = 3300)
ggdraw() +
  draw_plot(g_scale_r, 0, 0.66, 0.93, 0.33) +
  draw_plot(inc_plot_both, 0, 0.33, 1, 0.33) +
  draw_plot(rt_plot_both, 0, 0, 1, 0.33)
dev.off()
