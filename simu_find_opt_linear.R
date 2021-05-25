# Find the optimal linear ITR for each simulation setting by grid search, using
# the true mean potential outcomes of a large test sample.

source("simu_make_data.R")
Rcpp::sourceCpp("fit_linear_2dgridsearch.cpp")

n_test <- 1e5
p <- 4

param_cates <- c(0, 0.4)

opt <- list()

for (k in 1:2) {
   name <- sprintf("cate%d", k)
   print(name)
   set.seed(2)
   dat <- make_data(n_test, p, param_cate = param_cates[k])
   
   opt[[name]] <- list(beta = numeric(p),
                       beta_s = numeric(p))
   
   print(system.time({
      xt <- dat$x[dat$s == 0, 1:2]
      y1t <- dat$y1[dat$s == 0]
      y0t <- dat$y0[dat$s == 0]
      coef <- fit_linear_2dgridsearch(xt, y1t - y0t, by = .001)
      opt[[name]]$beta[1:3] <- coef
      d <- drop(xt %*% coef[-1]) + coef[1] > 0
      opt[[name]]$val <- mean(d * y1t + (1 - d) * y0t)
   }))
   print(system.time({
      xs <- dat$x[dat$s == 1, 1:2]
      y1s <- dat$y1[dat$s == 1]
      y0s <- dat$y0[dat$s == 1]
      opt[[name]]$beta_s[1:3] <- fit_linear_2dgridsearch(xs, y1s - y0s, by = .001)
   }))
}

saveRDS(opt, "saved_results/simu_opt.rds")
