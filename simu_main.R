library(dplyr)
library(foreach)
library(doParallel)
registerDoParallel()

source("utils.R")
source("method.R")
source("simu_make_data.R")

n <- 1600
p <- 4
n_simu <- 4
param_cates <- c(0, 0.4)

param_cate <- param_cates[2]
setting_trt <- 2

system.time(
   result <- foreach(1:n_simu) %do% {
      dat <- make_data(n, p, 
                       param_cate = param_cate,
                       setting_trt = setting_trt)
      with(
         dat,
         {
            xs <- x[s == 1,]
            xt <- x[s != 1,]
            trts <- trt[s == 1]
            ys <- y[s == 1]
            
            prob_trt_pred <- get_prop_score(xs, trts)$s
            prob_src_pred <- get_prob_source(xs, xt)$s
            
            kbal_pos <- kbal3way(xs, trts, ys, xt, type = "gaussian", echo = TRUE)
            kbal <- kbal3way(xs, trts, ys, xt, type = "gaussian", echo = TRUE, keep_w_positive = FALSE)
            
            wts <- 
               tibble(
                  kernel = kbal$weights,
                  kernel_pos = kbal_pos$weights,
                  ebal_s = ebal_generalize(xs, trts, xs, degree = 1),
                  ebal_t = ebal_generalize(xs, trts, xt, degree = 1),
                  prob_s = trts / prob_trt_pred + (1 - trts) / (1 - prob_trt_pred),
                  true_s = trts / prob_trt[s == 1] + (1 - trts) / (1 - prob_trt[s == 1]),
                  prob_t = prob_s / prob_src_pred * (1 - prob_src_pred),
                  true_t = true_s / prob_src[s == 1] * (1 - prob_src[s == 1]),
                  prob_o = trts * (1 - prob_trt_pred) + (1 - trts) * prob_trt_pred,
                  true_o = trts * (1 - prob_trt[s == 1]) + (1 - trts) * prob_trt[s == 1]
               )
            
            wts <- wts %>%
               mutate(trts = trts) %>%
               group_by(trts) %>%
               mutate_at(vars(-group_cols()), function(z) z / sum(z)) %>%
               ungroup() %>%
               select(-trts)
            
            betas <- lapply(wts, function(w) fit_linear_itr_owl(xs, trts, ys, w)$beta)
            
            list(wts = wts,
                 kbal_param = kbal$param,
                 kbal_pos_param = kbal_pos$param,
                 betas = betas,
                 cor_all = cor(wts),
                 cor_kbal = sapply(1:NCOL(kbal$tune_res$w_mat), 
                                   function(j) cor(kbal$tune_res$w_mat[,j], kbal_pos$tune_res$w_mat[,j])))
         }
      )
   }
)
