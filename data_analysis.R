#' This code runs the data analysis for a single round under one setting. Figure
#' 3 in the paper is based on 400 rounds under each of 6 settings.

library(foreach)
library(dplyr)

source("method.R")
source("utils.R")

p_trt_case <- 1 # choice of the g function: 0 or 1 or 2, corresponding to (a), (b), (c)
p_src_scale <- 1 # choice of the f function: 0 or 1, corresponding to (a), (b)

src_coefs <- p_src_scale * c(age = 1,
                             weight = -1,
                             genderF = .6)


# prepare data ------------------------------------------------------------

load("echo_data.rdata")

dat <- dat %>%
   mutate_at(vars(starts_with("lab_") & !ends_with("_flag")), function(v) log(1 + v)) %>% 
   mutate_at(all_of(features_cont), function(fea) {fea <- scale(fea); attributes(fea) <- NULL; fea})

fml_subgrp <- c(features_cont, features_disc) %>% paste(collapse = " + ") %>% sprintf(" ~ %s", .)
x <- model.matrix(as.formula(fml_subgrp), data = dat)[, -1]
probs <- pnorm(drop(x[, names(src_coefs)] %*% src_coefs) - .4 * p_src_scale)
probs <- .5 + .8 * (probs - .5)
s_idx <- sample(1:NROW(x), round(.4 * NROW(x)), prob = probs)
s <- rep(0, NROW(x))
s[s_idx] <- 1
trt <- as.integer(dat$echo == 1)
y <- dat$mort_28_day == 0

# split into source, target, test (1:1:2)

xs <- x[s == 1,]
trts <- trt[s == 1]
ys <- y[s == 1]
xt <- x[s != 1,]
trtt <- trt[s != 1]
yt <- y[s != 1]

trt_idx1 <- which(trts == 1)
trt_idx0 <- which(trts == 0)
if (p_trt_case == 0) {
   fz <- rep(0, NROW(xs))
}
if (p_trt_case == 1) {
   z <- xs[, c("saps", "sofa", "elix_score")]
   fz <- drop(1 * z[, 1] + .8 * z[,2] - .9 * z[,3])
}
if (p_trt_case == 2) {
   z <- xs[, c("saps", "sofa", "elix_score")]
   fz <- .3 * z[,1]^2 + .2 * z[,2]^2 + .2 * z[,3]^2 + .4 * z[,1] * z[,2] - .3 * z[,1] + .4* z[,2] + .2 * z[,3] - .6
}
prob1 <- pnorm(fz[trt_idx1])
prob0 <- pnorm(-fz[trt_idx0])
prob1 <- .5 + .8 * (prob1 - .5)
prob0 <- .5 + .8 * (prob0 - .5)
sub_idx1 <- sample(trt_idx1, round(.5 * length(trt_idx1)), prob = prob1)
sub_idx0 <- sample(trt_idx0, round(.5 * length(trt_idx0)), prob = prob0)
sub_idx <- sort(c(sub_idx1, sub_idx0))
xs <- xs[sub_idx, ]
trts <- trts[sub_idx]
ys <- ys[sub_idx]

in_test <- sample(c(rep(TRUE, round(NROW(xt) * 2/3)), rep(FALSE, NROW(xt) - round(NROW(xt) * 2/3))))
x_test <- xt[in_test,]
trt_test <- trtt[in_test]
y_test <- yt[in_test]
xt <- xt[!in_test,]
trtt <- trtt[!in_test]
yt <- yt[!in_test]


# obtain ITRs -------------------------------------------------------------

prob_trt <- get_prop_score(xs, trts)$s
prob_src <- get_prob_source(xs, xt)$s

system.time({
   kbal_kernel <- kbal3way(xs, trts, ys, xt, type = "gaussian", echo = TRUE)
})

wts <- tibble(kernel = kbal_kernel$weights,
              ebal_s = ebal_generalize(xs, trts, xs),
              ebal_t = ebal_generalize(xs, trts, xt),
              prob_s = trts / prob_trt + (1 - trts) / (1 - prob_trt),
              prob_t = (1 - prob_src) / prob_src * prob_s,
              prob_o = trts * (1 - prob_trt) + (1 - trts) * prob_trt) %>% 
   mutate(trts = trts) %>%
   group_by(trts) %>%
   mutate_at(vars(-group_cols()), function(z) z / sum(z)) %>%
   ungroup() %>%
   select(-trts)

itrs <- lapply(wts, function(w) fit_linear_itr_owl(xs, trts, ys, weights = w))
d_test <- lapply(itrs, function(itr) itr$pred_fn(x_test) > 0)
d_test$treat_all <- rep(TRUE, NROW(x_test))


# output ------------------------------------------------------------------

betas <- 
   lapply(itrs,
          function(itr) itr$beta) %>% 
   do.call(what = rbind)
colnames(betas) <- c("Intercept", colnames(x))

ps_test <- get_prop_score(x_test, trt_test)$s
wt_test <- trt_test / ps_test + (1 - trt_test) / (1 - ps_test)
wt_test[trt_test == 1] <- wt_test[trt_test == 1] / sum(wt_test[trt_test == 1]) * length(trt_test)
wt_test[trt_test != 1] <- wt_test[trt_test != 1] / sum(wt_test[trt_test != 1]) * length(trt_test)
eval_test <- 
   lapply(d_test,
          function(d) {
             mean(wt_test * (2 * trt_test - 1) * d * y_test)
          }) %>% 
   do.call(what = rbind)

eval_test
