make_data <- function(n, p = 4, 
                      param_src = .8, 
                      param_cate = .8, 
                      setting_trt = 1,
                      err_sd = .5)
{
   
   x_fn <- function(n) {
      matrix(4 * runif(n * p) - 2, n, p)
   }
   
   prob_src_fn <- function(x) {
      p <- pnorm(1 * x[,2] - 1.2 * x[,1])
      .5 + param_src * (p - .5)
   }
   
   eff_cate_fn <- function(x) {
      tau1 <- pnorm(1.5 * x[,2] + .8 * x[,1] - .4 * (x[,1] - x[,2])^2 - .3) - (x[,1] - x[,2])^2 / 14
      tau2 <- pnorm(.4 * x[,2] + .6 * x[,1]) - .5
      param_cate * tau1 + (1 - param_cate) * tau2
   }
   
   prob_trt_fn <- function(x) {
      p <- switch(
         setting_trt,
         `1` = pnorm(0.5 * x[,1] + 0.3 * x[,2] - .3),
         `2` = pnorm(1.6 * x[,1] + 1.3 * x[,2] - .8),
         `3` = pnorm(.4 * x[,1]^2 + .4 * x[,2]^2 + .5 * x[,1] * x[,2] - .4 * x[,1] + .4 * x[,2] - .9)
      )
      .5 + .8 * (p - .5)
   }
   
   eff_main_fn <- function(x) pnorm(-.6 * x[,1] - .6 * x[,2] + .2 * x[,3] + .5) + .5
   
   err_fn <- function(n) rnorm(n, 0, err_sd)
   
   .make_data_inner(n, x_fn, err_fn, prob_src_fn, prob_trt_fn, eff_main_fn, eff_cate_fn)
   
}

.make_data_inner <- function(n, x_fn, err_fn, prob_src_fn, prob_trt_fn, eff_main_fn, eff_cate_fn)
{
   x <- x_fn(n)
   prob_src <- prob_src_fn(x)
   prob_trt <- prob_trt_fn(x)
   eff_main <- eff_main_fn(x)
   eff_cate <- eff_cate_fn(x)
   err <- err_fn(n)
   
   s <- runif(n) <= prob_src
   trt <- runif(n) <= prob_trt
   
   y0 <- eff_main - eff_cate / 2
   y1 <- eff_main + eff_cate / 2
   
   y <- s * (trt * y1 + (1 - trt) * y0 + err)
   
   list(x   = x,
        s   = as.integer(s),
        trt = s * trt,
        y   = y,
        prob_src = prob_src,
        prob_trt = prob_trt,
        y0  = y0,
        y1  = y1)
}
