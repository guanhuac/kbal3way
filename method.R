library(kernlab)
library(foreach)
library(xgboost)
library(osqp)


# Kernel balancing --------------------------------------------------------

#' Main function for computing the proposed weights
#'
#' @param x matrix of covariates of the source sample.
#' @param trt vector of treatment indicators of the source sample.
#' @param y vector of observed outcomes of the source sample (higher means
#'   better).
#' @param x_target matrix of covariates for the target sample.
#' @param standardize_x a logical indicating whether to standardize the
#'   covariates, default as TRUE.
#' @param type choose what kernel/distance to use: gaussian, guassian_anova, or
#'   energy.
#' @param keep_w_positive a logical indicating whether the weights should be
#'   constrained to be non-negative.
#' @param alpha a numeric vector of alpha values to choose from. If not
#'   supplied, it will be automatically set as a sequence from 1 to 0.
#' @param alpha_length length of automatically generated alpha sequence.
#' @param alpha_val_tol tolerance level when selecting alpha based on the value
#'   function. 0 means the alpha correponding to the largest value function will
#'   be selected. A positive value allows for choosing other alpha value that
#'   leads to more regularized weights.
#' @param augment_y a logical indicating whether y should be replaced with
#'   residuals.
#' @param lambda lambda value. If not supplied, it will be adaptively determined
#'   for each alpha.
#' @param echo a logical indicating whether to echo result for each alpha.
#' @return A list of
#'   weights: the weights under the selected alpha and lambda;
#'   param: a numeric vector of the selected alpha and lambda;
#'   tune_res: results under other hyperparameters. Specifically, tune_log
#'     records the selected lambda for each alpha, the corresponding MMD between
#'     the treated and controls, and the norm of the corresponding weights;
#'     w_mat records the resulting weights under each alpha.
kbal3way <- function(x, trt, y, x_target = NULL, 
                     standardize_x = TRUE, 
                     type = c("gaussian", "gaussian_anova", "energy"),
                     keep_w_positive = TRUE,
                     alpha, alpha_length = 12, alpha_val_tol = 0,
                     augment_y = FALSE, lambda, echo = FALSE)
{
    
    stopifnot(length(trt) == NROW(x))
    if (!is.null(x_target)) stopifnot(NCOL(x) == NCOL(x_target))
    stopifnot(any(trt == 1) && any(trt != 1))
    trt[trt != 1] <- 0
    type <- match.arg(type)
    
    # Compute the kernel matrix
    n <- NROW(x)
    n_t <- NROW(x_target)
    x_all <- rbind(x, x_target)
    K_mat <- .compute_kernel(x_all, standardize_x, type)
    
    
    ## Tuning procedure
    
    # If alpha is not supplied, set it as a sequence from 1 to 0 (more densely
    # allocated near 0) with length equal to alpha_length
    if (missing(alpha))
    {
        alpha_grid <- c(seq(1, .1, length.out = alpha_length - 1)^4, 0)
    } else 
    {
        alpha_grid <- alpha
    }
    
    # w_mat, tune_log are for recording the tuning results. In tune_log, mmd_01
    # means the MMD between the reweighted treated and control groups, and
    # w_norm is the L2 norm of the weights.
    w_mat <- matrix(0, n, length(alpha_grid))
    tune_log <- matrix(0, length(alpha_grid), 4)
    colnames(tune_log) <- c("alpha", "lambda", "mmd_01", "w_norm")
    tune_log[, 1] <- alpha_grid
    
    # Initialize weights
    w <- ifelse(trt == 1, 1 / sum(trt == 1), 1 / sum(trt != 1))
    
    # Compute weights for each alpha, as well as the corresponding tune_log,
    # using the weights from the previous step as a warm start.
    for (j in seq_along(alpha_grid))
    {
        tmp <- .kbal3way_getw(K_mat, trt, type = type, 
                              alpha = alpha_grid[j], 
                              w_start = w,
                              lambda = lambda,
                              keep_w_positive = keep_w_positive)
        
        if (echo) cat(sprintf("j = %d, alpha = %.4f: lambda = %.4f\n", j, alpha_grid[j], tmp$lambda))
        
        w <- w_mat[, j] <- tmp$weights
        tune_log[j, -1] <- c(tmp$lambda, tmp$mmd_01, tmp$w_norm)
    }
    
    # Select alpha based on the value function on the target sample.
    if (length(alpha_grid) > 1)
    {
        if (augment_y) y <- y - .augment_y(x, trt, y)
        
        pred1 <- .train_predict_boot(x[trt == 1,], y[trt == 1], x_target)
        pred0 <- .train_predict_boot(x[trt != 1,], y[trt != 1], x_target)
        tau <- pred1 - pred0
        
        val_vec <- foreach(j = 1:NCOL(w_mat), .combine = c) %do% {
            pred_fn <- fit_linear_itr_owl(x, trt, y, w_mat[, j])$pred_fn
            d <- 2 * (pred_fn(x_target) > 0) - 1
            mean(tau * d)
        }
        
        tune_log <- cbind(tune_log, val = val_vec)
        
        if (alpha_val_tol == 0) {
            j <- which.max(val_vec)
        } else {
            has_good_val <- (max(val_vec) - val_vec) / abs(max(val_vec)) < alpha_val_tol
            w_norm_tmp <- ifelse(has_good_val, tune_log[, "w_norm"], Inf)
            j <- which.min(w_norm_tmp)
        }
        
    } else j <- 1
    
    # output
    list(weights = w_mat[, j],
         param = c(tune_log[j, 1], tune_log[j, 2]),
         tune_res = list(tune_log = tune_log,
                         w_mat = w_mat))
}

#' Inner function for computing the kernel matrix from the covariate matrix
.compute_kernel <- function(x, standardize_x, type)
{
    
    if (standardize_x) x <- scale(x)
    
    dist_x <- fields::rdist(x)
    
    if (type == "energy") return(-dist_x)
    
    # For gaussian/gaussian_anova kernel, set the bandwidth as median distance.
    median_dist <- median(dist_x)
    
    if (type == "gaussian")
    {
        K_mat <- exp(-0.5 * dist_x ^ 2 / median_dist ^ 2)
    } else if (type == "gaussian_anova")
    {
        K_mat <- kernelMatrix(anovadot(sigma = 0.5 / median_dist ^ 2, degree = 2), x)
    }
    
    K_mat
}

#' Inner function for Computing weights for a given alpha, as well as the
#' corresponding tune_log.
.kbal3way_getw <- function(K_mat, trt,
                           type, alpha, 
                           w_start, keep_w_positive,
                           lambda, lambda_length = 12, B = 50)
{
    
    n <- length(trt)
    n_t <- NROW(K_mat) - n
    
    # we will optimize weighted energy distance by minimizing f(w) = w'Qw + 2a'w
    # subject to constraints on w. constraints are that w >= 0 and Aw = b, where
    # b is a vector.
    
    # The following code constructs the Q matrix and 'a' vector corresponding to
    # the three-way kernel balancing or energy balancing.
    
    Q_x_x <- K_mat[1:n, 1:n]
    Q_x_xt <- K_mat[1:n, n + 1:n_t]
    
    Q_mat <- (trt %*% t(trt) + (1-trt) %*% t(1-trt) - 
                  (1-alpha) * trt %*% t(1-trt) - (1-alpha) * (1-trt) %*% t(trt)) * Q_x_x
    
    a_vec <- -alpha * rowSums(Q_x_xt) / n_t
    
    # Set up constraints on weights.
    A_mat <- matrix(0, 2, n)
    A_mat[1, trt == 1] <- 1
    A_mat[2, trt != 1] <- 1
    sum_constr <- c(1, 1)
    
    # Specify lambda_grid to choose from.
    if (missing(lambda))
    {
        if (type == "energy") {
            lambda_grid <- 0
        } else {
            lambda_grid <- exp(seq(log(1e-3), log(10), length.out = lambda_length))
        }
    } else
    {
        lambda_grid <- lambda
    }
    
    # Set up osqp.
    QP <- osqp(P = Q_mat + lambda_grid[1] * diag(n), q = a_vec, 
               A = rbind(A_mat, diag(n)), 
               l = c(sum_constr, rep(ifelse(keep_w_positive, 0, -Inf), n)), u = c(sum_constr, rep(Inf, n)),
               osqpSettings(verbose = FALSE, eps_rel = 1e-6, eps_abs = 1e-6))
    
    if (!missing(w_start)) QP$WarmStart(x = w_start)
    
    # Compute the weights for each lambda, using the result from the previous step as a warm start.
    w_mat <- foreach(j = 1:lambda_length, .combine = cbind) %do% 
        {
            if (j > 1) QP$Update(Px = Matrix::Matrix(Q_mat + lambda_grid[j] * diag(n), sparse = TRUE)@x)
            w <- QP$Solve()$x
            if (keep_w_positive) w[w < 0] <- 0
            w[trt == 1] <- w[trt == 1] / sum(w[trt == 1])
            w[trt != 1] <- w[trt != 1] / sum(w[trt != 1])
            w
        }
    
    # Select lambda based on the treated-control balance
    if (length(lambda_grid) > 1)
    {
        mmd01_mat <- foreach(1:B, .combine = rbind) %do% 
            {
                # boot_idx <- sample(1:n, replace = TRUE)
                boot_idx <- sample(1:n, floor(0.8 * n), replace = FALSE)
                w_mat_boot <- w_mat[boot_idx,]
                trt_boot <- trt[boot_idx]
                w_mat_boot[trt_boot == 1,] <- sweep(w_mat_boot[trt_boot == 1,], 2, colSums(w_mat_boot[trt_boot == 1,]), "/")
                w_mat_boot[trt_boot != 1,] <- -sweep(w_mat_boot[trt_boot != 1,], 2, colSums(w_mat_boot[trt_boot != 1,]), "/")
                sqrt(colSums(w_mat_boot * (K_mat[boot_idx, boot_idx] %*% w_mat_boot)))
            }
        
        j <- which.min(colMeans(mmd01_mat^2))
    } else j <- 1
    
    
    # Output
    weights <- w_mat[, j]
    lambda <- lambda_grid[j]
    weights2 <- ifelse(trt == 1, weights, -weights)
    
    list(weights = weights,
         lambda = lambda,
         mmd_01 = sqrt(sum(weights2 * (Q_x_x %*% weights2))),
         w_norm = sqrt(sum(weights2^2)))
}

#' Fit xgboost with (x, y) and get predictions on newx
.train_predict_boot <- function(x, y, newx) 
{
    obj <- switch(typeof(y),
                  logical = "binary:logistic",
                  double = "reg:squarederror")
    eval_metric <- switch(typeof(y),
                          logical = "logloss",
                          double = "rmse")
    
    cvfit <- xgb.cv(data = x, label = y,
                    objective = obj, metrics = eval_metric,
                    nfold = 5,
                    nrounds = 60, early_stopping_rounds = 3,
                    verbose = 0)
    
    best_iter <- cvfit$best_iteration
    mod <- xgboost(data = x, label = y,
                   objective = obj,
                   nrounds = best_iter,
                   verbose = 0)
    
    predict(mod, newdata = newx)
}


.augment_y <- function(x, trt, y) {
    pred1 <- .train_predict_boot(x[trt == 1,], y[trt == 1], newx = x)
    pred0 <- .train_predict_boot(x[trt != 1,], y[trt != 1], newx = x)
    (pred1 + pred0) / 2
}





# Entropy balancing -------------------------------------------------------

#' Obtain entropy balancing weights.
ebal_generalize <- function(x, trt, x_target = NULL, 
                            standardize_x = FALSE, degree = 1,
                            max.iterations = 600, constraint.tolerance = 1e-5)
{
    stopifnot(length(trt) == NROW(x))
    if (!is.null(x_target)) stopifnot(NCOL(x) == NCOL(x_target))
    stopifnot(degree %in% c(1, 2))
    stopifnot(any(trt == 1) && any(trt != 1))
    trt[trt != 1] <- 0
    
    if (is.null(x_target))
    {
        x_target <- x
    }
    
    n <- NROW(x)
    n_t <- NROW(x_target)
    
    if (standardize_x || degree == 2)
    {
        x_all <- rbind(x, x_target)
        if (standardize_x) x_all <- scale(x_all)
        if (degree == 2) 
        {
            x_all <- stats::poly(x_all, degree = 2)
            x_all <- unname(x_all)
            class(x_all) <- "matrix"
        }
        x <- x_all[1:n,]
        x_target <- x_all[n+1:n_t,]
    }
    
    moments <- colMeans(x_target)
    weights <- numeric(n)
    weights[trt == 1] <- .eb(x[trt == 1,], moments, 
                             max.iterations = max.iterations, 
                             constraint.tolerance = constraint.tolerance)
    weights[trt != 1] <- .eb(x[trt != 1,], moments,
                             max.iterations = max.iterations, 
                             constraint.tolerance = constraint.tolerance)
    
    weights
}

#' Helper function of entropy balancing.
.eb <- function(x, moments, sumw = 1,
                max.iterations, constraint.tolerance)
{
    stopifnot(NCOL(x) == length(moments))
    
    x <- cbind(x, 1)
    moments <- c(moments, 1)
    
    lambda <- rep(0, length(moments))
    i <- 0
    repeat {
        w <- exp(drop(x %*% lambda) - 1)
        grad <- moments - drop(w %*% x)
        if (sum(abs(grad)) < constraint.tolerance) break
        if (i == max.iterations) {
            warning("Max iterations without convergence!")
            break
        }
        hess <- -t(x) %*% (w * x)
        old <- lambda
        lambda <- old - solve(hess - 1e-5 * diag(NROW(hess)), grad)
        i <- i + 1
    }
    
    w / sum(w) * sumw
}




# Linear ITR learning -----------------------------------------------------

#' Outcome weighted learning for linear ITR.
fit_linear_itr_owl <- function(x, trt, y, weights, 
                               augment_y = FALSE, 
                               normalize_weight = TRUE,
                               normalize_beta = FALSE) 
{
    if (normalize_weight)
    {
        weights[trt == 1] <- weights[trt == 1] / sum(weights[trt == 1])
        weights[trt != 1] <- weights[trt != 1] / sum(weights[trt != 1])
    }
    
    if (augment_y) y <- y - .augment_y(x, trt, y)
    
    contr <- (2 * (trt == 1) - 1) * weights * y
    .fit_linear_itr_inner(x, contr, normalize_beta = normalize_beta)
}

#' Obtain linear ITR based on imputed outcomes on the target sample.
fit_linear_itr_contrast <- function(x, trt, y, x_target,
                                    normalize_beta = FALSE)
{
    if (missing(x_target)) x_target <- x
    
    pred1 <- .train_predict_boot(x[trt == 1,], y[trt == 1], x_target, B = 1)
    pred0 <- .train_predict_boot(x[trt != 1,], y[trt != 1], x_target, B = 1)
    
    .fit_linear_itr_inner(x_target, pred1 - pred0, 
                          normalize_beta = normalize_beta)
}

#' Helper function for ITR learning.
.fit_linear_itr_inner <- function(x, contr, normalize_beta)
{
    a <- contr >= 0
    w <- abs(contr)
    
    fit <- suppressWarnings(glm(a ~ x, weights = w, family = "binomial"))
    
    beta <- unname(coef(fit))
    
    if (normalize_beta) beta <- beta / sqrt(sum(beta[-1]^2))
    
    pred_fn <- function(x) drop(x %*% beta[-1]) + beta[1]
    
    list(beta = beta, 
         pred_fn = pred_fn)
}
