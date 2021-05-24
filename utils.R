#' Fit a linear or logistic regression and output predictions
#'
#' @param x covariates.
#' @param label labels, either numeric or logical.
#' @param x_extra additional covariates that we want to get predictions.
#' @return A list of predictions: pred(, pred_extra).
.fit <- function(x, label, x_extra = NULL)
{
   faml <- switch(typeof(label),
                  logical = "binomial",
                  double = "gaussian")
   
   mod <- suppressWarnings(glm(label ~ x, family = faml))
   
   .get_pred <- function(x) {
      coef <- mod$coefficients
      coef[is.na(coef)] <- 0
      z <- coef[1] + drop(x %*% coef[-1])
      if (is.logical(label)) z <- plogis(z)
      z
   }
   
   out <- list(pred = .get_pred(x))
   if (!is.null(x_extra)) out$pred_extra <- .get_pred(x_extra)
   out
}


#' Estimate the propensity score
#' 
#' @param xs covariates of the source sample.
#' @param as treatment assignment (binary) of the source sample.
#' @param xt covariates of the target sample.
#' @param x_test covariates of additional test data.
#' @return A list of estimated propensity scores: s(, t, test).
get_prop_score <- function(xs, as, xt = NULL, x_test = NULL)
{
   tmp <- .fit(xs, as == 1, rbind(xt, x_test))
   out <- list(s = tmp$pred)
   if (!is.null(xt)) out$t <- tmp$pred_extra[1:NROW(xt)]
   if (!is.null(x_test)) out$test <- tmp$pred_extra[NROW(xt) + 1:NROW(x_test)]
   out
}


#' Estimate the participation probability
#' 
#' @param xs,xt,x_test covariates.
#' @return A list of estimated participation probabilities: s, t(, test).
get_prob_source <- function(xs, xt, x_test = NULL)
{
   pred <- .fit(rbind(xs, xt), rep(c(TRUE, FALSE), c(NROW(xs), NROW(xt))), x_test)
   out <- list(s = pred$pred[1:NROW(xs)],
               t = pred$pred[NROW(xs) + 1:NROW(xt)])
   if (!is.null(x_test)) out$test <- pred$pred_extra
   out
}
