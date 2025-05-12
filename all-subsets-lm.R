# regsubsets does not work nicely with factor variables with more than 2 levels,
# because it treats the constrasts as separate variables. Here, we handle this
# manually adding and removing factors, and using regsubsets for the quantitative
# variables.
#
# lower_scope is a formula containing all terms that must be included in the model
# (defaults to ~ 1, just the intercept, which is always included anyway)
#
# Examples:
# # Consider all subsets of (x1, x2, x3)
# all_subsets_lm(lm(y ~ x1 + x2 + x3, data = data))
# # Consider all subsets of (x1, x2, x3), but x1 must be included
# all_subsets_lm(lm(y ~ x1 + x2 + x3, data = data), lower_scope = ~ x1)
# # Use BIC instead of AIC
# all_subsets_lm(lm(y ~ x1 + x2 + x3, data = data), criterion = 'BIC')
# # Limit the number of variables considered
# all_subsets_lm(lm(y ~ x1 + x2 + x3, data = data), p_max = 2)
all_subsets_lm <- function(
  model,
  lower_scope,
  criterion = c('AIC', 'AICc', 'BIC', 'adjR2', 'Cp'),
  p_min = 0,
  p_max = Inf
) {
  criterion <- match.arg(criterion)

  drop_scope <- if (missing(lower_scope)) {
    drop.scope(model)
  } else {
    drop.scope(model, lower_scope)
  }

  y <- model$residuals + model$fitted.values
  tss <- sum((y - mean(y)) ^ 2)
  n <- length(y)
  sigma_squared_full <- summary(model)$sigma ^ 2
  ic_fn <- function(rss, p) {
    k <- p + 1
    if (criterion == 'adjR2') {
      # NOTE: we use - because we want to maximise adjR2
      -(1 - (rss / (n - p - 1)) / (tss / (n - 1)))
    } else if (criterion == 'AIC') {
      n * log(rss / n) + 2 * k
    } else if (criterion == 'AICc') {
      n * log(rss / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)
    } else if (criterion == 'BIC') {
      n * log(rss / n) + k * log(n)
    } else if (criterion == 'Cp') {
      rss / sigma_squared_full + 2 * k - n
    }
  }

  term_names <- attr(model$terms, "term.labels")
  term_i <- seq_along(term_names)
  p_full <- length(term_names)

  if (p_min > p_full) {
    stop('p_min is greater than the number of terms in the model')
  }

  # Identify any factor variables with more than 2 levels
  X_all <- model.matrix(model)
  X_assign <- attr(X_all, 'assign')
  is_large_factor <- sapply(term_i, function(i) {
    sum(X_assign == i) > 1
  })
  large_factor_terms <- term_names[which(is_large_factor)]
  droppable_large_factor_terms <- intersect(large_factor_terms, drop_scope)
  droppable_large_factor_terms_i <- match(droppable_large_factor_terms, term_names)
  n_droppable_large_factors <- length(droppable_large_factor_terms_i)

  retain_scope <- setdiff(term_names, drop_scope)
  retain_terms_index <- match(retain_scope, term_names)

  terms_from_col_names <- function(col_names) {
    column_indices <- match(col_names, colnames(X_all))
    term_indices <- unique(sort(X_assign[column_indices]))
    term_names[term_indices]
  }

  # We will try dropping 0, 1, ..., n_large_factors factor variables,
  # and use regsubsets for the rest
  best_model <- list(ic = Inf, variables = NULL)
  for (n_drop_factors in 0 : n_droppable_large_factors) {
    cases <- combn(length(droppable_large_factor_terms_i), n_drop_factors)
    output <- lapply(seq_len(ncol(cases)), function(i) {
      case_i <- cases[, i]
      drop_terms_i <- droppable_large_factor_terms_i[case_i]
      retained_large_factor_terms_i <- setdiff(droppable_large_factor_terms_i, drop_terms_i)
      drop_cols_i <- which(X_assign %in% drop_terms_i)
      if (length(drop_cols_i) > 0) {
        X_dropped <- X_all[, -drop_cols_i, drop = FALSE]
      } else {
        X_dropped <- X_all
      }
      retained_colnames_i <- colnames(X_all)[
        X_assign %in% c(retain_terms_index, retained_large_factor_terms_i)
      ]
      new_retained_i <- match(retained_colnames_i, colnames(X_dropped))

      if (ncol(X_dropped) < p_min + 1L) {
        return(list(ic = Inf, variables = NULL))
      }

      p_max_remaining <- p_max - length(retained_colnames_i)
      if (p_max_remaining < 0) {
        return(list(ic = Inf, variables = NULL))
      } else if (p_max_remaining == 0) {
        X_retained <- X_dropped[, c(1, new_retained_i), drop = FALSE]
        rss <- sum(lm.fit(X_retained, y)$residuals ^ 2)
        ic <- ic_fn(rss, ncol(X_retained) - 1L)
        variable_names <- terms_from_col_names(colnames(X_retained))
        return(list(ic = ic, variables = variable_names))
      }

      regsubsets_output <- summary(leaps::regsubsets(
        X_dropped,
        y,
        force.in = c(1, new_retained_i),
        intercept = FALSE,  # Already included in X_dropped
        nvmax = p_max + 1L  # +1 for the intercept
      ))
      p_cases <- rowSums(regsubsets_output$which) - 1L
      rss_cases <- regsubsets_output$rss
      ic_cases <- ic_fn(rss_cases, p_cases)
      for (i in seq_along(p_cases)) {
        if (p_cases[i] < p_min) {
          ic_cases[i] <- Inf
        }
      }

      best_case_index <- which.min(ic_cases)
      best_case_ic <- ic_cases[best_case_index]

      best_column_names <- colnames(regsubsets_output$which)[
        regsubsets_output$which[best_case_index, ]
      ]
      best_variable_names <- terms_from_col_names(best_column_names)

      list(variables = best_variable_names, ic = best_case_ic)
    })
    best_case <- output[[which.min(sapply(output, function(m) m$ic))]]
    if (best_case$ic < best_model$ic) {
      best_model <- best_case
    }
  }
  best_model$model <- update(
    model,
    paste('~ ', paste(best_model$variables, collapse = ' + '))
  )
  if (criterion == 'adjR2') {
    best_model$ic <- -best_model$ic
  }
  best_model
}
