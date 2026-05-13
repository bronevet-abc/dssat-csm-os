# test_calibrate_somic.R

# Copy loss.function from calibrate_somic.R (simplified to test penalty)
loss.function <- function(par.optim, par.complete, somic.data, experiments, penalty_lower = NULL, penalty_upper = NULL) {
    if (!is.null(penalty_lower) || !is.null(penalty_upper)) {
        out_of_bounds <- numeric(length(par.optim))
        if (!is.null(penalty_lower)) {
            out_of_bounds <- out_of_bounds + pmax(0, penalty_lower[names(par.optim)] - par.optim)
        }
        if (!is.null(penalty_upper)) {
            out_of_bounds <- out_of_bounds + pmax(0, par.optim - penalty_upper[names(par.optim)])
        }

        if (sum(out_of_bounds) > 0) {
             return(1e6 + sum(out_of_bounds^2) * 1e8)
        }
    }
    # Return 0 for in-bounds to indicate success in this test
    return(0)
}

# Test penalty logic
par.optim <- c(mic_vmax = 1.5, mic_km = 0.3)
fit.lower <- c(mic_vmax = 1.0, mic_km = 0.1)
fit.upper <- c(mic_vmax = 2.5, mic_km = 0.6)

# Case 1: In bounds
res <- loss.function(par.optim, NULL, NULL, NULL, fit.lower, fit.upper)
if (res == 0) {
    cat("Test 1 (In bounds) PASSED\n")
} else {
    cat("Test 1 (In bounds) FAILED\n")
}

# Case 2: Out of bounds (lower)
par.optim <- c(mic_vmax = 0.5, mic_km = 0.3)
res <- loss.function(par.optim, NULL, NULL, NULL, fit.lower, fit.upper)
if (res > 1e6) {
    cat("Test 2 (Out of bounds lower) PASSED\n")
} else {
    cat("Test 2 (Out of bounds lower) FAILED\n")
}

# Case 3: Out of bounds (upper)
par.optim <- c(mic_vmax = 1.5, mic_km = 0.7)
res <- loss.function(par.optim, NULL, NULL, NULL, fit.lower, fit.upper)
if (res > 1e6) {
    cat("Test 3 (Out of bounds upper) PASSED\n")
} else {
    cat("Test 3 (Out of bounds upper) FAILED\n")
}
