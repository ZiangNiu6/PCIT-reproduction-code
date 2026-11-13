# Original pcmt 0.1.0 R implementation used by the chosen q=.10 run.
adaptive_fdr_bh <- function (T_orig, permute_stat, alpha = 0.05, h = 20L, B_max = 25000L, alternative = c("greater", "less", "two.sided")) 
{
    alternative <- match.arg(alternative)
    h <- as.integer(h)
    B_max <- as.integer(B_max)
    stopifnot(h >= 1L, alpha > 0, alpha < 1, B_max >= 1L)
    stopifnot(is.function(permute_stat))
    m <- length(T_orig)
    if (m == 0L) {
        return(list(rejected = integer(0), futility = integer(0), active = integer(0), p_value = numeric(0), K = integer(0), 
            stop_round = integer(0), t_total = 0L, R_threshold = 0))
    }
    S0 <- score_map(T_orig, alternative)
    K <- integer(m)
    p_val <- rep(1, m)
    status <- rep(0L, m)
    stop_round <- rep(NA_integer_, m)
    active <- seq_len(m)
    t <- 0L
    R_threshold <- 0
    while (length(active) > 0L && t < B_max) {
        t <- t + 1L
        z_perm <- permute_stat(active, t)
        if (length(z_perm) != length(active)) {
            stop("permute_stat must return a vector of length length(active_idx).")
        }
        S_perm <- score_map(z_perm, alternative)
        is_loss <- S_perm >= S0[active]
        K[active] <- K[active] + as.integer(is_loss)
        p_val[active] <- h/(t + h - K[active])
        futile_now <- active[K[active] >= h]
        if (length(futile_now) > 0L) {
            status[futile_now] <- 1L
            stop_round[futile_now] <- t
        }
        active <- active[status[active] == 0L]
        R_threshold <- bh_threshold(p_val, alpha)
        if (R_threshold > 0 && length(active) > 0L) {
            reject_now <- active[p_val[active] <= R_threshold]
            if (length(reject_now) > 0L) {
                status[reject_now] <- 2L
                stop_round[reject_now] <- t
                active <- setdiff(active, reject_now)
            }
        }
    }
    list(rejected = which(status == 2L), futility = which(status == 1L), active = active, p_value = p_val, K = K, stop_round = stop_round, 
        t_total = t, R_threshold = R_threshold)
}

score_map <- function (T_stat, alternative) 
{
    switch(alternative, greater = T_stat, less = -T_stat, two.sided = abs(T_stat))
}

bh_threshold <- function (p, alpha) 
{
    m <- length(p)
    if (m == 0L) 
        return(0)
    p_sorted <- sort(p)
    cuts <- seq_len(m) * alpha/m
    ks <- which(p_sorted <= cuts)
    if (length(ks) == 0L) 
        0
    else max(ks) * alpha/m
}

