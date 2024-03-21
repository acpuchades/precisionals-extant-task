tidy.coxme <- function(x, exponentiate = FALSE, conf.int = 0.95, ...) {
    beta <- x$coefficients
    nvar <- length(beta)
    nfrail <- nrow(x$var) - nvar
    nn <- c("estimate", "exp()", "std.error", "statistic", "p.value")
    se <- sqrt(diag(as.matrix(x$var))[nfrail + 1:nvar])
    z <- qnorm((1 + conf.int) / 2, 0, 1)
    ret <- data.frame(
        "term"      = names(beta),
        "estimate"  = beta,
        "std.error" = se,
        "statistic" = beta / se,
        "p.value"   = 1 - pchisq((beta / se)^2, 1),
        "conf.low"  =  beta - z * se,
        "conf.high" =  beta + z * se
    )
    if (exponentiate) {
        ret$estimate <- exp(ret$estimate)
        ret$conf.low <- exp(ret$conf.low)
        ret$conf.high <- exp(ret$conf.high)
    }
    rownames(ret) <- c(1:nrow(ret))
    ret
}

var.mean <- function(x, ...) {
    var(x, ...) / length(x)
}
