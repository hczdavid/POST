# Code from CKAT package
#Debashis Ghosh and Xiang Zhan. "CKAT software" (2016)
#Available at: http://works.bepress.com/debashis_ghosh/75/
#
#' @keywords internal
CKAT.bin <-
  function (y, K,  X=NULL) {
    n <- length(y)
    if (is.null(X)) {
      X1 <-  matrix(rep(1, length(y)), ncol=1)
    } else {
      X1 <- model.matrix(~. , as.data.frame(X))
    }

    glmfit <- glm(y ~ X1-1, family = binomial)

    betas <- glmfit$coef
    mu  <- glmfit$fitted.values
    eta <- glmfit$linear.predictors
    res.wk <- glmfit$residuals
    res 	<- y - mu

    w   <- mu*(1-mu)
    sqrtw <- sqrt(w)

    adj <- sum((sqrtw * res.wk)^2)

    DX12 <- sqrtw * X1


    qrX <- qr(DX12, tol = 1e-7)
    Q <- qr.Q(qrX)
    Q <- Q[, 1:qrX$rank, drop=FALSE]

    P0 <- diag(length(y)) - Q %*% t(Q)

    DKD <- tcrossprod(sqrtw) * K
    tQK <- t(Q) %*% DKD
    QtQK <- Q %*% tQK
    PKP1 <- DKD - QtQK - t(QtQK) + Q %*% (tQK %*% Q) %*% t(Q)
    q1 <- as.numeric(res %*% K %*% res)
    q1 = q1 / adj
    ee1 = eigen(PKP1 - q1 * P0, symmetric = T, only.values=T)
    lambda1 = ee1$values[abs(ee1$values) >= 1e-10]
    p1 <- davies(0, lambda=lambda1, acc=1e-6)$Qq

    return(p1)
  }



#' @keywords internal
CKAT.cont <-
  function (y, K, X=NULL) {
    n <- length(y)
    if (is.null(X)) {
      X1 <-  matrix(rep(1, length(y)), ncol=1)
    } else {
      X1 <- model.matrix(~. , as.data.frame(X))
    }

    px <- ncol(X1)
    mod <- lm(y ~ X1-1)
    res = resid(mod)
    s2 <- sum(res^2)

    D0  <- diag(length(y))  #diag(mu*(1-mu)) for logistic
    P0  <- D0 - X1 %*% solve(t(X1)%*%X1) %*% t(X1)
    PKP <- P0 %*% K %*% P0

    q <- as.numeric(res %*% K %*% res / s2)
    ee = eigen(PKP - q * P0, symmetric = T)
    lambda0 = ee$values[abs(ee$values) >= 1e-10]
    p1 <- davies(0, lambda=lambda0)$Qq
    return(p1)
  }
