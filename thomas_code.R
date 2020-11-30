library( VGAM )

example(acat)  # Get 'fit'; any VGAM model should do really


M <- npred(fit)  # M > 1 is best to see what's going on
nn <- nobs(fit)

wzedd <- weights(fit, type = "working")
# Nb. wzedd absorbs the prior weights in it.
# Few rows, many cols:
UU <- vchol(wzedd, M = M, n = nn, silent = TRUE)
# Note: UU is is matrix-band format (Sec.18.3.5 of my book)
X.vlm <- model.matrix(fit, type = "vlm")
UU.X.vlm <- mux111(cc = UU, xmat = X.vlm, M = M)
UU.X.vlm


wzar <- m2a(wzedd, M = M)
wzar
chol(wzar[,,1])  # Same as UU[, 1]


-------------------------------------------------------------




mux111 <- function(cc, xmat, M, upper = TRUE) {
# 19990813
#
# 20170605: some modifications to get hdeff() going.

# Input:
# 1. cc is dimm x n, where dimm is
#    usually M or M(M+1)/2, e.g., is output from vchol().
#    That is, cc is some matrix in matrix-band format.
# 2. xmat is n*M x R. Each sub-block is M x R and is
#    multiplied by each column of cc.
#
# Output:
# 1. (n*M) x R matrix, such that, in old notation,
#    ans[1:M, ] <- cc[, , 1] %*% xmat[1:M, ], etc.
#
# Notes:
# 1. It is specifically geared to vglm.fit(), and is a special case
#    of mux11().
# 2. If M = 1 then cc cannot be simply be a vector.
#    20170605; not sure about the above line.
# 3. NAs in xmat and/or cc are okay.
# 4. cc represents an upper-triangular matrix by default, but if not
#    then it is symmetric (set upper = FALSE).

   if (!is.matrix(xmat))
     xmat <- as.matrix(xmat)
   if (!is.matrix(cc))
     cc <- t(as.matrix(cc))
   R <- ncol(xmat)
   n <- nrow(xmat) / M
   index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
   dimm.value <- nrow(cc)  # Usually M or M(M+1)/2

   fred <- .C("mux111ccc",
              as.double(cc), b = as.double(t(xmat)),
              as.integer(M), as.integer(R), as.integer(n),
              wkcc = double(M * M), wk2 = double(M * R),
              as.integer(index$row), as.integer(index$col),
              as.integer(dimm.value),
              as.integer(upper), NAOK = TRUE)

   ans <- fred$b
   dim(ans) <- c(R, n * M)
   d <- dimnames(xmat)
   dimnames(ans) <- list(d[[2]], d[[1]])
   t(ans)
}