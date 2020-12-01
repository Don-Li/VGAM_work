library( VGAM )

source( "Thomas.R" )

example(acat)  # Get 'fit'; any VGAM model should do really

M <- npred(fit)
nn <- nobs(fit)

# getMethod( "weights", "vglm" )
# VGAM:::weightsvlm( fit, type = "working" )
wz_weights = wweights( fit, matrix.arg = TRUE, deriv.arg = FALSE )

wz_weights_array = m2a_R( wz_weights, M )

wz_weights_chol = apply( wz_weights_array, 3, function(matrix_){
    col_index = v_seq2( 1:M, 1 )
    row_index = rep( 1:M, 1:M )
    chol( matrix_ )[ cbind(col_index, row_index) ]
    } )

X.vlm <- model.matrix(fit, type = "vlm")

xmat = X.vlm
cc = wz_weights_chol
upper = TRUE

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
              as.integer(M), as.integer(R), 
              n = as.integer(n),
                # 2,
              wkcc = double(M * M), wk2 = double(M * R),
              row_index = as.integer(index$row), col_index = as.integer(index$col),
              dim = as.integer(dimm.value),
              as.integer(upper), NAOK = TRUE)

mux111ccc_R( wz_weights_chol, as.double(t(xmat)), M, R, n,
    index$row, index$col, dimm.value, upper = TRUE )
fred$b

wz_weights_chol2 = apply( wz_weights_array, 3, function(matrix_){
    chol( matrix_ )
    } )

matrix( wz_weights_chol2[,1], 2, 2 )
wz_weights_chol_array[,,1]

wz_weights_chol_array = m2a_R( t( wz_weights_chol ), M, TRUE )
xmat_array = design_matrix_to_array_R( xmat, y_dim = M, n_obs = n )



design_matrix = xmat
y_dim = M
n_obs = n

# n_ = 1
result_ = do.call( cbind, lapply( 1:n, function(n_){
    t( wz_weights_chol_array[,,n_] %*% xmat_array[,,n_] )
} ) )

abs( c(result_) - fred$b )



library( Rcpp )
sourceCpp("cpp_stuff.cpp")

# xmat_to_array( xmat, M ) == xmat_array



index <- iam(NA, NA, M, both = TRUE, diag = TRUE)
# 
# m2a_cpp( t( wz_weights_chol ), as.integer(M), 1L,
#     as.integer(index$row.index), as.integer(index$col.index) )

microbenchmark::microbenchmark(
a = {
# wz_weights_array = m2a_R( wz_weights, M, upper_triangular = T )
# wz_weights_array2 =  m2a_cpp( wz_weights, as.integer(M), 0L,
#     as.integer(index$row.index-1), as.integer(index$col.index-1) )
# 
# xx = xmat_to_array( xmat, M )
# e = lapply( 1:8, function(i){
#     t( chol(wz_weights_array[,,i]) %*% xx[,,i] )
# } )
# c( do.call( cbind, e ) )
    mux111ccc_cpp( wz_weights, xmat, as.integer(M), 0L,
        as.integer(index$row.index-1), as.integer(index$col.index-1 ) )
},
    b = {
           fred <- .C("mux111ccc",
              as.double(cc), b = as.double(t(xmat)),
              as.integer(M), as.integer(R), 
              n = as.integer(n),
                # 2,
              wkcc = double(M * M), wk2 = double(M * R),
              row_index = as.integer(index$row), col_index = as.integer(index$col),
              dim = as.integer(dimm.value),
              as.integer(upper), NAOK = TRUE)
    }, times = 100000 )
    


fred$b


# n = 1
# tlocal = 1
# for ( tlocal in 1:n ){
#     # The first ilocal loop
#     # wkcc[] = m2a_R( t(cc[ ,tlocal, drop = F ]), 2 )[,,1]
#     wkcc[] = c( cc[,tlocal][c(1)], 0, cc[,tlocal][c(3,2)] )
#     
#     # pd2 = (txmat)
#     # pd_index = 1
#     # for ( i in 1:M ){
#     #     for ( j in 1:R ){
#     #         wk2[ (i-1) + (j-1) * M + 1 ] = pd2[pd_index]
#     #         pd_index = pd_index + 1
#     #     }
#     # }
#     pd2 = (t(xmat))
#     wk2 = t(pd2)[ {
#         cbind(
#             rep( 1:M, times = R ) + (tlocal-1) * M,
#             rep( 1:R, each = M ) 
#         )
#     } ]
#     
#     y = 0
#     # txmat
#     ilocal = 1
#     for ( ilocal in 1:M ){
#         lowlim = ilocal
#         jlocal = 1
#         for ( jlocal in 1:R ){
#             slocal = 0
#             klocal = lowlim
#             for ( klocal in lowlim:M ){
#                 slocal[klocal] = 
#                     wk2[ (klocal-1) + (jlocal-1) * M + 1 ] *
#                         wkcc[(ilocal-1) + (klocal-1) * M + 1]
#             }
#             y = c(y, sum(slocal))
#             txmat[(jlocal-1) + (ilocal-1) * R + (tlocal-1) * n + 1 ] = sum(slocal)
#         }
#     }
# }
# 
# ( cbind( fred$b,
#     as.numeric((txmat)) ) )
# ( cbind( fred1$b,
#     as.numeric((txmat)) ) )







