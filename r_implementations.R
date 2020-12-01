v_seq = function( from, v_to ){
    c( unlist( sapply( v_to, seq, from = from ) ) )
}
v_seq2 = function( v_from, to ){
    c( unlist( sapply( v_from, seq, to = to ) ) )
}

m2a_R = function( matrix_, dim_matrix, upper_triangular = FALSE ){
    N = nrow(matrix_)
    return_array = array( 0, c(dim_matrix, dim_matrix, N) )
    row_index = v_seq( 1, dim_matrix:1 )
    col_index = rep( 1:dim_matrix, 1:dim_matrix )

    replacement_matrix = matrix( 0, dim_matrix, dim_matrix)
    
    i = 1
    for ( i in 1:N ){
        diag(replacement_matrix) = matrix_[ i, 1:dim_matrix ]
        j = 1
        for ( j in 1:length(row_index) ){
            if ( !upper_triangular ){
                if ( row_index[j] != col_index[j] ){
                    replacement_matrix[ row_index[j], col_index[j] ] = matrix_[i,j]
                    replacement_matrix[ col_index[j], row_index[j] ] = matrix_[i,j]
                }
            } else{
                replacement_matrix[ row_index[j], col_index[j] ] = matrix_[i,j]
            }
        }
        return_array[,,i] = replacement_matrix
    }
    # if ( !upper_triangular ){
    #     for ( i in 1:N ){
    #         # return_array[,,i] = matrix_[ 1, ]
    #     }
    # } else{
    #     id_matrix = cbind(row_index, col_index)
    #     diagonal_id = id_matrix[,1] != id_matrix[,2]
    #     id_matrix = rbind( id_matrix, id_matrix[ diagonal_id, 2:1 ] )
    #     for ( i in 1:N ){
    #         return_array[,,i][id_matrix] = matrix_[ i, ]
    #     }
    # 
    # }
    return_array
}

design_matrix_to_array_R = function( design_matrix, y_dim, n_obs ){
    index = rep( 1:n_obs, each = y_dim )
    output_array = array( 0, dim = c(y_dim, ncol(design_matrix), n_obs) )
    for ( ob in 1:n_obs ){
        output_array[,,ob] = design_matrix[ index == ob, ]
    }
    output_array
}

a2m_R = function( array_ ){
    if ( length(dim(array_)) == 2 ){
        stop()
        return_matrix_ncol =
            return_matrix_nrow = nrow(array_)-1
        N = 2
    } else{
        return_matrix_nrow = dim(array_)[3]
        return_matrix_ncol = dim(array_)[2] + 1
        if ( is.na(return_matrix_nrow) ){
            return_matrix_nrow = return_matrix_ncol
        }
        N = return_matrix_nrow
    }
    
    return_matrix = matrix( 0, nrow = return_matrix_nrow, 
        ncol = return_matrix_ncol )
    replacement_vector = numeric( return_matrix_ncol )
    
    
    
    if ( length(dim(array_)) == 3 ){ 
        row_index = v_seq( 1, return_matrix_nrow:1 )
        col_index = rep( 1:return_matrix_nrow, 1:return_matrix_nrow )
        
        i = 1
        for ( i in 1:N ){
            sub_matrix = array_[,,i]
            replacement_vector = sub_matrix[ cbind( row_index, col_index) ]
            return_matrix[ i, ] = replacement_vector
        }
    } else{
        row_index = v_seq( 1, return_matrix_nrow:1 )
        col_index = rep( 1:return_matrix_nrow, 1:return_matrix_nrow )
        
        return_matrix[ cbind( row_index, col_index ) ] = array_
        
    }
    return_matrix
}

mux111ccc_R = function( cc, txmat, M, R, n,
    row_index, col_index,
    dimm, upper = TRUE
    ){
    wkcc = double(M * M)
    wk2 = double(M * R)
    
    row_index = row_index - 1
    col_index = col_index - 1
    
    # for(ilocal in 1:dimm) {
    #     row_index[ilocal] = row_index[ilocal] - 1
    #     col_index[ilocal] = col_index[ilocal] - 1
    # }
    
    c_index = 1
    pd_index = 1
    tlocal = 1
    
    for ( tlocal in 1:n ){
        for ( ilocal in 1:dimm ){
            if ( upper == FALSE ){
                wkcc[row_index[ilocal] + col_index[ilocal] * M] =
                    wkcc[col_index[ilocal] + row_index[ilocal] * M] = cc[c_index]
            } else{
                wkcc[row_index[ilocal] + col_index[ilocal] * M + 1] = cc[c_index]
            }
            c_index = c_index + 1
        }
        
        pd2 = (t(xmat))
        for ( i in (1:M) ){
            for ( j in (1:R) ){
                wk2[ (i-1) + (j-1) * M + 1 ] = pd2[pd_index]
                pd_index = pd_index + 1
            }
        }
        
        ilocal = 1
        for ( ilocal in (1:M) ){
            lowlim = ifelse( upper == FALSE, 0, ilocal )
            jlocal = 1
            for ( jlocal in (1:R) ){
                slocal = 0
                for ( klocal in (lowlim:M) ){
                    slocal = slocal + 
                        wk2[(klocal-1) + (jlocal-1) * M + 1] *
                        wkcc[(ilocal-1) + (klocal-1) * M + 1]
                }
                txmat[ (jlocal-1) + (ilocal-1) * R + (tlocal-1) * n + 1] = sum(slocal, na.rm = T)
            }
        }
        txmat = txmat
        # txmat2 = txmat2
    }
    txmat
}










