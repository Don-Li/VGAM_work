#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube m2a_cpp( arma::dmat matrix_, int dim_matrix, int upper_triangular,
    arma::ucolvec row_index, arma::ucolvec col_index ){

    int N = matrix_.n_rows;
    arma::cube return_array(dim_matrix, dim_matrix, N, arma::fill::zeros);

    //  I think temporary assignment might be slow?
    arma::dmat temp_matrix( dim_matrix, dim_matrix, arma::fill::zeros );

    for ( int ilocal = 0; ilocal < N; ilocal++ ){
        temp_matrix.diag() = matrix_.row(ilocal).cols(0, dim_matrix-1);

        for ( int jlocal = 0; jlocal < row_index.n_rows; jlocal++ ){
            // Need to improve this somehow
            // Possible to use std::vector<int> row_index instead of arma::ucolvec row_index?
            // or maybe combine it with the matrix indexing from the other versions of this function 
            int row_id = arma::conv_to<int>::from( row_index.row(jlocal) );
            int col_id = arma::conv_to<int>::from( col_index.row(jlocal) );

            if ( row_id != col_id ){
                if ( upper_triangular == 0 ){
                    temp_matrix(row_id, col_id) =
                        temp_matrix(col_id, row_id) =
                            matrix_(ilocal, jlocal);
                } else{
                    temp_matrix(row_id, col_id) =
                        matrix_(ilocal, jlocal);
                }
            }
        }
        return_array.slice(ilocal) = temp_matrix;
    }
    return return_array;
}

// [[Rcpp::export]]
arma::dmat m2a_cpp_1( arma::dmat band_matrix, arma::umat index_matrix, 
    int upper_triangular, int n_triangular_elements, int band_row ){
    
    arma::dmat output_matrix( n_triangular_elements, n_triangular_elements, arma::fill::zeros );
    arma::uvec indices = arma::sub2ind( arma::size(output_matrix), index_matrix);
    output_matrix.elem(indices) = band_matrix.row(band_row);
    
    if ( upper_triangular == 0 ){
        output_matrix = arma::symmatu(output_matrix);
    }
    return output_matrix;
}

// [[Rcpp::export]]
arma::dcube m2a_cpp_x( arma::dmat band_matrix, arma::umat index_matrix, 
    int upper_triangular, int n_triangular_elements, int N ){
    
    arma::cube return_array(n_triangular_elements, n_triangular_elements, 
        N, arma::fill::zeros);
    arma::uvec indices = arma::sub2ind( arma::size(return_array.slice(0)), index_matrix);

    for ( arma::uword ilocal = 0; ilocal < N; ilocal++ ){
        return_array.slice(ilocal).elem(indices) = band_matrix.row(ilocal);
        
        if ( upper_triangular == 0 ){
            return_array.slice(ilocal) = arma::symmatu(return_array.slice(ilocal));
        }
    }

    return return_array;
}

// [[Rcpp::export]]
arma::dcube m2a_cpp_x2( arma::dmat band_matrix, arma::umat index_matrix, 
    int upper_triangular, int n_triangular_elements, int N ){
    
    arma::cube return_array(n_triangular_elements, n_triangular_elements, 
        N, arma::fill::zeros);
    arma::uvec indices = arma::sub2ind( arma::size(return_array), index_matrix );

    for ( arma::uword ilocal = 0; ilocal < N; ilocal++ ){
        index_matrix.row(2).fill(ilocal);
        indices = arma::sub2ind( arma::size(return_array), index_matrix );
        return_array.elem(indices) = band_matrix.row(ilocal);
        
        if ( upper_triangular == 0 ){
            return_array.slice(ilocal) = arma::symmatu(return_array.slice(ilocal));
        }
    }

    return return_array;
}

// [[Rcpp::export]]
arma::dcube m2a_cpp_x3( arma::dmat band_matrix, arma::umat index_matrix, 
    int upper_triangular, int n_triangular_elements, int N ){
    
    arma::cube return_array(n_triangular_elements, n_triangular_elements, 
        N, arma::fill::zeros);
    arma::uvec indices = arma::sub2ind( arma::size(return_array), index_matrix );

    for ( arma::uword ilocal = 0; ilocal < N; ilocal++ ){
        // index_matrix.row(2).fill(ilocal);
        for ( arma::uword jlocal = 0; jlocal < index_matrix.n_cols; jlocal ++ ){
            return_array(indices(jlocal)) = band_matrix(ilocal, jlocal);
        }

        if ( upper_triangular == 0 ){
            return_array.slice(ilocal) = arma::symmatu(return_array.slice(ilocal));
        }
        indices += n_triangular_elements * n_triangular_elements;
    }

    return return_array;
}

// [[Rcpp::export]]
arma::cube xmat_to_array( arma::dmat xmat, int mat_cols ){
    int return_slices = xmat.n_rows / mat_cols;
    int return_rows = xmat.n_rows / return_slices;
    arma::cube return_array( return_rows, xmat.n_cols, return_slices, arma::fill::zeros );
    for ( int slice_ = 0; slice_ < return_slices; slice_++ ){
        // Rcout << slice_ * return_rows << " " << (slice_ + 1) * return_rows - 1 << "\n" ;
        return_array.slice(slice_) = xmat.rows( slice_ * return_rows, (slice_ + 1) * return_rows - 1 );
    }
    return return_array;
}

// [[Rcpp::export]]
arma::dmat mux111ccc_cpp( arma::dmat matrix_, arma::dmat xmat, int dim_matrix, int upper_triangular,
    arma::ucolvec row_index, arma::ucolvec col_index ){
    
    arma::cube weight_array = m2a_cpp( matrix_, dim_matrix, upper_triangular,
        row_index, col_index);
    arma::cube x_array = xmat_to_array( xmat, dim_matrix );
    
    arma::dmat return_matrix( xmat.n_cols, xmat.n_rows );
    
    for ( int i = 0; i < matrix_.n_rows; i++ ){
        arma::mat K = weight_array.slice(i);
        return_matrix.cols( i * dim_matrix, (i+1) * dim_matrix - 1 ) = 
            arma::trans( arma::chol(K) * x_array.slice(i) );
    }
    
    return return_matrix;
}

