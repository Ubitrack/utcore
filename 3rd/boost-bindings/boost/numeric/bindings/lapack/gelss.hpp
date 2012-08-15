/*
 * 
 * Copyright (c) Toon Knapen & Kresimir Fresl 2003
 *
 * Permission to copy, modify, use and distribute this software 
 * for any non-commercial or commercial purpose is granted provided 
 * that this license appear on all copies of the software source code.
 *
 * Authors assume no responsibility whatsoever for its use and makes 
 * no guarantees about its quality, correctness or reliability.
 *
 * KF acknowledges the support of the Faculty of Civil Engineering, 
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_GELSS_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_GELSS_HPP

#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>
#include <boost/numeric/bindings/traits/detail/array.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
#  include <boost/static_assert.hpp>
#  include <boost/type_traits/is_same.hpp>
#endif 

#include <cassert>


namespace boost { namespace numeric { namespace bindings { 

  namespace lapack {

    ///////////////////////////////////////////////////////////////////
    //
    // solving a system of linear equations A * X = B in least-squares fashion 
    // 
    ///////////////////////////////////////////////////////////////////

    /* 
     * gelss() computes the solution to a system of linear equations 
     * A * X = B, where A is an M-by-N matrix and X and B are N-by-NRHS 
     * matrices.
     *
     * The singular value decomposition is used to compute the result.
     *
     * Note: The result X is returned in B, A is overwritten by the
     * decomposition.
     *
     * The parameters are:
	 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
	 *          On entry, the M-by-N matrix A.
	 *          On exit, the first min(m,n) rows of A are overwritten with
	 *          its right singular vectors, stored rowwise.
	 *
	 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
	 *          On entry, the M-by-NRHS right hand side matrix B.
	 *          On exit, B is overwritten by the N-by-NRHS solution
	 *          matrix X.  If m >= n and RANK = n, the residual
	 *          sum-of-squares for the solution in the i-th column is given
	 *          by the sum of squares of elements n+1:m in that column.
	 *
	 *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
	 *          The singular values of A in decreasing order.
	 *          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
	 *
	 *  RCOND   (input) DOUBLE PRECISION
	 *          RCOND is used to determine the effective rank of A.
	 *          Singular values S(i) <= RCOND*S(1) are treated as zero.
	 *          If RCOND < 0, machine precision is used instead.
	 *
	 *  RANK    (output) INTEGER
	 *          The effective rank of A, i.e., the number of singular values
	 *          which are greater than RCOND*S(1).
	 *
	 *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
	 *          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )
	 *          For good performance, LWORK should generally be larger.
	 */

    namespace detail {

      inline 
      void gelss(int const m, int const n, int const nrhs,
                 float* a, int const lda, float* b, int const ldb,
                 float* s, const float rcond, int* rank,
                 float* work, int const lwork, int* info) 
      {
        LAPACK_SGELSS (&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, info);
      }

      inline 
      void gelss(int const m, int const n, int const nrhs,
                 double* a, int const lda, double* b, int const ldb,
                 double* s, const double rcond, int* rank,
                 double* work, int const lwork, int* info) 
      {
        LAPACK_DGELSS (&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, info);
      }

      inline 
      void gelss(int const m, int const n, int const nrhs,
                 traits::complex_f* a, int const lda, traits::complex_f* b, int const ldb,
                 float* s, const float rcond, int* rank,
                 traits::complex_f* work, int const lwork, float* rwork, int* info) 
      {
        LAPACK_CGELSS (&m, &n, &nrhs, traits::complex_ptr( a ), &lda, 
            traits::complex_ptr( b ), &ldb, s, &rcond, rank, 
            traits::complex_ptr( work ), &lwork, rwork, info);
      }

      inline 
      void gelss(int const m, int const n, int const nrhs,
                 traits::complex_d* a, int const lda, traits::complex_d* b, int const ldb,
                 double* s, const double rcond, int* rank,
                 traits::complex_d* work, int const lwork, double* rwork, int* info) 
      {
        LAPACK_ZGELSS (&m, &n, &nrhs, traits::complex_ptr( a ), &lda, 
            traits::complex_ptr( b ), &ldb, s, &rcond, rank, 
            traits::complex_ptr( work ), &lwork, rwork, info);
      }
    }

    template <typename MatrA, typename MatrB, typename SVec, typename T, typename WVec>
    inline
    int gelss(MatrA& a, MatrB& b, SVec& s, const T rcond, int& rank, WVec& work) {
      //TODO: same for complex types (with rwork array)
#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure, 
        traits::general_t
      >::value)); 
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrB>::matrix_structure, 
        traits::general_t
      >::value)); 
#endif 

      int const m = traits::matrix_size1 (a);
      int const n = traits::matrix_size2 (a);
      int const nrhs = traits::matrix_size2 (b);
      assert (traits::matrix_size1 (b) == m); 
      int const mn = std::min( m, n );
      assert (traits::vector_size (s) >= mn);
      assert (traits::vector_size (work) >= 3*mn + std::max( std::max( 2*mn, std::max(m, n) ), nrhs ) ) ; 

      int info; 
      detail::gelss (m, n, nrhs, 
                    traits::matrix_storage (a), 
                    traits::leading_dimension (a),
                    traits::matrix_storage (b),
                    traits::leading_dimension (b),
                    traits::vector_storage (s),
                    rcond, &rank,
                    traits::vector_storage (work),  
                    traits::vector_size (work),
                    &info);
      return info; 
    }

    template <typename MatrA, typename MatrB, typename SVec, typename T>
    inline
    int gelss ( MatrA& a, MatrB& b, SVec& s, const T rcond, int& rank ) {
      // with 'internal' work vector 
      
#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS
      typedef typename traits::matrix_traits<MatrA>::value_type val_t; 
#else 
      typedef typename MatrA::value_type val_t; 
#endif 
      int const m = traits::matrix_size1 (a);
      int const n = traits::matrix_size2 (a);
      int const nrhs = traits::matrix_size2 (b);
      int const mn = std::min( m, n );
      traits::detail::array <val_t> wvec( 3*mn + std::max( std::max( 2*mn, std::max(m, n) ), nrhs ) ); 
  
      return gelss( a, b, s, rcond, rank, wvec );
    }
	
}}}}

#endif 
