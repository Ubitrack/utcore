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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_GELS_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_GELS_HPP

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
     * gels() computes the solution to a system of linear equations 
     * A * X = B, where A is an M-by-N matrix and X and B are N-by-NRHS 
     * matrices.
     *
     * The QR or LQ decomposition is used to compute the result.
	 *
	 * trans is either 'N' which solves A * X = B 
	 * or 'T' which solves A^T * X = B
	 *
	 * Note: The result X is returned in B, A is overwritten by the
	 * decomposition.
     */ 

    namespace detail {

      inline 
      void gels (char const trans, int const m, int const n, int const nrhs,
                 float* a, int const lda, float* b, int const ldb, 
				 float* work, int const lwork, int* info) 
      {
        LAPACK_SGELS (&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, info);
      }

      inline 
      void gels (char const trans, int const m, int const n, int const nrhs,
                 double* a, int const lda, double* b, int const ldb, 
				 double* work, int const lwork, int* info) 
      {
        LAPACK_DGELS (&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, info);
      }

      inline 
      void gels (char const trans, int const m, int const n, int const nrhs,
                 traits::complex_f* a, int const lda, traits::complex_f* b, int const ldb, 
				 traits::complex_f* work, int const lwork, int* info) 
      {
        LAPACK_CGELS (&trans, &m, &n, &nrhs, traits::complex_ptr( a ), &lda, traits::complex_ptr( b ), &ldb, 
			traits::complex_ptr( work ), &lwork, info);
      }

      inline 
      void gels (char const trans, int const m, int const n, int const nrhs,
                 traits::complex_d* a, int const lda, traits::complex_d* b, int const ldb, 
				 traits::complex_d* work, int const lwork, int* info) 
      {
        LAPACK_ZGELS (&trans, &m, &n, &nrhs, traits::complex_ptr( a ), &lda, traits::complex_ptr( b ), &ldb, 
			traits::complex_ptr( work ), &lwork, info);
      }
      
    } 

    template <typename MatrA, typename MatrB, typename WVec>
    inline
    int gels (char const trans, MatrA& a, MatrB& b, WVec& work) {

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

      assert( trans == 'N' || trans == 'T' );
      int const m = traits::matrix_size1 (a);
      int const n = traits::matrix_size2 (a);
      assert (traits::matrix_size1 (b) == ( trans == 'T' ? n : m ) ); 
	  int const mn = std::min( m, n );
      assert (traits::vector_size (work) >= mn + std::max( mn, traits::matrix_size2( b ) ) ) ; 

      int info; 
      detail::gels (trans, m, n, traits::matrix_size2 (b), 
                    traits::matrix_storage (a), 
                    traits::leading_dimension (a),
                    traits::matrix_storage (b),
                    traits::leading_dimension (b),
                    traits::vector_storage (work),  
					traits::vector_size (work),
                    &info);
      return info; 
    }

    template <typename MatrA, typename MatrB>
    inline
    int gels ( char trans, MatrA& a, MatrB& b ) {
      // with 'internal' work vector 
      
      // gels() errors: 
      //   if (info == 0), successful

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS
      typedef typename traits::matrix_traits<MatrA>::value_type val_t; 
#else 
      typedef typename MatrA::value_type val_t; 
#endif 
      int const m = traits::matrix_size1 (a);
      int const n = traits::matrix_size2 (a);
      int const nrhs = traits::matrix_size2 (b);
	  int const mn = std::min( m, n );
	  traits::detail::array <val_t> wvec( mn + std::max( mn, nrhs ) ); 
	  
	  return gels( trans, a, b, wvec );
    }
	
}}}}

#endif 
