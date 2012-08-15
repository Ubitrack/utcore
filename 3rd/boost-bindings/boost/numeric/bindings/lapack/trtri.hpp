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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_TRTRI_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_TRTRI_HPP

#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
#  include <boost/numeric/bindings/traits/detail/symm_herm_traits.hpp>
#  include <boost/static_assert.hpp>
#  include <boost/type_traits/is_same.hpp>
#endif 


namespace boost { namespace numeric { namespace bindings { 

  namespace lapack {

    /////////////////////////////////////////////////////////////////////
    //
    // invert an upper or lower diagonal matrix.
    //
    /////////////////////////////////////////////////////////////////////

    namespace detail {

      inline 
      void trtri (char const uplo, char const diag, int const n, 
                  float* a, int const lda, int* info) 
      {
        LAPACK_STRTRI (&uplo, &diag, &n, a, &lda, info);
      }

      inline 
      void trtri (char const uplo, char const diag, int const n, 
                  double* a, int const lda, int* info) 
      {
        LAPACK_DTRTRI (&uplo, &diag, &n, a, &lda, info);
      }

      inline 
      void trtri (char const uplo, char const diag, int const n, 
                  traits::complex_f* a, int const lda, int* info) 
      {
        LAPACK_CTRTRI (&uplo, &diag, &n, traits::complex_ptr (a), &lda, info);
      }

      inline 
      void trtri (char const uplo, char const diag, int const n, 
                  traits::complex_d* a, int const lda, int* info) 
      {
        LAPACK_ZTRTRI (&uplo, &diag, &n, traits::complex_ptr (a), &lda, info);
      }

      template <typename SymmMatrA> 
      inline
      int trtri (char const uplo, char const diag, SymmMatrA& a) {
        int const n = traits::matrix_size1 (a);
        assert (n == traits::matrix_size2 (a));
        int info; 
        trtri (uplo, diag, n, traits::matrix_storage (a), 
               traits::leading_dimension (a), &info);
        return info; 
      }

    }

    template <typename SymmMatrA> 
    inline
    int trtri (char const uplo, char const diag, SymmMatrA& a) {

      assert (uplo == 'U' || uplo == 'L'); 
      assert (diag == 'N' || uplo == 'U'); 

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<SymmMatrA>::matrix_structure, 
        traits::general_t
      >::value));
#endif

      return detail::trtri (uplo, diag, a); 
    }

    template <typename SymmMatrA> 
    inline
    int trtri (char const uplo, SymmMatrA& a) {
	  return trtri( uplo, 'N', a );
    }

  }

}}}

#endif 
