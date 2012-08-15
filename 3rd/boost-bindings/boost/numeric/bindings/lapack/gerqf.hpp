/*
 * 
 * Copyright (c) Toon Knapen, Karl Meerbergen & Kresimir Fresl 2003
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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_GERQF_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_GERQF_HPP

#include <complex>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/traits/detail/array.hpp>
// #include <boost/numeric/bindings/traits/std_vector.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
#  include <boost/static_assert.hpp>
#  include <boost/type_traits.hpp>
#endif 


namespace boost { namespace numeric { namespace bindings { 

  namespace lapack {

    ///////////////////////////////////////////////////////////////////
    //
    // RQ factorization of a general m x n matrix  A = R * Q
    // 
    ///////////////////////////////////////////////////////////////////

    /* 
     * gerqf() computes the RQ factorization of a rectangular matrix
     * A = R * Q.
     * TODO: Explain more details
     */ 

    namespace detail {
      /* Crossover point for blocked algorithm.
       * Directly call non-blocked version if the problem is smaller 
       * to avoid some expensive checks in lapack.
       */
      static const int gerqfXOver = 128;

      inline 
      void gerqf (int const m, int const n,
                 float* a, int const lda,
                 float* tau, float* work, int const lwork, int& info) 
      {
        LAPACK_SGERQF (&m, &n, a, &lda, tau, work, &lwork, &info);
      }

      inline 
      void gerqf (int const m, int const n,
                 double* a, int const lda,
                 double* tau, double* work, int const lwork, int& info) 
      {
        LAPACK_DGERQF (&m, &n, a, &lda, tau, work, &lwork, &info);
      }

      inline 
      void gerqf (int const m, int const n,
                  std::complex<float>* a, int const lda,
                  std::complex<float>* tau, std::complex<float>* work,
		  int const lwork, int& info) 
      {
        LAPACK_CGERQF (&m, &n,
                      reinterpret_cast<fcomplex_t*> (a), &lda,
                      reinterpret_cast<fcomplex_t*> (tau),
                      reinterpret_cast<fcomplex_t*> (work), &lwork, &info );
      }
      

      inline 
      void gerqf (int const m, int const n,
                  std::complex<double>* a, int const lda,
                  std::complex<double>* tau, std::complex<double>* work,
		  int const lwork, int& info) 
      {
        LAPACK_ZGERQF (&m, &n,
                      reinterpret_cast<dcomplex_t*> (a), &lda,
                      reinterpret_cast<dcomplex_t*> (tau),
                      reinterpret_cast<dcomplex_t*> (work), &lwork, &info );
      }
      
      inline 
      void gerq2 (int const m, int const n,
                 float* a, int const lda,
                 float* tau, float* work, int& info) 
      {
        LAPACK_SGERQ2 (&m, &n, a, &lda, tau, work, &info);
      }

      inline 
      void gerq2 (int const m, int const n,
                 double* a, int const lda,
                 double* tau, double* work, int& info) 
      {
        LAPACK_DGERQ2 (&m, &n, a, &lda, tau, work, &info);
      }

      inline 
      void gerq2 (int const m, int const n,
                  std::complex<float>* a, int const lda,
                  std::complex<float>* tau, std::complex<float>* work,
                  int& info) 
      {
        LAPACK_CGERQ2 (&m, &n,
                      reinterpret_cast<fcomplex_t*> (a), &lda,
                      reinterpret_cast<fcomplex_t*> (tau),
                      reinterpret_cast<fcomplex_t*> (work), &info );
      }
      

      inline 
      void gerq2 (int const m, int const n,
                  std::complex<double>* a, int const lda,
                  std::complex<double>* tau, std::complex<double>* work,
		          int& info) 
      {
        LAPACK_ZGERQ2 (&m, &n,
                      reinterpret_cast<dcomplex_t*> (a), &lda,
                      reinterpret_cast<dcomplex_t*> (tau),
                      reinterpret_cast<dcomplex_t*> (work), &info );
      }
      
    } 

    template <typename A, typename Tau, typename Work>
    inline
    int gerqf (A& a, Tau& tau, Work& work) {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<A>::matrix_structure, 
        traits::general_t
      >::value)); 
#endif 

      int const m = traits::matrix_size1 (a);
      int const n = traits::matrix_size2 (a);
      int const o = std::min(m,n);
      assert (o <= traits::vector_size (tau)); 
      assert (m <= traits::vector_size (work)); 

      int info; 
      if ( o >= detail::gerqfXOver )
          detail::gerqf (m, n,
                         traits::matrix_storage (a), 
                         traits::leading_dimension (a),
                         traits::vector_storage (tau),  
                         traits::vector_storage (work),
                         traits::vector_size (work),
                         info);
      else
          detail::gerq2 (m, n,
                         traits::matrix_storage (a), 
                         traits::leading_dimension (a),
                         traits::vector_storage (tau),  
                         traits::vector_storage (work),
                         info);
      return info; 
    }

    // Computation of the RQ factorization.
    // Workspace is allocated dynamically so that the optimization of
    // blas 3 calls is optimal.
    template <typename A, typename Tau>
    inline
    int gerqf (A& a, Tau& tau, optimal_workspace ) {
       typedef typename A::value_type value_type ;
       const int m = traits::matrix_size1 (a);
       traits::detail::array<value_type> work(std::max(1, m*32));
       return gerqf( a, tau, work );
    }

    // Computation of the RQ factorization.
    // Workspace is allocated dynamically to its minimum size.
    // Blas 3 calls are not optimal.
    template <typename A, typename Tau>
    inline
    int gerqf (A& a, Tau& tau, minimal_workspace ) {
       typedef typename A::value_type value_type ;
       const int m = traits::matrix_size1 (a);
       traits::detail::array<value_type> work(std::max(1, m));
       return gerqf( a, tau, work );
    }

    // Computation of the RQ factorization.
    // Workspace is taken from the array in workspace.
    // The calling sequence is
    // gerqf( a, tau, workspace( work ) ) where work is an array with the same value_type
    // as a.
    template <typename A, typename Tau, typename Work>
    inline
    int gerqf (A& a, Tau& tau, detail::workspace1<Work> workspace ) {
       return gerqf( a, tau, workspace.w_ );
    }

    // Function without workarray as argument
    template <typename A, typename Tau>
    inline
    int gerqf (A& a, Tau& tau) {
       return gerqf( a, tau, optimal_workspace() );
    }

  }

}}}

#endif 
