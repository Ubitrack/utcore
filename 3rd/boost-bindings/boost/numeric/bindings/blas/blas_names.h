//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//  Copyright (C) 2002, 2003 Si-Lab b.v.b.a and Toon Knapen 

#ifndef BOOST_NUMERIC_BINDINGS_BLAS_BLAS_NAMES_H
#define BOOST_NUMERIC_BINDINGS_BLAS_BLAS_NAMES_H

#ifndef ANDROID
#include <boost/numeric/bindings/traits/fortran.h>

/*
#ifndef BOOST_NUMERIC_BINDINGS_USE_CLAPACK
#include <boost/numeric/bindings/traits/fortran.h>
#define FORTRAN_ID2( id ) FORTRAN_ID( id )
#else
#define FORTRAN_ID2( id ) f2c_##id
#endif 
*/


//
// level 1
//
#define BLAS_SSCAL FORTRAN_ID( sscal )
#define BLAS_DSCAL FORTRAN_ID( dscal )
#define BLAS_CSCAL FORTRAN_ID( cscal )
#define BLAS_ZSCAL FORTRAN_ID( zscal )

#define BLAS_SAXPY FORTRAN_ID( saxpy )
#define BLAS_DAXPY FORTRAN_ID( daxpy )
#define BLAS_CAXPY FORTRAN_ID( caxpy )
#define BLAS_ZAXPY FORTRAN_ID( zaxpy )

#define BLAS_SDOT  FORTRAN_ID( sdot )
#define BLAS_DDOT  FORTRAN_ID( ddot )

#define BLAS_CDOTU FORTRAN_ID( cdotu )
#define BLAS_ZDOTU FORTRAN_ID( zdotu )

#define BLAS_CDOTC FORTRAN_ID( cdotc )
#define BLAS_ZDOTC FORTRAN_ID( zdotc )

#define BLAS_SNRM2 FORTRAN_ID( snrm2 )
#define BLAS_DNRM2 FORTRAN_ID( dnrm2 )
#define BLAS_SCNRM2 FORTRAN_ID( scnrm2 )
#define BLAS_DZNRM2 FORTRAN_ID( dznrm2 )

#define BLAS_SASUM FORTRAN_ID( sasum )
#define BLAS_DASUM FORTRAN_ID( dasum )
#define BLAS_SCASUM FORTRAN_ID( scasum )
#define BLAS_DZASUM FORTRAN_ID( dzasum )

#define BLAS_SCOPY FORTRAN_ID( scopy )
#define BLAS_DCOPY FORTRAN_ID( dcopy )
#define BLAS_CCOPY FORTRAN_ID( ccopy )
#define BLAS_ZCOPY FORTRAN_ID( zcopy )

//
// level 2
//
#define BLAS_SGEMV FORTRAN_ID( sgemv )
#define BLAS_DGEMV FORTRAN_ID( dgemv )
#define BLAS_CGEMV FORTRAN_ID( cgemv )
#define BLAS_ZGEMV FORTRAN_ID( zgemv )

#define BLAS_SGER  FORTRAN_ID( sger )
#define BLAS_DGER  FORTRAN_ID( dger )

#define BLAS_CGERU FORTRAN_ID( cgeru )
#define BLAS_ZGERU FORTRAN_ID( zgeru )

#define BLAS_CGERC FORTRAN_ID( cgerc )
#define BLAS_ZGERC FORTRAN_ID( zgerc )

//
// level 3
//
#define BLAS_SGEMM FORTRAN_ID( sgemm )
#define BLAS_DGEMM FORTRAN_ID( dgemm )
#define BLAS_CGEMM FORTRAN_ID( cgemm )
#define BLAS_ZGEMM FORTRAN_ID( zgemm )

#define BLAS_SSYMM FORTRAN_ID( ssymm )
#define BLAS_DSYMM FORTRAN_ID( dsymm )
#define BLAS_CSYMM FORTRAN_ID( csymm )
#define BLAS_ZSYMM FORTRAN_ID( zsymm )

#define BLAS_SSYRK FORTRAN_ID( ssyrk )
#define BLAS_DSYRK FORTRAN_ID( dsyrk )
#define BLAS_CSYRK FORTRAN_ID( csyrk )
#define BLAS_ZSYRK FORTRAN_ID( zsyrk )
#define BLAS_CHERK FORTRAN_ID( cherk )
#define BLAS_ZHERK FORTRAN_ID( zherk )

#define BLAS_STRSM FORTRAN_ID( strsm )
#define BLAS_DTRSM FORTRAN_ID( dtrsm )
#define BLAS_CTRSM FORTRAN_ID( ctrsm )
#define BLAS_ZTRSM FORTRAN_ID( ztrsm )

#else



#ifndef BOOST_NUMERIC_BINDINGS_USE_CLAPACK
#include <boost/numeric/bindings/traits/fortran.h>
#define FORTRAN_ID2( id ) FORTRAN_ID( id )
#else
#define FORTRAN_ID2( id ) f2c_##id
#endif 



//
// level 1
//
#define BLAS_SSCAL FORTRAN_ID2( sscal )
#define BLAS_DSCAL FORTRAN_ID2( dscal )
#define BLAS_CSCAL FORTRAN_ID2( cscal )
#define BLAS_ZSCAL FORTRAN_ID2( zscal )

#define BLAS_SAXPY FORTRAN_ID2( saxpy )
#define BLAS_DAXPY FORTRAN_ID2( daxpy )
#define BLAS_CAXPY FORTRAN_ID2( caxpy )
#define BLAS_ZAXPY FORTRAN_ID2( zaxpy )

#define BLAS_SDOT  FORTRAN_ID2( sdot )
#define BLAS_DDOT  FORTRAN_ID2( ddot )

#define BLAS_CDOTU FORTRAN_ID2( cdotu )
#define BLAS_ZDOTU FORTRAN_ID2( zdotu )

#define BLAS_CDOTC FORTRAN_ID2( cdotc )
#define BLAS_ZDOTC FORTRAN_ID2( zdotc )

#define BLAS_SNRM2 FORTRAN_ID2( snrm2 )
#define BLAS_DNRM2 FORTRAN_ID2( dnrm2 )
#define BLAS_SCNRM2 FORTRAN_ID2( scnrm2 )
#define BLAS_DZNRM2 FORTRAN_ID2( dznrm2 )

#define BLAS_SASUM FORTRAN_ID2( sasum )
#define BLAS_DASUM FORTRAN_ID2( dasum )
#define BLAS_SCASUM FORTRAN_ID2( scasum )
#define BLAS_DZASUM FORTRAN_ID2( dzasum )

#define BLAS_SCOPY FORTRAN_ID2( scopy )
#define BLAS_DCOPY FORTRAN_ID2( dcopy )
#define BLAS_CCOPY FORTRAN_ID2( ccopy )
#define BLAS_ZCOPY FORTRAN_ID2( zcopy )

//
// level 2
//
#define BLAS_SGEMV FORTRAN_ID2( sgemv )
#define BLAS_DGEMV FORTRAN_ID2( dgemv )
#define BLAS_CGEMV FORTRAN_ID2( cgemv )
#define BLAS_ZGEMV FORTRAN_ID2( zgemv )

#define BLAS_SGER  FORTRAN_ID2( sger )
#define BLAS_DGER  FORTRAN_ID2( dger )

#define BLAS_CGERU FORTRAN_ID2( cgeru )
#define BLAS_ZGERU FORTRAN_ID2( zgeru )

#define BLAS_CGERC FORTRAN_ID2( cgerc )
#define BLAS_ZGERC FORTRAN_ID2( zgerc )

//
// level 3
//
#define BLAS_SGEMM FORTRAN_ID2( sgemm )
#define BLAS_DGEMM FORTRAN_ID2( dgemm )
#define BLAS_CGEMM FORTRAN_ID2( cgemm )
#define BLAS_ZGEMM FORTRAN_ID2( zgemm )

#define BLAS_SSYMM FORTRAN_ID2( ssymm )
#define BLAS_DSYMM FORTRAN_ID2( dsymm )
#define BLAS_CSYMM FORTRAN_ID2( csymm )
#define BLAS_ZSYMM FORTRAN_ID2( zsymm )

#define BLAS_SSYRK FORTRAN_ID2( ssyrk )
#define BLAS_DSYRK FORTRAN_ID2( dsyrk )
#define BLAS_CSYRK FORTRAN_ID2( csyrk )
#define BLAS_ZSYRK FORTRAN_ID2( zsyrk )
#define BLAS_CHERK FORTRAN_ID2( cherk )
#define BLAS_ZHERK FORTRAN_ID2( zherk )

#define BLAS_STRSM FORTRAN_ID2( strsm )
#define BLAS_DTRSM FORTRAN_ID2( dtrsm )
#define BLAS_CTRSM FORTRAN_ID2( ctrsm )
#define BLAS_ZTRSM FORTRAN_ID2( ztrsm )


#endif

#endif // BOOST_NUMERIC_BINDINGS_BLAS_BLAS_NAMES_H
