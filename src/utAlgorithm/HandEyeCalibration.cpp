/*
 * Ubitrack - Library for Ubiquitous Tracking
 * Copyright 2006, Technische Universitaet Muenchen, and individual
 * contributors as indicated by the @authors tag. See the
 * copyright.txt in the distribution for a full listing of individual
 * contributors.
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA, or see the FSF site: http://www.fsf.org.
 */

/**
 * @ingroup tracking_algorithms
 * @file
 * Implementation of  HandEye Calibration
 *
 * @author Daniel Muhra <muhra@in.tum.de>
 */


#include "HandEyeCalibration.h"

#ifdef HAVE_LAPACK

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/lapack/gels.hpp>

#include <log4cpp/Category.hh>
#include <utUtil/Exception.h>
#include <utUtil/Logging.h>
#include <utMath/MatrixOperations.h>

//shortcuts to namespaces
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

namespace Ubitrack { namespace Algorithm {


//helper class, which is used to mainly store 2 vectors.
template< typename T >
class TransformCont
{
	private:
		std::vector<Math::Matrix< T, 4, 4 > > m_hg;
		std::vector<Math::Matrix< T, 4, 4 > > m_hc;
		
	public:
		void setNumber( const std::size_t i, bool bUseAllPairs )
		{
			const std::size_t p = bUseAllPairs ? i * i / 2 : i;
			m_hg.reserve( p );
			m_hc.reserve( p );
		}
		
		std::size_t getNumber() const
		{ return m_hg.size(); }
		
		void push_back_hg( const Math::Matrix< T, 4, 4 > &hg )
		{ m_hg.push_back( hg ); }
		
		void push_back_hc( const Math::Matrix< T, 4, 4 > &hc )
		{ m_hc.push_back( hc ); }
			
		const Math::Matrix< T, 4, 4 >& hg_at( const std::size_t i ) const
		{ return m_hg[ i ]; }
		
		const Math::Matrix< T, 4, 4 >& hc_at( const std::size_t i ) const
		{ return m_hc[ i ]; }
};


/** \internal */
template<typename T> 
Math::Vector< T, 3 > computeSidesTrans( const Math::Matrix< T, 4, 4 >& hgij, const Math::Matrix< T, 4, 4 >& hcij, Math::Matrix< T, 3, 3 > rcg, 
	Math::Matrix< T, 3, 3 > & leftT)
{

	Math::Matrix< T, 3, 3 > rgij = ublas::subrange(hgij, 0, 3, 0, 3 );
	Math::Vector< T, 3 > tgij = ublas::subrange( ublas::column(hgij, 3), 0, 3 );
	Math::Vector< T, 3 > tcij = ublas::subrange( ublas::column(hcij, 3), 0, 3 );

	leftT = rgij - Math::Matrix< T, 3, 3 >::identity();

	return (ublas::prod(rcg, tcij) - tgij);
}


template<typename T> 
Math::Vector< T, 3 > computeTcg( TransformCont< T >& tc, const Math::Matrix< T, 3, 3 >& rcg)
{
	Math::Matrix< T, 0, 0 > tA(3*tc.getNumber(), 3);
	Math::Matrix< T, 0, 0 > tB(3*tc.getNumber(), 1);

	Math::Matrix< T, 4, 4 > hgij;
	Math::Matrix< T, 4, 4 > hcij;
	Math::Matrix< T, 3, 3 > leftT;

	Math::Vector< T, 3 > rightT;

	for(int i=0;i<tc.getNumber();i++)
	{
		hgij = tc.hg_at(i);
		hcij = tc.hc_at(i);

		rightT = computeSidesTrans(hgij, hcij, rcg, leftT);
		
		//right side
		tB(3*i, 0) = rightT(0);
		tB(3*i+1, 0) = rightT(1);
		tB(3*i+2, 0) = rightT(2);

		//left side

		tA(3*i, 0)	 = leftT(0, 0);
		tA(3*i+1, 0) = leftT(1, 0);
		tA(3*i+2, 0) = leftT(2, 0);
		tA(3*i, 1)	 = leftT(0, 1);
		tA(3*i+1, 1) = leftT(1, 1);
		tA(3*i+2, 1) = leftT(2, 1);
		tA(3*i, 2)   = leftT(0, 2);
		tA(3*i+1, 2) = leftT(1, 2);
		tA(3*i+2, 2) = leftT(2, 2);
	}

	lapack::gels('N', tA, tB);			//pcg_ will be written to tB, solves tA*tcg = tB and tcg -> tB

	Math::Vector< T, 3 > tcg;

	tcg(0) = tB(0, 0);
	tcg(1) = tB(1, 0);
	tcg(2) = tB(2, 0);

	return tcg;	
}


template< typename T > 
Math::Matrix< T, 3, 3 > skew(const Math::Vector< T, 3 >& rotVec)
{
	Math::Matrix< T, 3, 3 > skew;
	skew(0, 0) = 0.0;
	skew(0, 1) = - rotVec(2);
	skew(0, 2) = rotVec(1);
	skew(1, 0) = rotVec(2);
	skew(1, 1) = 0.0;
	skew(1, 2) = - rotVec(0);
	skew(2, 0) = - rotVec(1);
	skew(2, 1) = rotVec(0);
	skew(2, 2) = 0.0;
	return skew;
}


template< typename T > 
Math::Matrix< T, 3, 3 > getMatrix( const Math::Vector< T, 3 >& source)
{   
	T v[3];	
	v[0] = source(0) * source(0);
	v[1] = source(1) * source(1);
	v[2] = source(2) * source(2);

	T length = (v[0] + v[1] + v[2]);										//|Pr|²
	T a = (T)1.0 - length/(T)2.0;											//(1 - |Pr|²/2)
	
	// ublas::identity_matrix<T> id(3);
	Math::Matrix< T, 3, 3 > identity = Math::Matrix< T, 3, 3 >::identity() * a;

	Math::Matrix< T, 3, 3 > skewT = skew(source);

	T alpha = std::sqrt(4 - length);
	
	skewT *= T (alpha);
	
	Math::Matrix< T, 3, 3 > right;

	right = outer_prod(source, source) + skewT;

	right *= T (0.5);

	return (identity + right);	
} 


template< typename T > 
Math::Matrix< T, 3, 3 > getRcg(const Math::Vector< T, 3 >& pcg_)
{
	double v[3];
	double twice[3];

	Math::Vector< T, 3 > pcg;

	for(int i=0; i<3; i++)
	{
		v[i] = pcg_(i) * pcg_(i); 
		twice[i] = 2*pcg_(i);
	}

	double divisor = std::sqrt( 1.0 + (v[0] + v[1] + v[2]));

	for(int i=0; i<3; i++)
	{
		pcg(i) = T (twice[i]/divisor);
	}	

	return getMatrix(pcg);
}


template< typename T > 
Math::Vector< T, 3 > getQuaternion(const Math::Matrix< T, 3, 3 >& source)
{
    T quat[4];

	T q[4];
    T qoff[6];
    std::size_t c( 0 );

    q[0] = (1 + source(0, 0) + source(1, 1) + source(2, 2)) / 4;
    q[1] = (1 + source(0, 0) - source(1, 1) - source(2, 2)) / 4;
	q[2] = (1 - source(0, 0) + source(1, 1) - source(2, 2)) / 4;
	q[3] = (1 - source(0, 0) - source(1, 1) + source(2, 2)) / 4;
    if (q[c] < q[1])
        c = 1;   
    if (q[c] < q[2])
        c = 2;
    if (q[c] < q[3])
        c = 3;


    qoff[0] = (source(2, 1) - source(1, 2)) / 4;
    qoff[1] = (source(0, 2) - source(2, 0)) / 4;
    qoff[2] = (source(1, 0) - source(0, 1)) / 4;
    qoff[3] = (source(1, 0) + source(0, 1)) / 4;
    qoff[4] = (source(0, 2) + source(2, 0)) / 4;
    qoff[5] = (source(2, 1) + source(1, 2)) / 4;

    if (c==0)
    {
        quat[3] = std::sqrt(q[c]);
        quat[0] = qoff[0] / quat[3];
        quat[1] = qoff[1] / quat[3];
        quat[2] = qoff[2] / quat[3];
    } else if (c==1)
    {
        quat[0] = std::sqrt(q[c]);
        quat[3] = qoff[0] / quat[0];
        quat[1] = qoff[3] / quat[0];
        quat[2] = qoff[4] / quat[0];
    } else if (c==2)
    {
        quat[1] = std::sqrt(q[c]);
        quat[3] = qoff[1] / quat[1];
        quat[0] = qoff[3] / quat[1];
        quat[2] = qoff[5] / quat[1];
    } else if (c==3)
    {
        quat[2] = std::sqrt(q[c]);
        quat[3] = qoff[2] / quat[2];
        quat[0] = qoff[4] / quat[2];
        quat[1] = qoff[5] / quat[2];
    }

    // ensure w to be positive
    if (quat[3] < 0.0)
    {
        quat[0] = -quat[0];
        quat[1] = -quat[1];
        quat[2] = -quat[2];
        quat[3] = -quat[3];
    }

	Math::Vector< T, 3 > destination;

	destination(0) = quat[0];
	destination(1) = quat[1];
	destination(2) = quat[2];
	//destination(3) = 0.0;						//will be dropped

	return destination;
} 


template< typename T > 
Math::Vector< T, 3 > computeSidesRot( const Math::Matrix< T, 4, 4 >& hgij, const Math::Matrix< T, 4, 4 >& hcij, Math::Matrix< T, 3, 3 >& skewP)	//computes  P'cg
{
	Math::Matrix< T, 3, 3 > source;
	Math::Vector< T, 3> pgij;
	Math::Vector< T, 3 > pcij;

	//copy rgij to source
	source = ublas::subrange(hgij, 0, 3, 0, 3);

	//convert rgij to pgij
	pgij = getQuaternion(source);

	//copy rcij to source
	source = ublas::subrange(hcij, 0, 3, 0, 3);

	//convert rcij to pcij
	pcij = getQuaternion(source);

	Math::Vector< T, 3 > temp = (pgij + pcij);

	skewP = skew(temp);
	return (pcij - pgij);	
}


template< typename T > 
Math::Matrix< T, 3, 3 > computePcg( TransformCont< T >& tc)
{
	Math::Matrix< T, 0, 0 > tA(3*tc.getNumber(), 3);
	Math::Matrix< T, 0, 0 > tB(3*tc.getNumber(), 1);

	Math::Matrix< T, 4, 4 > hgij;
	Math::Matrix< T, 4, 4 > hcij;
	Math::Matrix< T, 3, 3 > skew;
	
	Math::Vector< T, 3 > rightR;

	for(int i=0;i<tc.getNumber();i++)
	{
		hgij = tc.hg_at(i);
		hcij = tc.hc_at(i);

		rightR = computeSidesRot(hgij, hcij, skew);

		//right side
		tB(3*i, 0) = rightR(0);
		tB(3*i+1, 0) = rightR(1);
		tB(3*i+2, 0) = rightR(2);

		//left side
		tA(3*i, 0)	 = skew(0, 0);
		tA(3*i+1, 0) = skew(1, 0);
		tA(3*i+2, 0) = skew(2, 0);
		tA(3*i, 1)	 = skew(0, 1);
		tA(3*i+1, 1) = skew(1, 1);
		tA(3*i+2, 1) = skew(2, 1);
		tA(3*i, 2)   = skew(0, 2);
		tA(3*i+1, 2) = skew(1, 2);
		tA(3*i+2, 2) = skew(2, 2);
	}

	lapack::gels('N', tA, tB);			//pcg_ will be written to tB, solves tA*Pcg_ = tB and Pcg_ -> tB

	Math::Vector< T, 3 > pcg_;

	pcg_(0) = tB(0, 0);
	pcg_(1) = tB(1, 0);
	pcg_(2) = tB(2, 0);

	return getRcg(pcg_);
}


template<typename T>
Math::Matrix< T, 4, 4 > computeTransformation(const Math::Matrix< T, 4, 4 >& hi, const Math::Matrix< T, 4, 4 >& hj, int mode) //mode: 0 to compute Hgij, 1 to compute Hcij
{
	Math::Matrix< T, 4, 4 > inverted;

	if(mode == 0) //Hgij
	{
		inverted = Math::invert_matrix(hj);				//invert matrix need right function
		return ublas::prod(inverted, hi);
	}
	else //Hcij
	{
		inverted = Math::invert_matrix(hi);				//invert matrix need right function
		return ublas::prod(hj, inverted);
	}
}


template<typename T> 
void fillTransformationVectors( TransformCont<T>& tc, const std::vector<Math::Matrix< T, 4, 4 > >& hand, const std::vector<Math::Matrix< T, 4, 4 > >& eye, bool bUseAllPairs )
{
	Math::Matrix< T, 4, 4 > hij;
	Math::Matrix< T, 4, 4 > hi;
	Math::Matrix< T, 4, 4 > hj;
	const std::size_t n_hands( hand.size() ) ;
	for( std::size_t i( 0 ); i<(n_hands-1); ++i )
	{
		const std::size_t to = bUseAllPairs ? n_hands : i + 2;
		for( std::size_t k( i+1 ); k<to; ++k )
		{
			//first we compute the hg
			hi = hand.at(i);
			hj = hand.at(k);

			hij = computeTransformation(hi, hj, 0);
			tc.push_back_hg(hij);

			//second we compute the hc
			hi = eye.at(i);
			hj = eye.at(k);

			hij = computeTransformation(hi, hj, 1);
			tc.push_back_hc(hij);
		}
	}
}


template<typename T> 
Math::Pose performHandEyeCalibrationImp ( const std::vector<Math::Matrix< T, 4, 4 > >& hand,  const std::vector<Math::Matrix< T, 4, 4 > >& eye, bool bUseAllPairs )
{
	static log4cpp::Category& logger(log4cpp::Category::getInstance( "Ubitrack.Calibration.HandEyeCalibration" )); 
	const std::size_t n_eyes( eye.size() );
	if( n_eyes != hand.size())
	{
		LOG4CPP_ERROR ( logger, "Input sizes of the vectors do not match ");
		UBITRACK_THROW ( "Input sizes do not match" );		
	}

	if( n_eyes <= 2 )
	{
		Math::Vector< T, 3 >* v = new Math::Vector< T, 3 >( 0, 0, 0 );
		Math::Quaternion* q = new Math::Quaternion();
		return Math::Pose( *q, *v);
	}
	
	TransformCont< T > tc;

	tc.setNumber( n_eyes, bUseAllPairs );

	fillTransformationVectors( tc, hand, eye, bUseAllPairs );			//readies values //ai = eye, bi = hand
	Math::Matrix< T, 3, 3 > rcg = computePcg< T >( tc );		//returns Rcg
	Math::Vector< T, 3 > tcg = computeTcg( tc, rcg );

	return Math::Pose(Math::Quaternion(rcg), tcg);
}


Math::Pose performHandEyeCalibration ( const std::vector< Math::Matrix< float, 4, 4 > >& hand,  const std::vector< Math::Matrix< float, 4, 4 > >& eye, bool bUseAllPairs )
{
	return performHandEyeCalibrationImp (hand, eye, bUseAllPairs );
}


Math::Pose performHandEyeCalibration ( const std::vector< Math::Matrix< double, 4, 4 > >& hand,  const std::vector< Math::Matrix< double, 4, 4 > >& eye, bool bUseAllPairs )
{
	return performHandEyeCalibrationImp (hand, eye, bUseAllPairs );
}


void fillTransformationVectors( TransformCont<double>& tc, const std::vector<Math::Pose>& hand, const std::vector<Math::Pose>& eye, bool bUseAllPairs )
{
	Math::Matrix< double, 4, 4 > hij;
	Math::Matrix< double, 4, 4 > hi;
	Math::Matrix< double, 4, 4 > hj;
	const std::size_t n_hands( hand.size() );
	for( std::size_t i( 0 ) ; i<(n_hands-1); ++i )
	{
		const std::size_t to = bUseAllPairs ? n_hands : i + 2;
		for( std::size_t k( i+1 ); k<to; ++k )
		{
			//first we compute the hg
			hi = Math::Matrix< double, 4, 4 >(hand.at(i));
			hj = Math::Matrix< double, 4, 4 >(hand.at(k));

			hij = computeTransformation(hi, hj, 0);
			tc.push_back_hg(hij);

			//second we compute the hc
			hi =  Math::Matrix< double, 4, 4 >(eye.at(i));
			hj =  Math::Matrix< double, 4, 4 >(eye.at(k));

			hij = computeTransformation( hi, hj, 1 );
			tc.push_back_hc(hij);
		}
	}
}


Math::Pose performHandEyeCalibration ( const std::vector< Math::Pose >& hand,  const std::vector< Math::Pose >& eye, bool bUseAllPairs )
{
	static log4cpp::Category& logger(log4cpp::Category::getInstance( "Ubitrack.Calibration.HandEyeCalibration" )); 
	const std::size_t n_eyes( eye.size() );
	if( n_eyes != hand.size())
	{
		LOG4CPP_ERROR ( logger, "Input sizes of the vectors do not match ");
		UBITRACK_THROW ( "Input sizes do not match" );		
	}

	if( n_eyes <= 2 )
	{
		Math::Vector< double, 3 >* v = new Math::Vector< double, 3 >( 0, 0, 0 );
		Math::Quaternion* q = new Math::Quaternion();
		return Math::Pose( *q, *v);
	}

	TransformCont< double > tc;

	tc.setNumber( n_eyes, bUseAllPairs );

	fillTransformationVectors( tc, hand, eye, bUseAllPairs );			//readies values //ai = eye, bi = hand
	Math::Matrix< double, 3, 3 > rcg = computePcg< double >( tc );									//returns also Rcg
	Math::Vector< double, 3 > tcg = computeTcg( tc, rcg );

	return Math::Pose(Math::Quaternion( rcg ), tcg);
}

}}

#endif // HAVE_LAPACK

