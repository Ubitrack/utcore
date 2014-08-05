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


#include "ErrorPose.h"
#include "ErrorVector.h"
#include "Stochastic/CovarianceTransform.h"

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>


namespace ublas = boost::numeric::ublas;

namespace Ubitrack { namespace Math {

namespace {

/**
 * \internal
 * computes the jacobian of pose inversion
 * See Daniel Pustka's Diplomarbeit, page 62 ff. for details.
 */
void inversionJacobian( Matrix< double, 6, 6 >& j, const Pose& p )
{
	double t5364 = p.rotation().x();
	double t5362 = p.rotation().w();
	double t5369 = p.rotation().y();
	double t5371 = p.rotation().z();
	double t5378 = p.translation()( 0 );
	double t5380 = p.translation()( 1 );
	double t5365 = t5364 * t5364;
	double t5386 = p.translation()( 2 );
	double t5382 = t5369 * t5378;
	double t5363 = t5362 * t5362;
	double t5388 = t5369 * t5369;
	double t5379 = 2. * t5364 * t5371 * t5378;
	double t5381 = 2. * t5369 * t5371 * t5380;
	double t5383 =  - ( t5364 * t5380 );
	double t5384 = t5382 + t5383;
	double t5385 = 2. * t5362 * t5384;
	double t5387 =  - 2. * t5365 * t5386;
	double t5389 =  - 2. * t5386 * t5388;
	double t5390 = t5379 + t5381 + t5385 + t5386 + t5387 + t5389;
	double t5394 = t5371 * t5371;
	double t5392 =  - 2. * t5362 * t5371 * t5378;
	double t5393 =  - 2. * t5365 * t5380;
	double t5395 =  - 2. * t5380 * t5394;
	double t5396 = 2. * t5369 * t5371 * t5386;
	double t5397 = t5362 * t5386;
	double t5398 = t5382 + t5397;
	double t5399 = 2. * t5364 * t5398;
	double t5400 = t5380 + t5392 + t5393 + t5395 + t5396 + t5399;
	double t5417 = t5371 * t5380;
	double t5418 =  - ( t5369 * t5386 );
	double t5419 = t5417 + t5418;
	double t5420 = t5362 * t5419;
	double t5421 = t5369 * t5380;
	double t5422 = t5371 * t5386;
	double t5423 = t5421 + t5422;
	double t5424 = t5364 * t5423;
	double t5425 = t5420 + t5424;
	double t5428 = t5362 * t5369;
	double t5429 = t5364 * t5371;
	double t5430 = t5428 + t5429;
	double t5370 = t5364 * t5369;
	double t5372 = t5362 * t5371;
	double t5373 = t5370 + t5372;
	double t5446 =  - 2. * t5394;
	double t5408 = t5362 * t5364;
	double t5409 = t5369 * t5371;
	double t5410 = t5408 + t5409;
	double t5453 =  - 2. * t5365;
	double t5445 =  - 2. * t5388;
	j( 0, 0 ) = 1. - 2. * ( t5363 + t5365 );
	j( 0, 1 ) =  - 2. * t5373;
	j( 0, 2 ) = 2. * t5362 * t5369 - 2. * t5364 * t5371;
	j( 0, 3 ) = 0.;
	j( 0, 4 ) = 2. * t5390;
	j( 0, 5 ) =  - 2. * t5400;
	j( 1, 0 ) =  - 2. * t5364 * t5369 + 2. * t5362 * t5371;
	j( 1, 1 ) = 1. - 2. * ( t5363 + t5388 );
	j( 1, 2 ) =  - 2. * t5410;
	j( 1, 3 ) =  - 2. * t5390;
	j( 1, 4 ) = 0.;
	j( 1, 5 ) = t5378 * ( 2. - 4. * t5388 - 4. * t5394 ) + 4. * t5425;
	j( 2, 0 ) =  - 2. * t5430;
	j( 2, 1 ) = 2. * t5362 * t5364 - 2. * t5369 * t5371;
	j( 2, 2 ) = 1. - 2. * ( t5363 + t5394 );
	j( 2, 3 ) = 2. * t5400;
	j( 2, 4 ) = t5378 * (  - 2. + 4. * t5388 + 4. * t5394 ) - 4. * t5425;
	j( 2, 5 ) = 0.;
	j( 3, 0 ) = 0.;
	j( 3, 1 ) = 0.;
	j( 3, 2 ) = 0.;
	j( 3, 3 ) = 1. + t5445 + t5446;
	j( 3, 4 ) = 2. * t5364 * t5369 - 2. * t5362 * t5371;
	j( 3, 5 ) = 2. * t5430;
	j( 4, 0 ) = 0.;
	j( 4, 1 ) = 0.;
	j( 4, 2 ) = 0.;
	j( 4, 3 ) = 2. * t5373;
	j( 4, 4 ) = 1. + t5446 + t5453;
	j( 4, 5 ) =  - 2. * t5362 * t5364 + 2. * t5369 * t5371;
	j( 5, 0 ) = 0.;
	j( 5, 1 ) = 0.;
	j( 5, 2 ) = 0.;
	j( 5, 3 ) =  - 2. * t5362 * t5369 + 2. * t5364 * t5371;
	j( 5, 4 ) = 2. * t5410;
	j( 5, 5 ) = 1. + t5445 + t5453;
}


void multiplicationJacobians( Matrix< double, 6, 6 >& j1, Matrix< double, 6, 6 >& j2, const Pose& p1, const Pose& p2 )
{
	double t5559 = p1.rotation().y();
	double t5557 = p2.translation()( 1 );
	double t5556 = p1.rotation().z();
	double t5560 = p2.translation()( 2 );
	double t5564 = p1.rotation().w();
	double t5555 = p1.rotation().x();
	double t5573 = p2.translation()( 0 );
	double t5576 = t5555 * t5555;
	double t5588 = t5556 * t5556;
	double t5589 = 2. * t5588;
	double t5566 = t5556 * t5560;
	double t5575 = t5564 * t5564;
	double t5586 = t5559 * t5559;
	double t5582 = t5555 * t5559;
	double t5609 = 4. * t5575;
	double t5594 = t5556 * t5559;
	double t5610 = 4. * t5586;
	double t5571 = t5555 * t5556;
	double t5605 = t5555 * t5573;
	double t5565 = t5557 * t5559;
	double t5638 = p2.rotation().y();
	double t5641 = p2.rotation().z();
	double t5647 = p2.rotation().w();
	double t5645 = p2.rotation().x();
	double t5642 = t5641 * t5641;
	double t5643 =  - 2. * t5642;
	double t5657 = t5645 * t5645;
	double t5658 =  - 2. * t5657;
	double t5639 = t5638 * t5638;
	double t5640 =  - 2. * t5639;
	double t5577 = t5575 + t5576;
	double t5578 = 2. * t5577;
	double t5579 =  - 1. + t5578;
	double t5570 = t5559 * t5564;
	double t5572 = t5570 + t5571;
	double t5613 = t5556 * t5564;
	double t5614 = t5582 + t5613;
	double t5620 = t5555 * t5564;
	double t5621 = t5594 + t5620;
	
	j1( 0, 0 ) = 1.;
	j1( 0, 1 ) = 0.;
	j1( 0, 2 ) = 0.;
	j1( 0, 3 ) = 4. * t5555 * ( t5556 * t5557 - t5559 * t5560 ) + 4. * t5564 * ( t5565 + t5566 );
	j1( 0, 4 ) =  - 4. * t5572 * t5573 + 2. * t5560 * t5579;
	j1( 0, 5 ) = 4. * t5573 * (  - ( t5556 * t5564 ) + t5582 ) + 2. * t5557 * (  - 1. + 2. * t5586 + t5589 );
	j1( 1, 0 ) = 0.;
	j1( 1, 1 ) = 1.;
	j1( 1, 2 ) = 0.;
	j1( 1, 3 ) = 2. * t5560 * (  - 1. + 2. * t5576 + t5589 ) + 4. * t5557 * (  - ( t5555 * t5564 ) + t5594 );
	j1( 1, 4 ) = 4. * t5559 * ( t5555 * t5560 - t5556 * t5573 ) + 4. * t5564 * ( t5566 + t5605 );
	j1( 1, 5 ) = t5573 * (  - 2. + t5609 + t5610 ) - 4. * t5557 * t5614;
	j1( 2, 0 ) = 0.;
	j1( 2, 1 ) = 0.;
	j1( 2, 2 ) = 1.;
	j1( 2, 3 ) = t5557 * (  - 2. + 4. * t5588 + t5609 ) - 4. * t5560 * t5621;
	j1( 2, 4 ) = 4. * t5560 * (  - ( t5559 * t5564 ) + t5571 ) + t5573 * (  - 2. + 4. * t5576 + t5610 );
	j1( 2, 5 ) = 4. * t5556 * (  - ( t5555 * t5557 ) + t5559 * t5573 ) + 4. * t5564 * ( t5565 + t5605 );
	j1( 3, 0 ) = 0.;
	j1( 3, 1 ) = 0.;
	j1( 3, 2 ) = 0.;
	j1( 3, 3 ) = 1. + t5640 + t5643;
	j1( 3, 4 ) = 2. * ( t5638 * t5645 + t5641 * t5647 );
	j1( 3, 5 ) = 2. * t5641 * t5645 - 2. * t5638 * t5647;
	j1( 4, 0 ) = 0.;
	j1( 4, 1 ) = 0.;
	j1( 4, 2 ) = 0.;
	j1( 4, 3 ) = 2. * t5638 * t5645 - 2. * t5641 * t5647;
	j1( 4, 4 ) = 1. + t5643 + t5658;
	j1( 4, 5 ) = 2. * ( t5638 * t5641 + t5645 * t5647 );
	j1( 5, 0 ) = 0.;
	j1( 5, 1 ) = 0.;
	j1( 5, 2 ) = 0.;
	j1( 5, 3 ) = 2. * ( t5641 * t5645 + t5638 * t5647 );
	j1( 5, 4 ) = 2. * t5638 * t5641 - 2. * t5645 * t5647;
	j1( 5, 5 ) = 1. + t5640 + t5658;
	
	j2( 0, 0 ) = t5579;
	j2( 0, 1 ) = 2. * t5555 * t5559 - 2. * t5556 * t5564;
	j2( 0, 2 ) = 2. * t5572;
	j2( 0, 3 ) = 0.;
	j2( 0, 4 ) = 0.;
	j2( 0, 5 ) = 0.;
	j2( 1, 0 ) = 2. * t5614;
	j2( 1, 1 ) =  - 1. + 2. * ( t5575 + t5586 );
	j2( 1, 2 ) = 2. * t5556 * t5559 - 2. * t5555 * t5564;
	j2( 1, 3 ) = 0.;
	j2( 1, 4 ) = 0.;
	j2( 1, 5 ) = 0.;
	j2( 2, 0 ) = 2. * t5555 * t5556 - 2. * t5559 * t5564;
	j2( 2, 1 ) = 2. * t5621;
	j2( 2, 2 ) =  - 1. + 2. * ( t5575 + t5588 );
	j2( 2, 3 ) = 0.;
	j2( 2, 4 ) = 0.;
	j2( 2, 5 ) = 0.;
	j2( 3, 0 ) = 0.;
	j2( 3, 1 ) = 0.;
	j2( 3, 2 ) = 0.;
	j2( 3, 3 ) = 1.;
	j2( 3, 4 ) = 0.;
	j2( 3, 5 ) = 0.;
	j2( 4, 0 ) = 0.;
	j2( 4, 1 ) = 0.;
	j2( 4, 2 ) = 0.;
	j2( 4, 3 ) = 0.;
	j2( 4, 4 ) = 1.;
	j2( 4, 5 ) = 0.;
	j2( 5, 0 ) = 0.;
	j2( 5, 1 ) = 0.;
	j2( 5, 2 ) = 0.;
	j2( 5, 3 ) = 0.;
	j2( 5, 4 ) = 0.;
	j2( 5, 5 ) = 1.;
}

template< class MT >
void errorPoseTimesVectorJacobian( MT& j, const Pose& p, const Vector< double, 3 >& v )
{
	typedef typename MT::value_type VType;
	VType t216 = p.rotation().y();
	VType t214 = v( 1 );
	VType t213 = p.rotation().z();
	VType t217 = v( 2 );
	VType t221 = p.rotation().w();
	VType t212 = p.rotation().x();
	VType t228 = v( 0 );
	VType t231 = t221 * t221;
	VType t233 = t212 * t212;
	VType t235 = t216 * t216;
	VType t236 = t213 * t213;
	VType t237 = t235 + t236;
	VType t223 = t213 * t217;
	VType t264 =  - ( t228 * t231 );
	VType t265 = t228 * t233;
	VType t259 = t212 * t228;
	VType t222 = t214 * t216;
	j( 0, 0 ) = 1;
	j( 0, 1 ) = 0;
	j( 0, 2 ) = 0;
	j( 0, 3 ) = 4 * ( t212 * ( t213 * t214 - t216 * t217 ) + t221 * ( t222 + t223 ) );
	j( 0, 4 ) = 2 * (  - 2 * t212 * t213 * t228 - 2 * t216 * t221 * t228 + t217 * t231 + t217 * t233 - t217 * t237 );
	j( 0, 5 ) = 4 * t212 * t216 * t228 - 4 * t213 * t221 * t228 - 2 * t214 * t231 - 2 * t214 * t233 + 2 * t214 * t237;
	j( 1, 0 ) = 0;
	j( 1, 1 ) = 1;
	j( 1, 2 ) = 0;
	j( 1, 3 ) = 2 * ( 2 * t213 * t214 * t216 - 2 * t212 * t214 * t221 - t217 * t231 - t217 * t235 + t217 * ( t233 + t236 ) );
	j( 1, 4 ) = 4 * ( t216 * ( t212 * t217 - t213 * t228 ) + t221 * ( t223 + t259 ) );
	j( 1, 5 ) =  - 2 * ( 2 * t212 * t214 * t216 + 2 * t213 * t214 * t221 + t228 * (  - t235 + t236 ) + t264 + t265 );
	j( 2, 0 ) = 0;
	j( 2, 1 ) = 0;
	j( 2, 2 ) = 1;
	j( 2, 3 ) =  - 2 * ( 2 * t213 * t216 * t217 + 2 * t212 * t217 * t221 - t214 * t231 + t214 * t233 + t214 * t235 - t214 * t236 );
	j( 2, 4 ) = 2 * ( 2 * t212 * t213 * t217 - 2 * t216 * t217 * t221 + t228 * ( t235 - t236 ) + t264 + t265 );
	j( 2, 5 ) = 4 * ( t213 * (  - ( t212 * t214 ) + t216 * t228 ) + t221 * ( t222 + t259 ) );
}

void invertMultiplyJacobians( Matrix< double, 6, 6 >& j1, Matrix< double, 6, 6 >& j2, const Pose& p1, const Pose& p2 )
{
	double t701 = p2.rotation().y();
	double t704 = p2.rotation().z();
	double t710 = p2.rotation().w();
	double t708 = p2.rotation().x();
	double t705 = t704 * t704;
	double t706 =  - 2 * t705;
	double t720 = t708 * t708;
	double t721 =  - 2 * t720;
	double t702 = t701 * t701;
	double t703 =  - 2 * t702;
	j1( 0, 0 ) = 1 + t703 + t706;
	j1( 0, 1 ) = 2 * ( t701 * t708 + t704 * t710 );
	j1( 0, 2 ) = 2 * t704 * t708 - 2 * t701 * t710;
	j1( 0, 3 ) = 0.;
	j1( 0, 4 ) = 0.;
	j1( 0, 5 ) = 0.;
	j1( 1, 0 ) = 2 * t701 * t708 - 2 * t704 * t710;
	j1( 1, 1 ) = 1 + t706 + t721;
	j1( 1, 2 ) = 2 * ( t701 * t704 + t708 * t710 );
	j1( 1, 3 ) = 0.;
	j1( 1, 4 ) = 0.;
	j1( 1, 5 ) = 0.;
	j1( 2, 0 ) = 2 * ( t704 * t708 + t701 * t710 );
	j1( 2, 1 ) = 2 * t701 * t704 - 2 * t708 * t710;
	j1( 2, 2 ) = 1 + t703 + t721;
	j1( 2, 3 ) = 0.;
	j1( 2, 4 ) = 0.;
	j1( 2, 5 ) = 0.;
	j1( 3, 0 ) = 0.;
	j1( 3, 1 ) = 0.;
	j1( 3, 2 ) = 0.;
	j1( 3, 3 ) = 1.;
	j1( 3, 4 ) = 0.;
	j1( 3, 5 ) = 0.;
	j1( 4, 0 ) = 0.;
	j1( 4, 1 ) = 0.;
	j1( 4, 2 ) = 0.;
	j1( 4, 3 ) = 0.;
	j1( 4, 4 ) = 1.;
	j1( 4, 5 ) = 0.;
	j1( 5, 0 ) = 0.;
	j1( 5, 1 ) = 0.;
	j1( 5, 2 ) = 0.;
	j1( 5, 3 ) = 0.;
	j1( 5, 4 ) = 0.;
	j1( 5, 5 ) = 1.;

	double t788 = p2.rotation().y();
	double t791 = p2.rotation().z();
	double t797 = p2.rotation().w();
	double t795 = p2.rotation().x();
	double t804 = p2.translation()( 0 );
	double t805 = p1.translation()( 0 );
	double t806 =  - t805;
	double t807 = t804 + t806;
	double t810 = p2.translation()( 1 );
	double t812 = p1.translation()( 1 );
	double t818 = p2.translation()( 2 );
	double t820 = p1.translation()( 2 );
	double t817 = t795 * t795;
	double t824 =  - t812;
	double t825 = t810 + t824;
	double t827 =  - t818;
	double t828 = t820 + t827;
	double t792 = t791 * t791;
	double t793 = 2 * t792;
	double t834 =  - 4 * t804;
	double t835 = 4 * t805;
	double t836 = t834 + t835;
	double t837 = t788 * t836;
	double t864 =  - t820;
	double t865 = t818 + t864;
	double t789 = t788 * t788;
	double t857 = 4 * t795 * t825;
	double t838 = 4 * t797 * t828;
	double t850 = 2 * t817;
	double t790 = 2 * t789;
	double t860 =  - 4 * t810;
	double t861 = 4 * t812;
	double t862 = t860 + t861;
	double t809 = 4 * t788 * t807;
	double t863 = t791 * t862;
	double t866 = 4 * t788 * t865;
	double t889 = 4 * t797 * t865;
	double t906 = p1.rotation().y();
	double t907 = t906 * t906;
	double t909 = p1.rotation().z();
	double t910 = t909 * t909;
	double t912 =  - 4 * t907;
	double t913 =  - 4 * t910;
	double t914 = 2 + t912 + t913;
	double t919 = p1.rotation().w();
	double t917 = p1.rotation().x();
	double t918 = t906 * t917;
	double t920 = t909 * t919;
	double t921 = t918 + t920;
	double t923 = t906 * t919;
	double t924 =  - ( t909 * t917 );
	double t925 = t923 + t924;
	double t908 = 2 * t907;
	double t911 = 2 * t910;
	double t946 = 4 * t907;
	double t947 = 4 * t910;
	double t948 =  - 2 + t946 + t947;
	double t949 = t788 * t948;
	double t794 =  - 1 + t790 + t793;
	double t976 = t917 * t919;
	double t977 = t906 * t909;
	double t978 = t976 + t977;
	double t980 = t917 * t917;
	double t981 = 4 * t980;
	double t982 =  - 2 + t947 + t981;
	double t992 =  - 4 * t980;
	double t993 = 2 + t913 + t992;
	double t973 =  - ( t909 * t919 );
	double t974 = t918 + t973;
	double t991 = 2 * t980;
	double t1000 =  - ( t906 * t917 );
	double t1001 = t1000 + t920;
	double t938 = t909 * t917;
	double t1025 =  - ( t906 * t909 );
	double t1026 = t1025 + t976;
	double t1023 =  - 2 + t946 + t981;
	double t1021 = t923 + t938;
	double t851 =  - 1 + t793 + t850;
	double t963 =  - 2 * t907;
	double t1030 = 2 + t912 + t992;
	j2( 0, 0 ) = t794;
	j2( 0, 1 ) =  - 2 * ( t788 * t795 + t791 * t797 );
	j2( 0, 2 ) =  - 2 * t791 * t795 + 2 * t788 * t797;
	j2( 0, 3 ) = 0.;
	j2( 0, 4 ) = 4 * t791 * t795 * t807 + t797 * ( t809 + 4 * t795 * (  - t810 + t812 ) ) + t817 * (  - 4 * t818 + 4 * t820 ) + t788 * ( 4 * t791 * t825 + 4 * t788 * t828 );
	j2( 0, 5 ) = 4 * t817 * t825 + 4 * t791 * ( t797 * t807 + t791 * t825 + t788 * t828 ) + t795 * ( t837 + t838 );
	j2( 1, 0 ) =  - 2 * t788 * t795 + 2 * t791 * t797;
	j2( 1, 1 ) = t851;
	j2( 1, 2 ) =  - 2 * ( t788 * t791 + t795 * t797 );
	j2( 1, 3 ) = t791 * t795 * t836 + t797 * ( t837 + t857 ) + 4 * t817 * t865 + t788 * ( t863 + t866 );
	j2( 1, 4 ) = 0.;
	j2( 1, 5 ) = t789 * t836 + t788 * ( t838 + t857 ) + t791 * ( 4 * t797 * t825 + t791 * t836 + 4 * t795 * t865 );
	j2( 2, 0 ) =  - 2 * ( t791 * t795 + t788 * t797 );
	j2( 2, 1 ) =  - 2 * t788 * t791 + 2 * t795 * t797;
	j2( 2, 2 ) =  - 1 + t790 + t850;
	j2( 2, 3 ) = t817 * t862 + t791 * ( t797 * t836 + t863 + t866 ) + t795 * ( t809 + t889 );
	j2( 2, 4 ) = 4 * t789 * t807 + t791 * ( 4 * t791 * t807 + 4 * t795 * t828 + t797 * t862 ) + t788 * ( t795 * t862 + t889 );
	j2( 2, 5 ) = 0.;
	j2( 3, 0 ) = 0.;
	j2( 3, 1 ) = 0.;
	j2( 3, 2 ) = 0.;
	j2( 3, 3 ) =  - 1 + t908 + t911 + t789 * t914 + t792 * t914 + 4 * t791 * (  - ( t797 * t921 ) + t795 * t925 ) - 4 * t788 * ( t795 * t921 + t797 * t925 );
	j2( 3, 4 ) = 4 * t817 * t921 - 2 * ( t791 * t797 * (  - 1 + t908 + t911 ) - (  - 1 + t793 ) * t921 + 2 * t788 * t791 * (  - ( t906 * t919 ) + t938 ) ) + t795 * ( 4 * t797 * t925 + t949 );
	j2( 3, 5 ) = t797 * ( 4 * t795 * t921 + t949 ) - 2 * ( ( 2 * t788 * t791 * t906 + t909 - 2 * t789 * t909 ) * t917 + ( (  - 1 + t790 ) * t906 + 2 * t788 * t791 * t909 ) * t919 + 2 * t817 * t925 + t791 * t795 * ( 1 - 2 * t910 + t963 ) );
	j2( 4, 0 ) = 0.;
	j2( 4, 1 ) = 0.;
	j2( 4, 2 ) = 0.;
	j2( 4, 3 ) = 2 * t794 * t974 + t795 * (  - 4 * t791 * t978 + t788 * t982 ) + t797 * ( 4 * t788 * t978 + t791 * t982 );
	j2( 4, 4 ) =  - 1 + t911 - 4 * t791 * ( t1001 * t797 + t788 * t978 ) - 4 * t795 * ( t788 * t974 + t797 * t978 ) + t991 + t792 * t993 + t817 * t993;
	j2( 4, 5 ) =  - 2 * t906 * t909 + 4 * t817 * t906 * t909 - 4 * t791 * t795 * t906 * t917 + 4 * t791 * t795 * t909 * t919 - 2 * t917 * t919 + 4 * t817 * t917 * t919 + 4 * t789 * t978 + 2 * t788 * t791 * (  - 1 + t911 + t991 ) + t797 * ( 4 * t1001 * t788 + t795 * t993 );
	j2( 5, 0 ) = 0.;
	j2( 5, 1 ) = 0.;
	j2( 5, 2 ) = 0.;
	j2( 5, 3 ) = 2 * t1021 * t794 + ( 4 * t1026 * t788 + t1023 * t791 ) * t795 + ( t1030 * t788 + 4 * t1026 * t791 ) * t797;
	j2( 5, 4 ) = ( 4 * t1021 * t791 + t1023 * t795 ) * t797 - 2 * ( t1026 * t851 + t788 * ( 2 * t1021 * t795 + t791 * ( 1 + t963 - 2 * t980 ) ) );
	j2( 5, 5 ) =  - 1 + t1030 * t789 + 4 * t788 * ( t1026 * t791 - t1021 * t797 ) - 4 * t795 * ( t1021 * t791 + t1026 * t797 ) + t1030 * t817 + t908 + t991;
}


/**
 * Helper function for creating the additive<->multiplicative error conversion jacobian
 */
struct ErrorConversion
{
	Math::Quaternion r;

	ErrorConversion( const Math::Quaternion& q )
		: r( q )
	{}

	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& jacobian ) const
	{
		result = input;
		jacobian( 0, 0 ) = jacobian( 1, 1 ) = jacobian( 2, 2 ) = jacobian( 3, 3 ) = r.w();
		jacobian( 0, 1 ) = jacobian( 3, 2 ) = r.z();
		jacobian( 0, 2 ) = jacobian( 1, 3 ) = -r.y();
		jacobian( 0, 3 ) = jacobian( 2, 1 ) = -r.x();
		jacobian( 1, 0 ) = jacobian( 2, 3 ) = -r.z();
		jacobian( 1, 2 ) = jacobian( 3, 0 ) = r.x();
		jacobian( 2, 0 ) = jacobian( 3, 1 ) = r.y();
	}

};

} // anonymous namespace


ErrorPose ErrorPose::operator~() const
{
	// covariance transform
	Matrix< double, 6, 6 > jacobian;
	Matrix< double, 6, 6 > tmp;
	Matrix< double, 6, 6 > newCovariance;
	inversionJacobian( jacobian, *this );
	noalias( tmp ) = ublas::prod( jacobian, m_covariance );
	noalias( newCovariance ) = ublas::prod( tmp, ublas::trans( jacobian ) );

	return ErrorPose( Pose::operator~(), newCovariance );
}


ErrorPose operator*( const ErrorPose& a, const ErrorPose& b )
{
	// covariance transform
	Matrix< double, 6, 6 > jacobian1;
	Matrix< double, 6, 6 > jacobian2;
	Matrix< double, 6, 6 > tmp;
	Matrix< double, 6, 6 > newCovariance;
	multiplicationJacobians( jacobian1, jacobian2, a, b );
	
	// here we could save some time by regarding the zero and identity blocks of the jacobians...
	noalias( tmp ) = ublas::prod( jacobian1, a.covariance() );
	noalias( newCovariance ) = ublas::prod( tmp, ublas::trans( jacobian1 ) );
	noalias( tmp ) = ublas::prod( jacobian2, b.covariance() );
	noalias( newCovariance ) += ublas::prod( tmp, ublas::trans( jacobian2 ) );

	return ErrorPose( static_cast< const Pose& >( a ) * static_cast< const Pose& >( b ), newCovariance );
}


ErrorPose operator*( const Pose& a, const ErrorPose& b )
{
	// covariance transform
	Matrix< double, 6, 6 > jacobian1;
	Matrix< double, 6, 6 > jacobian2;
	Matrix< double, 6, 6 > tmp;
	Matrix< double, 6, 6 > newCovariance;
	multiplicationJacobians( jacobian1, jacobian2, a, b );
	
	// here we could save some time by regarding the zero and identity blocks of the jacobians...
	noalias( tmp ) = ublas::prod( jacobian2, b.covariance() );
	noalias( newCovariance ) = ublas::prod( tmp, ublas::trans( jacobian2 ) );

	return ErrorPose( static_cast< const Pose& >( a ) * static_cast< const Pose& >( b ), newCovariance );
}


ErrorPose operator*( const ErrorPose& a, const Pose& b )
{
	// covariance transform
	Matrix< double, 6, 6 > jacobian1;
	Matrix< double, 6, 6 > jacobian2;
	Matrix< double, 6, 6 > tmp;
	Matrix< double, 6, 6 > newCovariance;
	multiplicationJacobians( jacobian1, jacobian2, a, b );

	// here we could save some time by regarding the zero and identity blocks of the jacobians...
	noalias( tmp ) = ublas::prod( jacobian1, a.covariance() );
	noalias( newCovariance ) = ublas::prod( tmp, ublas::trans( jacobian1 ) );

	return ErrorPose( static_cast< const Pose& >( a ) * static_cast< const Pose& >( b ), newCovariance );
}


ErrorVector< double, 3 > operator*( const ErrorPose& a, const Math::Vector< double, 3 >& b )
{
	// covariance transform
	Matrix< double, 3, 6 > jacobian;
	errorPoseTimesVectorJacobian( jacobian, a, b );

	Matrix< double, 3, 6 > tmp;
	Matrix< double, 3, 3 > newCovariance;

	noalias( tmp ) = ublas::prod( jacobian, a.covariance() );
	noalias( newCovariance ) = ublas::prod( tmp, ublas::trans( jacobian ) );

	return ErrorVector< double, 3 >( static_cast< const Pose& >( a ) * b, newCovariance );
}

//TODO

ErrorVector< double, 3 > operator*( const ErrorPose& a, const Math::ErrorVector< double, 3 >& b )
{
	
	// covariance transform
	Matrix< double, 3, 6 > jacobian;
	errorPoseTimesVectorJacobian( jacobian, a, b.value );

	Matrix< double, 3, 6 > tmp;
	Matrix< double, 3, 3 > newCovariance;

	noalias( tmp ) = ublas::prod( jacobian, a.covariance() );
	noalias( newCovariance ) = ublas::prod( tmp, ublas::trans( jacobian ) );

	return ErrorVector< double, 3 >( static_cast< const Pose& >( a ) * b.value, newCovariance );
}


ErrorPose invertMultiply( const ErrorPose& a, const ErrorPose& b )
{
	// covariance transform
	Matrix< double, 6, 6 > jacobian1;
	Matrix< double, 6, 6 > jacobian2;
	Matrix< double, 6, 6 > tmp;
	Matrix< double, 6, 6 > newCovariance;
	multiplicationJacobians( jacobian1, jacobian2, a, b );
	
	// here we could save some time by regarding the zero and identity blocks of the jacobians...
	noalias( tmp ) = ublas::prod( jacobian1, a.covariance() );
	noalias( newCovariance ) = ublas::prod( tmp, ublas::trans( jacobian1 ) );
	noalias( tmp ) = ublas::prod( jacobian2, b.covariance() );
	noalias( newCovariance ) += ublas::prod( tmp, ublas::trans( jacobian2 ) );

	return ErrorPose( ~static_cast< const Pose& >( a ) * static_cast< const Pose& >( b ), newCovariance );
}


void ErrorPose::toAdditiveErrorVector( ErrorVector< double, 7 >& v )
{
	Pose::toVector( v.value );

	ublas::subrange( v.covariance, 0, 6, 0, 6 ) = m_covariance;
	for ( unsigned i = 0; i < 7; i++ )
		v.covariance( 6, i ) = v.covariance( i, 6 ) = 0;
	// set the w variance to the sum of the x, y, and z rotation errors
	v.covariance( 6, 6 ) = v.covariance( 3, 3 ) + v.covariance( 4, 4 ) + v.covariance( 5, 5 );

	Stochastic::transformRangeInternalWithCovariance< 7 >( ErrorConversion( rotation() ), v, 3, 7, 3, 7 );
}


ErrorPose ErrorPose::fromAdditiveErrorVector( const ErrorVector< double, 7 >& v )
{
	Math::ErrorVector< double, 7 > p( v );
	Quaternion q( Quaternion::fromVector( ublas::subrange( v.value, 3, 7 ) ) );
	Stochastic::transformRangeInternalWithCovariance< 7 >( ErrorConversion( ~q ), p, 3, 7, 3, 7 );
	
	return ErrorPose( Pose::fromVector( p.value ), ublas::subrange( p.covariance, 0, 6, 0, 6 ) );
}

ErrorPose linearInterpolate( const ErrorPose& x, const ErrorPose& y, double t )
{
	double a = ( 1 - t ) * ( 1 - t );
	double b = t * t;
	return ErrorPose( 
		linearInterpolate( static_cast< const Pose& >( x ), static_cast< const Pose& >( y ), t ),
		a * x.covariance() + b * y.covariance() );
}

std::ostream& operator<<( std::ostream& s, const ErrorPose& ep )
{
	s << ep.translation() << " " << ep.rotation() << "\n" << ep.covariance();
	{	// remove this block, if too much info annoys you ;)
		const Math::Matrix< double, 6, 6 >& covar = ep.covariance();
		const double stdDevPos = std::sqrt ( covar ( 0, 0 ) + covar ( 1, 1 ) + covar ( 2, 2 ) );
		const double norm = std::sqrt ( covar ( 3, 3 ) + covar ( 4, 4 ) + covar ( 5, 5 ) );
		const double phiRad = std::asin ( norm ) * 2;
		const double phiDeg = phiRad * ( 180 / boost::math::constants::pi< double >() );
		s << "\nStd dev position: " << stdDevPos << " [m], Std dev orientation: " << phiRad << "/" << phiDeg << " [rad/deg]\n";
	}
	
	return s;
}

} } // namespaceUbitrack::Math

