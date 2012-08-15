#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <Math/Quaternion.h>
#include <Math/Matrix.h>
#include "../tools.h"

#ifndef M_PI
#define _USE_MATH_DEFINES //for having PI
#include <math.h>
#endif

using namespace Ubitrack::Math;
namespace ublas = boost::numeric::ublas;


// check conversion quaternion <-> Euler angles 
void EulerCheck( Quaternion& rot0 ) {

	rot0.normalize();
	Vector<3> angles = rot0.getEulerAngles();

	double rx = angles[0];
	double ry = angles[1];
	double rz = angles[2];

	if (rx < 0) rx += 2*M_PI;
	if (ry < 0) ry += 2*M_PI;
	if (rz < 0) rz += 2*M_PI;

	Quaternion qx( Vector<3>(1,0,0), rx );
	Quaternion qy( Vector<3>(0,1,0), ry );
	Quaternion qz( Vector<3>(0,0,1), rz );

	Quaternion rot1 = qz*qy; rot1 = rot1*qx;
	rot1 = rot1.negateIfCloser( rot0 );

	BOOST_CHECK_CLOSE( rot0.x(), rot1.x(), 10e-6 );
	BOOST_CHECK_CLOSE( rot0.y(), rot1.y(), 10e-6 );
	BOOST_CHECK_CLOSE( rot0.z(), rot1.z(), 10e-6 );
	BOOST_CHECK_CLOSE( rot0.w(), rot1.w(), 10e-6 );
}

void TestQuaternionConversion() {

	// check 10 random rotations
	for (int i = 0; i < 10; i++) {
		Quaternion tmp( randomQuaternion() );
		
		EulerCheck( tmp );
	}

	// check pathological case: r_y == pi/2
	Quaternion y90( Vector<3>(0,1,0), M_PI/2.0 );
	EulerCheck( y90 );
}

