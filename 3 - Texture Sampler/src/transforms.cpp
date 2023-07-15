#include "transforms.h"

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include <cmath>

#define MY_PI 3.1415926536

namespace CGL {

Vector2D operator*(const Matrix3x3 &m, const Vector2D &v) {
	Vector3D mv = m * Vector3D(v.x, v.y, 1);
	return Vector2D(mv.x / mv.z, mv.y / mv.z);
}

Matrix3x3 translate(float dx, float dy) {
	// Part 3: Fill this in.
	Matrix3x3 T = Matrix3x3::identity();
	T(0,2) = dx;
	T(1,2) = dy;
	return T;
}

Matrix3x3 scale(float sx, float sy) {
	// Part 3: Fill this in.
	Matrix3x3 S = Matrix3x3::identity();
	S(0,0) = sx;
	S(1,1) = sy;
	return S;
}

// The input argument is in degrees counterclockwise
Matrix3x3 rotate(float deg) {
	// Part 3: Fill this in.
	Matrix3x3 R = Matrix3x3::identity();
	R(0,0) = cos(MY_PI * deg/180);
	R(0,1) = -sin(MY_PI * deg/180);
	R(1,0) = sin(MY_PI * deg/180);
	R(1,1) = cos(MY_PI * deg/180);
	return R;
}

}
