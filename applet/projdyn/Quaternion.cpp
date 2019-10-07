//****************************************************
//* quaternion.c++                                   *
//*                                                  *
//* Implementaion for a generalized quaternion class *
//*                                                  *
//* Written 1.25.00 by Angela Bennett                *
//****************************************************


#include "Quaternion.h"

//Quaternion
// -default constructor
// -creates a new quaternion with all parts equal to zero
template<class _Tp>
Quaternion<_Tp>::Quaternion(void)
{
    x = 0;
    y = 0;
    z = 0;
    w = 0;
}


//Quaternion
// -constructor
// -parametes : x, y, z, w elements of the quaternion
// -creates a new quaternion based on the elements passed in
template<class _Tp>
Quaternion<_Tp>::Quaternion(_Tp wi, _Tp xi, _Tp yi, _Tp zi)
{
    w = wi;
    x = xi;
    y = yi;
    z = zi;
}


//Quaternion
// -constructor
// -parameters : vector/array of four elements
// -creates a new quaternion based on the elements passed in
template<class _Tp>
Quaternion<_Tp>::Quaternion(_Tp v[4])
{
    w = v[0];
    x = v[1];
    y = v[2];
    z = v[3];
}


//Quaternion
// -copy constructor
// -parameters : const quaternion q
// -creates a new quaternion based on the quaternion passed in
template<class _Tp>
Quaternion<_Tp>::Quaternion(const Quaternion<_Tp>& q)
{
    w = q.w;
    x = q.x;
    y = q.y;
    z = q.z;
}

//~Quaternion
// -destructor
// -deleted dynamically allocated memory
template<class _Tp>
Quaternion<_Tp>::~Quaternion()
{
}


//operator=
// -parameters : q1 - Quaternion object
// -return value : Quaternion
// -when called on quaternion q2 sets q2 to be an object of  q3
template<class _Tp>
Quaternion<_Tp> Quaternion<_Tp>::operator = (const Quaternion<_Tp>& q)
{
    w = q.w;
    x = q.x;
    y = q.y;
    z = q.z;

    return (*this);
}

//operator+
// -parameters : q1 - Quaternion object
// -return value : Quaternion
// -when called on quaternion q2 adds q1 + q2 and returns the sum in a new quaternion
template<class _Tp>
Quaternion<_Tp> Quaternion<_Tp>::operator + (const Quaternion<_Tp>& q)
{
    return Quaternion(w+q.w, x+q.x, y+q.y, z+q.z);
}

//operator-
// -parameters : q1- Quaternion object
// -return values : Quaternion
// -when called on q1 subtracts q1 - q2 and returns the difference as a new quaternion
template<class _Tp>
Quaternion<_Tp> Quaternion<_Tp>::operator - (const Quaternion<_Tp>& q)
{
    return Quaternion(w-q.w, x-q.x, y-q.y, z-q.z);
}


//operator*
// -parameters : q1 - Quaternion object
// -return values : Quaternion
// -when called on a quaternion q2, multiplies q2 *q1  and returns the product in a new quaternion
template<class _Tp>
Quaternion<_Tp> Quaternion<_Tp>::operator * (const Quaternion<_Tp>& q)
{
    return Quaternion(
            w*q.w - x*q.x - y*q.y - z*q.z,
            w*q.x + x*q.w + y*q.z - z*q.y,
            w*q.y + y*q.w + z*q.x - x*q.z,
            w*q.z + z*q.w + x*q.y - y*q.x);
}

//operator/
// -parameters : q1 and q2- Quaternion objects
// -return values : Quaternion
// -divide q1 by q2 and returns the quotient q1
template<class _Tp>
Quaternion<_Tp> Quaternion<_Tp>::operator / (Quaternion<_Tp>& q)
{
    return ((*this) * (q.inverse()));
}


//operator+=
// -parameters : q1- Quaternion object
// -return values : Quaternion
// -when called on quaternion q3, adds q1 and q3 and returns the sum as q3
template<class _Tp>
Quaternion<_Tp>& Quaternion<_Tp>::operator += (const Quaternion<_Tp>& q)
{
    w += q.w;
    x += q.x;
    y += q.y;
    z += q.z;

    return (*this);
}


//operator-=
// -parameters : q1- Quaternion object
// -return values : Quaternion
// -when called on quaternion q3, subtracts q1 from q3 and returns the difference as q3
template<class _Tp>
Quaternion<_Tp>& Quaternion<_Tp>::operator -= (const Quaternion<_Tp>& q)
{
    w -= q.w;
    x -= q.x;
    y -= q.y;
    z -= q.z;

    return (*this);
}


//operator*=
// -parameters : q1- Quaternion object
// -return values : Quaternion
// -when called on quaternion q3, multiplies q3 by q1 and returns the product as q3
template<class _Tp>
Quaternion<_Tp>& Quaternion<_Tp>::operator *= (const Quaternion<_Tp>& q)
{
    _Tp w_val = w*q.w - x*q.x - y*q.y - z*q.z;
    _Tp x_val = w*q.x + x*q.w + y*q.z - z*q.y;
    _Tp y_val = w*q.y + y*q.w + z*q.x - x*q.z;
    _Tp z_val = w*q.z + z*q.w + x*q.y - y*q.x;

    w = w_val;
    x = x_val;
    y = y_val;
    z = z_val;

    return (*this);
}


//operator/=
// -parameters : q1- Quaternion object
// -return values : quaternion
// -when called on quaternion q3, divides q3 by q1 and returns the quotient as q3
template<class _Tp>
Quaternion<_Tp>& Quaternion<_Tp>::operator /= (Quaternion<_Tp>& q)
{
    (*this) = (*this)*q.inverse();
    return (*this);
}


//operator!=
// -parameters : q1 and q2- Quaternion objects
// -return value : bool
// -determines if q1 and q2 are not equal
template<class _Tp>
bool Quaternion<_Tp>::operator != (const Quaternion<_Tp>& q)
{
    return (w!=q.w || x!=q.x || y!=q.y || z!=q.z) ? true : false;
}

//operator==
// -parameters : q1 and q2- Quaternion objects
// -return value : bool
// -determines if q1 and q2 are equal
template<class _Tp>
bool Quaternion<_Tp>::operator == (const Quaternion<_Tp>& q)
{
    return (w==q.w && x==q.x && y==q.y && z==q.z) ? true : false;
}

//norm
// -parameters : none
// -return value : _Tp
// -when called on a quaternion object q, returns the norm of q
template<class _Tp>
_Tp Quaternion<_Tp>::norm()
{
    return (w*w + x*x + y*y + z*z);
}

//magnitude
// -parameters : none
// -return value : _Tp
// -when called on a quaternion object q, returns the magnitude q
template<class _Tp>
_Tp Quaternion<_Tp>::magnitude()
{
    return sqrt(norm());
}

//scale
// -parameters :  s- a value to scale q1 by
// -return value: quaternion
// -returns the original quaternion with each part, w,x,y,z, multiplied by some scalar s
template<class _Tp>
Quaternion<_Tp> Quaternion<_Tp>::scale(_Tp  s)
{
    return Quaternion(w*s, x*s, y*s, z*s);
}

// -parameters : none
// -return value : quaternion
// -when called on a quaternion object q, returns the inverse of q
template<class _Tp>
Quaternion<_Tp> Quaternion<_Tp>::inverse()
{
    return conjugate().scale(1/norm());
}

//conjugate
// -parameters : none
// -return value : quaternion
// -when called on a quaternion object q, returns the conjugate of q
template<class _Tp>
Quaternion<_Tp> Quaternion<_Tp>::conjugate()
{
    return Quaternion(w, -x, -y, -z);
}

//UnitQuaternion
// -parameters : none
// -return value : quaternion
// -when called on quaterion q, takes q and returns the unit quaternion of q
template<class _Tp>
Quaternion<_Tp> Quaternion<_Tp>::UnitQuaternion()
{
    return (*this).scale(1/(*this).magnitude());
}

// -parameters : vector of type _Tp
// -return value : void
// -when given a 3D vector, v, rotates v by this quaternion
template<class _Tp>
void Quaternion<_Tp>::QuatRotation(_Tp v[3])
{
    Quaternion <_Tp> qv(0, v[0], v[1], v[2]);
    Quaternion <_Tp> qm = (*this) * qv * (*this).inverse();

    v[0] = qm.x;
    v[1] = qm.y;
    v[2] = qm.z;
}
