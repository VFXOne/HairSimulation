//
// Created by lucas on 06.10.19.
//

#ifndef APPLET_QUATERNION_H
#define APPLET_QUATERNION_H


#include <math.h>

template<class _Tp>
class Quaternion
{

public:

    //Quaternion
    // -default constructor
    // -creates a new quaternion with all parts equal to zero
    Quaternion(void);

    //Quaternion
    // -constructor
    // -parametes : w, x, y, z elements of the quaternion
    // -creates a new quaternion based on the elements passed in
    Quaternion(_Tp wi, _Tp xi, _Tp yi, _Tp zi);

    //Quaternion
    // -constructor
    // -parameters : 4D vector
    // -creates a new quaternion based on the elements passed in
    Quaternion(_Tp v[4]);

    //Quaternion
    // -copy constructor
    // -parameters : const quaternion q
    // -creates a new quaternion based on the quaternion passed in
    Quaternion(const Quaternion<_Tp>& q);

    //~Quaternion
    // -default destructor
    ~Quaternion();

    //operator=
    // -parameters : q1- Quaternion object
    // -return values : Quaternion
    // -when called on quaternion q2 sets q2 to be an object of  q3
    Quaternion<_Tp> operator = (const Quaternion<_Tp>& q);

    //operator+
    // -parameters : q1 - Quaternion object
    // -return value : Quaternion
    // -when called on quaternion q2 adds q1 + q2 and returns the sum in a new quaternion
    Quaternion<_Tp> operator + (const Quaternion<_Tp>& q);

    //operator-
    // -parameters : q1- Quaternion object
    // -return values : Quaternion
    // -when called on q1 subtracts q1 - q2 and returns the difference as a new quaternion
    Quaternion<_Tp> operator - (const Quaternion<_Tp>& q);

    //operator*
    // -parameters : q1 - Quaternion object
    // -return values : Quaternion
    // -when called on a quaternion q2, multiplies q2 *q1  and returns the product in a new quaternion
    Quaternion<_Tp> operator * (const Quaternion<_Tp>& q);

    //operator/
    // -parameters : q1 and q2- Quaternion objects
    // -return values : Quaternion
    // -divide q1 by q2 and returns the quotient as q1
    Quaternion<_Tp> operator / (Quaternion<_Tp>& q);

    //operator+=
    // -parameters : q1- Quaternion object
    // -return values : Quaternion
    // -when called on quaternion q3 adds q1 and q3 and returns the sum as q3
    Quaternion<_Tp>& operator += (const Quaternion<_Tp>& q);

    //operator-=
    // -parameters : q1- Quaternion object
    // -return values : Quaternion
    // -when called on quaternion q3, subtracts q1 from q3 and returns the difference as q3
    Quaternion<_Tp>& operator -= (const Quaternion<_Tp>& q);

    //operator*=
    // -parameters : q1- Quaternion object
    // -return values : Quaternion
    // -when called on quaternion q3, multiplies q3 by q1 and returns the product as q3
    Quaternion<_Tp>& operator *= (const Quaternion<_Tp>& q);

    //operator/=
    // -parameters : q1- Quaternion object
    // -return values : quaternion
    // -when called on quaternion q3, divides q3 by q1 and returns the quotient as q3
    Quaternion<_Tp>& operator /= (Quaternion<_Tp>& q);

    /*//operator<<
    // -parameters : ostream o, quaternion q
    // -return values :
    // -prints out a quaternion by it's components
    friend inline ostream& operator << (ostream& output, const Quaternion<_Tp>& q)
    {
        output << "[" << q.w << ", " << "(" << q.x << ", " << q.y << ", " << q.z << ")]";
        return output;
    }*/

    //operator!=
    // -parameters : q1 and q2- Quaternion objects
    // -return value : bool
    // -determines if q1 and q2 and equal
    bool operator != (const Quaternion<_Tp>& q);

    //operator==
    // -parameters : q1 and q2- Quaternion objects
    // -return value : bool
    // -determines if q1 and q2 and equal
    bool operator == (const Quaternion<_Tp>& q);

    //other methods: norm, inverse, conjugate, toEuler

    //norm
    // -parameters : none
    // -return value : _Tp
    // -when called on a quaternion object q, returns the norm of q
    _Tp norm();

    //magnitude
    // -parameters : none
    // -return value : _Tp
    // -when called on a quaternion object q, returns the magnitude q
    _Tp magnitude();

    //scale
    // -parameters :  s- a value to scale q1 by
    // -return value: quaternion
    // -returns the original quaternion with each part, w,x,y,z, multiplied by some scalar s
    Quaternion<_Tp> scale(_Tp s);

    //inverse
    // -parameters : none
    // -return value : quaternion
    // -when called on a quaternion object q, returns the inverse of q
    Quaternion<_Tp> inverse();

    //conjugate
    // -parameters : none
    // -return value : quaternion
    // -when called on a quaternion object q, returns the conjugate of q
    Quaternion<_Tp> conjugate();

    //UnitQuaternion
    // -parameters : none
    // -return value : quaternion
    // -when called on quaterion q, takes q and returns the unit quaternion of q
    Quaternion<_Tp> UnitQuaternion();

    // -parameters : 3D vector of type _Tp
    // -return value : void
    // -when given a  3D vector, v, rotates v by the quaternion
    void QuatRotation(_Tp v[3]);


private:
    // [w, (x, y, z)]
    _Tp w, x, y, z;

};


#endif //APPLET_QUATERNION_H
