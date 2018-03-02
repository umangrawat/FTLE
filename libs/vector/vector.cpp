#include "vector.h"

Vector3 Vector3::getRoots(Vector3 a)
{
    Vector3 r;
    double c1 = a.y - a.z*a.z/3.0;
    double c0 = a.x - a.y*a.z/3.0 + 2.0/27.0*a.z*a.z*a.z;

    //Make cubic coefficient 4 and linear coefficient +- 3
    //by substituting y = z*k and multiplying with 4/k^3

    if (c1 == 0) {
        if (c0 == 0)     r.x =  0;
        else if (c0 > 0) r.x = -pow(c0, 1./3.);
        else             r.x =  pow(-c0, 1./3.);
    }
    else {
        bool negc1 = c1 < 0;
        double absc1 = negc1 ? -c1 : c1;
        double k = sqrt(4./3.*absc1);
        double d0 = c0 * 4./(k*k*k);

        //Find the first solution

        if (negc1) {
            if (d0 > 1)       r.x = -cosh(acosh(d0)/3);
            else if (d0 > -1) r.x = -cos(acos(d0)/3);
            else              r.x =  cosh(acosh(-d0)/3);
        }
        else {
            r.x = -sinh(asinh(d0)/3);
        }

        //Transform back
        r.x *= k;
    }
    r.x -= a.z/3;

    //Other two solutions
    double p = r.x + a.z;
    double q = r.x * p + a.y;

    double discrim = p*p - 4*q;
    //if (forceReal && discrim < 0.0) discrim = 0.0;

    if (discrim >= 0) {
        double root = sqrt(discrim);
        r.y = (-p - root)/2.;
        r.z = (-p + root)/2.;
        return 3;
    }
    else {
        double root = sqrt(-discrim);
        r.y = -p/2;
        r.z = root/2.;
        return 1;
    }
    return v;
}
