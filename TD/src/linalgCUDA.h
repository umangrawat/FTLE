
#ifndef __LINEALGCUDA_H__
#define __LINEALGCUDA_H__

#ifndef MAKEDEVICE
#define MAKEDEVICE __device__
#endif

#include <limits>
const double LIMIT_DOUBLE = 1000* std::numeric_limits<double>::epsilon();

////changed 3 to 2
/// Vec lib for GPU (device)
/*
MAKEDEVICE void  d_vec3sub( vec3 a,  vec3 b,  vec3 c)     { c[0] = a[0] - b[0]; c[1] = a[1] - b[1]; c[2] = a[2] - b[2]; }
MAKEDEVICE void d_vec3print(vec3 v, char* name) {printf("%c  %f  %f  %f \n", name,v[0],v[1],v[2]);}
MAKEDEVICE void  d_vec3set( vec3 a, double a0, double a1, double a2) { a[0] = a0; a[1] = a1; a[2] = a2; }

MAKEDEVICE  void  d_vec3copy( vec3 a,  vec3 b) { b[0] = a[0]; b[1] = a[1]; b[2] = a[2];}

MAKEDEVICE double d_vec3dot(vec3 a, vec3 b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

MAKEDEVICE double d_vec3sqr(vec3 a) { return d_vec3dot(a, a); }

MAKEDEVICE double d_vec3mag(vec3 a) { return sqrt(d_vec3sqr(a)); }

*/

MAKEDEVICE void  d_vec2sub( vec2 a,  vec2 b,  vec2 c)     { c[0] = a[0] - b[0]; c[1] = a[1] - b[1];}
MAKEDEVICE void d_vec2print(vec2 v, char* name) {printf("%c  %f  %f \n", name,v[0],v[1]);}
MAKEDEVICE void  d_vec2set( vec2 a, double a0, double a1) { a[0] = a0; a[1] = a1;}

MAKEDEVICE  void  d_vec2copy( vec2 a,  vec2 b) { b[0] = a[0]; b[1] = a[1];}

MAKEDEVICE double d_vec2dot(vec2 a, vec2 b) { return a[0] * b[0] + a[1] * b[1];}

MAKEDEVICE double d_vec2sqr(vec2 a) { return d_vec2dot(a, a); }

MAKEDEVICE double d_vec2mag(vec2 a) { return sqrt(d_vec2sqr(a)); }

////changed 3D to 2D
/*
MAKEDEVICE void d_vec3scal(vec3 a, double b, vec3 c) {
    c[0] = a[0] * b;
    c[1] = a[1] * b;
    c[2] = a[2] * b;
}

MAKEDEVICE void d_vec3add(vec3 a, vec3 b, vec3 c) {
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

MAKEDEVICE void d_vec3nrm(vec3 a, vec3 b) {
    double l = d_vec3mag(a);
    if (l == 0) l = 1;
    b[0] = a[0] / l;
    b[1] = a[1] / l;
    b[2] = a[2] / l;
}
MAKEDEVICE void  d_vec3zero( vec3 a) { a[0] = a[1] = a[2] = 0.0; }
*/

MAKEDEVICE void d_vec2scal(vec2 a, double b, vec2 c) {
    c[0] = a[0] * b;
    c[1] = a[1] * b;
    ////c[2] = a[2] * b;
}

MAKEDEVICE void d_vec2add(vec2 a, vec2 b, vec2 c) {
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    ////c[2] = a[2] + b[2];
}

MAKEDEVICE void d_vec2nrm(vec2 a, vec2 b) {
    double l = d_vec2mag(a);
    if (l == 0) l = 1;
    b[0] = a[0] / l;
    b[1] = a[1] / l;
    ////b[2] = a[2] / l;
}
MAKEDEVICE void  d_vec2zero( vec2 a) { a[0] = a[1] = 0.0; }

////changed 3D to 2D
/*

MAKEDEVICE void d_vec3lint(vec3 v0, vec3 v1, double t, vec3 v)
// Do a linear interpolation given the 2 corner values and the variable
{
    vec3 tmp;

    d_vec3zero(v);
    d_vec3scal(v0, (1-t), tmp);  d_vec3add(v, tmp, v);
    d_vec3scal(v1,   t  , tmp);  d_vec3add(v, tmp, v);
}

MAKEDEVICE void d_vec3bilint(vec3 v00, vec3 v10, vec3 v01, vec3 v11, double s, double t, vec3 v)
// Do a bilinear interpolation of a VECTOR in a quadrangle
// given the 4 corner values and the two parameters
{
    vec3 tmp;

    d_vec3zero(v);
    d_vec3scal(v00, (1-s)*(1-t), tmp);  d_vec3add(v, tmp, v);
    d_vec3scal(v10,   s  *(1-t), tmp);  d_vec3add(v, tmp, v);
    d_vec3scal(v01, (1-s)*  t  , tmp);  d_vec3add(v, tmp, v);
    d_vec3scal(v11,   s  *  t  , tmp);  d_vec3add(v, tmp, v);
}

MAKEDEVICE void d_vec3trilint(vec3 v000, vec3 v100, vec3 v010, vec3 v110, vec3 v001, vec3 v101, vec3 v011, vec3 v111,
                              double r, double s, double t, vec3 v)
// trilinear interpolation inside a cube
{
    vec3 w0, w1;
    d_vec3bilint(v000, v100, v010, v110, r, s, w0);
    d_vec3bilint(v001, v101, v011, v111, r, s, w1);
    d_vec3lint(w0, w1, t, v);
}

MAKEDEVICE void d_mat3setrows(mat3 a, vec3 a0, vec3 a1, vec3 a2)
{
    a[0][0] = a0[0]; a[0][1] = a0[1]; a[0][2] = a0[2];
    a[1][0] = a1[0]; a[1][1] = a1[1]; a[1][2] = a1[2];
    a[2][0] = a2[0]; a[2][1] = a2[1]; a[2][2] = a2[2];
}
MAKEDEVICE void  d_mat3copy( mat3 a,  mat3 b) { memcpy(b, a, sizeof( mat3)); }

MAKEDEVICE void  d_mat3mul( mat3 a,  mat3 b,  mat3 c)
{
    mat3 d;
    d[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
    d[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
    d[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
    d[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
    d[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
    d[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
    d[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
    d[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
    d[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
    d_mat3copy(d, c);
}
MAKEDEVICE double  d_mat3det( mat3 a)
{
    return  a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1]) +
            a[0][1]*(a[1][2]*a[2][0] - a[1][0]*a[2][2]) +
            a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0]);
}
*/

MAKEDEVICE void d_vec2lint(vec2 v0, vec2 v1, double t, vec2 v)
// Do a linear interpolation given the 2 corner values and the variable
{
    vec2 tmp;

    d_vec2zero(v);
    d_vec2scal(v0, (1-t), tmp);  d_vec2add(v, tmp, v);
    d_vec2scal(v1,   t  , tmp);  d_vec2add(v, tmp, v);
}

MAKEDEVICE void d_vec2bilint(vec2 v00, vec2 v10, vec2 v01, vec2 v11, double s, double t, vec2 v)
// Do a bilinear interpolation of a VECTOR in a quadrangle
// given the 4 corner values and the two parameters
{
    vec2 tmp;

    d_vec2zero(v);
    d_vec2scal(v00, (1-s)*(1-t), tmp);  d_vec2add(v, tmp, v);
    d_vec2scal(v10,   s  *(1-t), tmp);  d_vec2add(v, tmp, v);
    d_vec2scal(v01, (1-s)*  t  , tmp);  d_vec2add(v, tmp, v);
    d_vec2scal(v11,   s  *  t  , tmp);  d_vec2add(v, tmp, v);
}

MAKEDEVICE void d_vec2trilint(vec2 v000, vec2 v100, vec2 v010, vec2 v110, vec2 v001, vec2 v101, vec2 v011, vec2 v111,
                              double r, double s, double t, vec2 v)
// trilinear interpolation inside a cube
{
    vec2 w0, w1;
    d_vec2bilint(v000, v100, v010, v110, r, s, w0);
    d_vec2bilint(v001, v101, v011, v111, r, s, w1);
    d_vec2lint(w0, w1, t, v);
}

MAKEDEVICE void d_mat2setrows(mat2 a, vec2 a0, vec2 a1)
{
    a[0][0] = a0[0]; a[0][1] = a0[1];
    a[1][0] = a1[0]; a[1][1] = a1[1];
    ////a[2][0] = a2[0]; a[2][1] = a2[1]; a[2][2] = a2[2];
}
MAKEDEVICE void  d_mat2copy( mat2 a,  mat2 b) { memcpy(b, a, sizeof( mat2)); }

MAKEDEVICE void  d_mat2mul( mat2 a,  mat2 b,  mat2 c)
{
    mat2 d;
    d[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0];
    d[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1];
    ////d[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
    d[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0];
    d[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1];
    ////d[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
    ////d[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
    ////d[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
    ////d[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
    d_mat2copy(d, c);
}
MAKEDEVICE double  d_mat2det( mat2 a)
{
    return  a[0][0]*(a[1][1])-
            a[0][1]*(a[1][0]);
}

/*

MAKEDEVICE  void d_mat3invariants(mat3 m, vec3 pqr)
{
    // invariant0 = -det(M)
    pqr[0] = - d_mat3det(m);

    // invariant1 = det2(M#0) + det2(M#1) + det2(M#2)
    pqr[1] = m[1][1]*m[2][2] - m[1][2]*m[2][1]
             + m[2][2]*m[0][0] - m[2][0]*m[0][2]
             + m[0][0]*m[1][1] - m[0][1]*m[1][0];

    // invariant2 = -trace M
    pqr[2] = -(m[0][0] + m[1][1] + m[2][2]);
}
MAKEDEVICE int d_vec3cubicroots(vec3 a, vec3 r, bool forceReal=false)
//  Cubic equation (multiple solutions are returned several times)
//
//	Solves equation
//	    1 * x^3 + a[2]*x^2 + a[1]*x + a[0] = 0
//
//	On output,
//	    r[0], r[1], r[2], or
//	    r[0], r[1] +- i*r[2] are the roots
//
//	returns number of real solutions

{
    // Eliminate quadratic term by substituting
    // x = y - a[2] / 3

    double c1 = a[1] - a[2]*a[2]/3.;
    double c0 = a[0] - a[1]*a[2]/3. + 2./27.*a[2]*a[2]*a[2];

    // Make cubic coefficient 4 and linear coefficient +- 3
    // by substituting y = z*k and multiplying with 4/k^3

    if (c1 == 0) {
        if (c0 == 0)     r[0] =  0;
        else if (c0 > 0) r[0] = -pow(c0, 1./3.);
        else             r[0] =  pow(-c0, 1./3.);
    }
    else {
        bool negc1 = c1 < 0;
        double absc1 = negc1 ? -c1 : c1;

        double k = sqrt(4./3.*absc1);

        double d0 = c0 * 4./(k*k*k);

        // Find the first solution

        if (negc1) {
            if (d0 > 1)       r[0] = -cosh(acosh(d0)/3);
            else if (d0 > -1) r[0] = -cos(acos(d0)/3);
            else              r[0] =  cosh(acosh(-d0)/3);
        }
        else {
            r[0] = -sinh(asinh(d0)/3);
        }

        // Transform back
        r[0] *= k;
    }
    r[0] -= a[2]/3;

    // Other two solutions
    double p = r[0] + a[2];
    double q = r[0]*p + a[1];

    double discrim = p*p - 4*q;
    if (forceReal && discrim < 0.0) discrim = 0.0;

    if (discrim >= 0) {
        double root = sqrt(discrim);
        r[1] = (-p - root)/2.;
        r[2] = (-p + root)/2.;
        return 3;
    }
    else {
        double root = sqrt(-discrim);
        r[1] = -p/2;
        r[2] = root/2.;
        return 1;
    }
}


MAKEDEVICE int d_mat3eigenvalues(mat3 m, vec3 lambda)
// calculate eigenvalues in lambda, return number of real eigenvalues.
// either returnval==1, lambda[0]=real ev, l[1] real part+-l[2] imag part
// or     returnval==3, lambda[0-2] = eigenvalues
{
    vec3 pqr;
    d_mat3invariants(m, pqr);

    // force real solutions for symmetric matrices
    bool forceReal = false;
    if (m[1][0] == m[0][1] && m[2][0] == m[0][2] && m[2][1] == m[1][2])
        forceReal = true;

    return (d_vec3cubicroots(pqr, lambda, forceReal));
}

MAKEDEVICE  void  d_mat3trp( mat3 a,  mat3 b)
{
    if (a != b) d_mat3copy(a, b);
    double x;
    x = b[0][1]; b[0][1] = b[1][0]; b[1][0] = x;
    x = b[0][2]; b[0][2] = b[2][0]; b[2][0] = x;
    x = b[1][2]; b[1][2] = b[2][1]; b[2][1] = x;
}

#endif
*/

MAKEDEVICE  void d_mat2invariants(mat2 m, vec2 pqr)
{
    // invariant0 = -det(M)
    pqr[0] = - d_mat2det(m);

    //// invariant1 = det2(M#0) + det2(M#1) + det2(M#2)
    ////pqr[1] = m[0][0]*m[1][1] - m[0][1]*m[1][0];

    ////changed 2 to 1 in pqr[2]
    // invariant2 = -trace M
    pqr[1] = -(m[0][0] + m[1][1]);
}

////used 2D square roots from linalg.h
MAKEDEVICE int d_vec2squareroots(vec2 a, vec2 r)
/*
 *	Solves equation
 *	    1 * x^2 + a[1]*x + a[0] = 0
 *
 *	On output, 
 *	    r[0], r[1] or
 *	    r[0] +- i*r[1] are the roots 
 *	   
 *	returns number of real solutions
 */
{
  double discrim, root;
  
  discrim = a[1] * a[1] - 4 * a[0];
  
  if (discrim >= 0) {
    root = sqrt(discrim);
    r[0] = (-a[1] - root) / 2.0;
    r[1] = (-a[1] + root) / 2.0;
    return(2);
  }
  else {
    root = sqrt(-discrim);
    r[0] = -a[1] / 2.0;
    r[1] = root / 2.0;
    return(0);
  }
}

/*
MAKEDEVICE int d_vec3cubicroots(vec3 a, vec3 r, bool forceReal=false)
//  Cubic equation (multiple solutions are returned several times)
//
//	Solves equation
//	    1 * x^3 + a[2]*x^2 + a[1]*x + a[0] = 0
//
//	On output,
//	    r[0], r[1], r[2], or
//	    r[0], r[1] +- i*r[2] are the roots
//
//	returns number of real solutions

{
    // Eliminate quadratic term by substituting
    // x = y - a[2] / 3

    double c1 = a[1] - a[2]*a[2]/3.;
    double c0 = a[0] - a[1]*a[2]/3. + 2./27.*a[2]*a[2]*a[2];

    // Make cubic coefficient 4 and linear coefficient +- 3
    // by substituting y = z*k and multiplying with 4/k^3

    if (c1 == 0) {
        if (c0 == 0)     r[0] =  0;
        else if (c0 > 0) r[0] = -pow(c0, 1./3.);
        else             r[0] =  pow(-c0, 1./3.);
    }
    else {
        bool negc1 = c1 < 0;
        double absc1 = negc1 ? -c1 : c1;

        double k = sqrt(4./3.*absc1);

        double d0 = c0 * 4./(k*k*k);

        // Find the first solution

        if (negc1) {
            if (d0 > 1)       r[0] = -cosh(acosh(d0)/3);
            else if (d0 > -1) r[0] = -cos(acos(d0)/3);
            else              r[0] =  cosh(acosh(-d0)/3);
        }
        else {
            r[0] = -sinh(asinh(d0)/3);
        }

        // Transform back
        r[0] *= k;
    }
    r[0] -= a[2]/3;

    // Other two solutions
    double p = r[0] + a[2];
    double q = r[0]*p + a[1];

    double discrim = p*p - 4*q;
    if (forceReal && discrim < 0.0) discrim = 0.0;

    if (discrim >= 0) {
        double root = sqrt(discrim);
        r[1] = (-p - root)/2.;
        r[2] = (-p + root)/2.;
        return 3;
    }
    else {
        double root = sqrt(-discrim);
        r[1] = -p/2;
        r[2] = root/2.;
        return 1;
    }
}
*/





////changed 3 to 2
/*
MAKEDEVICE int d_mat2eigenvalues(mat2 m, vec2 lambda)

// calculate eigenvalues in lambda, return number of real eigenvalues.
// either returnval==1, lambda[0]=real ev, l[1] real part+-l[2] imag part
// or     returnval==3, lambda[0-2] = eigenvalues
{
    vec2 pqr;
    d_mat2invariants(m, pqr);

	
    // force real solutions for symmetric matrices
    bool forceReal = false;
    if (m[1][0] == m[0][1])
        forceReal = true;

    ////changed cubic to square
    return (d_vec2squareroots(pqr, lambda));
}
*/

MAKEDEVICE int d_mat2eigenvalues(mat2 m, vec2 lambda)
{
  vec2 pqr;

  d_mat2invariants(m, pqr);

  return (d_vec2squareroots(pqr, lambda));
}

MAKEDEVICE  void  d_mat2trp( mat2 a,  mat2 b)
{
    if (a != b) d_mat2copy(a, b);
    double x;
    x = b[0][1]; b[0][1] = b[1][0]; b[1][0] = x;
    ////x = b[0][2]; b[0][2] = b[2][0]; b[2][0] = x;
    ////x = b[1][2]; b[1][2] = b[2][1]; b[2][1] = x;
}

#endif
