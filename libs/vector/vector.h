#ifndef VECTOR_H
#define VECTOR_H

#pragma once

class Vector3;

class Vector4 {
public:
    double x;
    double y;
    double z;
    double w;

    Vector4();
    Vector4(const Vector4 &v);
    Vector4(double, double, double, double);
    Vector4(const Vector3 &v, double);
    Vector4 &operator=(const Vector4 &rhs); 
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      this = this + vector
    /// @param rhs                  vector to be added
    /// @return                     result of addition
    ///////////////////////////////////////////////////////////////////////////////////////////     
    Vector4 &operator=(const Vector3 &rhs); 
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      this = this + vector
    /// @param rhs                  vector to be added
    /// @return                     result of addition
    ///////////////////////////////////////////////////////////////////////////////////////////     
    Vector4 &operator+=(const Vector4 &rhs); 
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      this = this * scalar
    /// @param scalar               scalar value to be multiplied to vector
    /// @return                     result of multiplication with a scalar
    ///////////////////////////////////////////////////////////////////////////////////////////     
    Vector4 &operator*=(const double &scalar); 
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      negate vector
    /// @return                     negative vector (all components *(-1))
    ///////////////////////////////////////////////////////////////////////////////////////////     
    Vector4 operator-();
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      subtraction
    /// @return                     result of subtraction
    ///////////////////////////////////////////////////////////////////////////////////////////     
    Vector4 operator-(Vector4);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      addition
    /// @return                     result of addition
    ///////////////////////////////////////////////////////////////////////////////////////////     
    Vector4 operator+(Vector4);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      same as *= operator
    /// @param scalar               scalar value to be multiplied to vector
    /// @return                     result of multiplication with a scalar
    ///////////////////////////////////////////////////////////////////////////////////////////     
    Vector4 operator*(double scalar);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      this = this / scalar
    /// @param scalar       scalar value by which components are divided
    /// @return                     result of division by scalar
    ///////////////////////////////////////////////////////////////////////////////////////////     
    Vector4 operator/(double scalar);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      scalar / dot product
    /// @param v            vector to be multiplied with    
    /// @return                     double precision dot product result
    ///////////////////////////////////////////////////////////////////////////////////////////     
    double operator*(Vector4 v);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      cross product
    /// @param v            vector to be multiplied with    
    /// @return                     cross product result, right handed
    ///////////////////////////////////////////////////////////////////////////////////////////     
    Vector4 operator&(Vector4);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      normalize the vector
    ///////////////////////////////////////////////////////////////////////////////////////////     
    void normalize(void);
        ///////////////////////////////////////////////////////////////////////////////////////////
        /// @brief                      lenght of the vector
        /// @return                     lenght, double precision
        ///////////////////////////////////////////////////////////////////////////////////////////     
    double getLength(void);
};

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief 3 dimensional vector class with all neccessary methods and properties
///////////////////////////////////////////////////////////////////////////////////////////
class Vector3 {
public:
    double x;
    double y;
    double z;
    inline Vector3();
    inline Vector3(const Vector3 &v);
    inline Vector3(double xIn, double yIn, double zIn);

    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  this = rhs
    /// @param rhs              source value
    /// @return                 rhs
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 &operator=(const Vector3 &rhs); 
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  this = this + vector
    /// @param rhs              vector to be added
    /// @return                 result of addition
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 &operator+=(const Vector3 &rhs); 
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  this.x = this-x + scalar
    /// @param rhs              scalar to be added to each component
    /// @return                 result of addition
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 &operator+=(const double &rhs); 
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  this = this * scalar
    /// @param scalar   scalar value to be multiplied to vector
    /// @return                 result of multiplication with a scalar
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 &operator*=(const double &scalar); 
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  negate vector
    /// @return                 negative vector (all components *(-1))
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 operator-();
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  subtraction
    /// @return                 result of subtraction
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 operator-(Vector3);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  addition
    /// @return                 result of addition
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 operator+(Vector3);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  addition of scalar to each component
    /// @return                 result of addition
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 operator+(double);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  same as *= operator
    /// @param scalar           scalar value to be multiplied to vector
    /// @return                 result of multiplication with a scalar
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 operator*(double scalar);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  this = this / scalar
    /// @param scalar           scalar value by which components are divided
    /// @return                 result of division by scalar
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 operator/(double scalar);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  scalar / dot product
    /// @param v                vector to be multiplied with    
    /// @return                 double precision dot product result
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline double operator*(Vector3 v);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  cross product
    /// @param v                vector to be multiplied with    
    /// @return                 cross product result, right handed
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 operator&(Vector3);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  normalize the vector
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline Vector3 normalize(void);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  lenght of the vector
    /// @return                 lenght, double precision
    /////////////////////////////////////////////////////////////////////////////////////////// 
    inline double getLength(void);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                  lenght of the vector
    /// @return                 lenght, double precision
    /////////////////////////////////////////////////////////////////////////////////////////// 
    static Vector3 getRoots(Vector3);
};

///////////////////////////////////////////////////////////////////////////////////////////
/// @brief 4 by 4 matrix class with all methods usually implemented by direct x D3DXMATRIX 
///  extended library, but not platform or API specific
///////////////////////////////////////////////////////////////////////////////////////////
class Matrix4x4 {
public:
    double _11;
    double _12;
    double _13;
    double _14;

    double _21;
    double _22;
    double _23;
    double _24;

    double _31;
    double _32;
    double _33;
    double _34;

    double _41;
    double _42;
    double _43;
    double _44;

    Matrix4x4();
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      transposes the matrix
    ///////////////////////////////////////////////////////////////////////////////////////////     
    void transpose(void);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      sets matrix to 1 on diagonal 0 everywhere else
    ///////////////////////////////////////////////////////////////////////////////////////////     
    void identity(void);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      Multiplies Matrix times Vector4
    ///////////////////////////////////////////////////////////////////////////////////////////     
    Vector4 operator*(Vector4);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      multiplies matrices, made to mimic direct x's MatrixMultiply 
    /// @param pout         result of the multiplication
    /// @param pm1          first matrix to be multiplied
    /// @param pm2          second matrix to be multiplied
    /// @return                     same as in pout
    ///////////////////////////////////////////////////////////////////////////////////////////     
    static Matrix4x4* MatrixMultiply(Matrix4x4* pout, const Matrix4x4* pm1, const Matrix4x4* pm2);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      creates a matrix which applies the given translation when used on a vector
    /// @param mOut         resulting matrix
    /// @param x            x coordinate
    /// @param y            y coordinate
    /// @param z            z coordinate
    /// @return                     translation matrix
    ///////////////////////////////////////////////////////////////////////////////////////////     
    static Matrix4x4* MatrixTranslation(Matrix4x4* mOut, double x, double y, double z);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      creates rotate matrix for rotation around x axis, made to mimic direct x's MatrixRotationX
    /// @param mOut         result of the rotation
    /// @param angle        angle in radians
    /// @return                     result of the rotation
    ///////////////////////////////////////////////////////////////////////////////////////////     
    static Matrix4x4* MatrixRotationX(Matrix4x4* mOut, double angle);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      creates rotate matrix for rotation around y axis, made to mimic direct x's MatrixRotationY
    /// @param mOut         result of the rotation
    /// @param angle        angle in radians
    /// @return                     result of the rotation
    ///////////////////////////////////////////////////////////////////////////////////////////     
    static Matrix4x4* MatrixRotationY(Matrix4x4* mOut, double angle);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      creates rotate matrix for rotation around z axis, made to mimic direct x's MatrixRotationZ
    /// @param mOut         result of the rotation
    /// @param angle        angle in radians
    /// @return                     result of the rotation
    ///////////////////////////////////////////////////////////////////////////////////////////     
    static Matrix4x4* MatrixRotationZ(Matrix4x4* mOut, double angle);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      creates rotate matrix for rotation using euler angles, made to mimic direct x version
    ///                 http://msdn.microsoft.com/en-us//library/windows/desktop/bb205361 for details
    /// @param pOut         resulting matrix
    /// @param y            yaw
    /// @param p            pitch
    /// @param r            roll
    /// @return                     result of the rotation
    ///////////////////////////////////////////////////////////////////////////////////////////     
    static Matrix4x4* MatrixRotationYawPitchRoll(Matrix4x4* pOut, double y, double p, double r);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      access members
    /// @param index        desired element (0-15 are valid indices)
    /// @return                     (const) rhs version of the element
    ///////////////////////////////////////////////////////////////////////////////////////////     
    const double& operator[](int const index) const {return *((double*)this+index);}
    ///////////////////////////////////////////////////////////////////////////////////////////
    /// @brief                      access members
    /// @param index        desired element (0-15 are valid indices)
    /// @return                     lhs version of the element
    ///////////////////////////////////////////////////////////////////////////////////////////     
    double& operator[](int const index) {return *((double*)this+index);}
};

///////////////definitions of inlined functions
inline Vector3::Vector3( )
{
}

inline Vector3::Vector3(const Vector3 &v)
{
        x=v.x; y=v.y; z=v.z;
}

inline Vector3::Vector3(double xIn, double yIn, double zIn)
{
        x=xIn; y=yIn; z=zIn;
}

inline Vector3 &Vector3::operator=(const Vector3 &v)
{
        x = v.x;
        y = v.y; 
        z = v.z;
        return *this;
}

inline Vector3 &Vector3::operator+=(const Vector3 &v)
{
        Vector3 ret;
        x = x+v.x;
        y = y+v.y; 
        z = z+v.z;
        return *this;
}

inline Vector3 &Vector3::operator+=(const double &s)
{
        Vector3 ret;
        x = x+s;
        y = y+s;
        z = z+s;
        return *this;
}

inline Vector3 &Vector3::operator*=(const double &s)
{
        Vector3 ret;
        x = x*s;
        y = y*s; 
        z = z*s;
        return *this;
}

inline Vector3 Vector3::operator-( )
{
        Vector3 ret;
        ret.x = -x;
        ret.y = -y;
        ret.z = -z;
        return ret;
}

inline Vector3 Vector3::operator-(Vector3 w)
{
        Vector3 ret;
        ret.x = x - w.x;
        ret.y = y - w.y;
        ret.z = z - w.z;
        return ret;
}

inline Vector3 Vector3::operator+(Vector3 w)
{
        Vector3 ret;
        ret.x = x + w.x;
        ret.y = y + w.y;
        ret.z = z + w.z;
        return ret;
}

inline Vector3 Vector3::operator+(double s)
{
        Vector3 ret;
        ret.x = x + s;
        ret.y = y + s;
        ret.z = z + s;
        return ret;
}

inline Vector3 Vector3::operator*(double s)
{
        Vector3 ret(0.0,0.0,0.0);
        ret.x = x*s;
        ret.y = y*s;
        ret.z = z*s;
        return ret;
}

inline Vector3 Vector3::operator/(double s)
{
        Vector3 ret(0.0,0.0,0.0);
        ret.x = x/s;
        ret.y = y/s;
        ret.z = z/s;
        return ret;
}

inline double Vector3::operator*(Vector3 v)
{
        double ret = x*v.x+y*v.y+z*v.z;
        return ret;
}

inline Vector3 Vector3::operator&(Vector3 v)
{
        Vector3 ret(0.0,0.0,0.0);
        ret.x = y*v.z - z*v.y;
        ret.y = z*v.x - x*v.z;
        ret.z = x*v.y - y*v.x;
        return ret;
}

inline Vector3 Vector3::normalize(void) 
{
    double vlen = sqrt(x*x+y*y+z*z);
    x/=vlen; y/=vlen; z/=vlen;
    return *this;
}

inline double Vector3::getLength(void) {
    return sqrt(x*x+y*y+z*z);
}

#endif
