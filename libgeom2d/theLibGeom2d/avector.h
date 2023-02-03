#pragma once

#include "point.h"
#include "theLibGeom2d.h"

namespace geom2d
{

  struct THELIBGEOM2D_API vector {

    vector() = default;
    vector(const double xx, const double yy);
    vector(const point p, const point q);
    vector(const point p);

    double length() const;
    
    void normalize();
    const vector normalized() const;

    void reverse();
    const vector reversed() const;

    THELIBGEOM2D_API friend vector operator+ (const vector p, const vector q);
    THELIBGEOM2D_API friend vector operator- (const vector p, const vector q);
    THELIBGEOM2D_API friend vector operator* (const vector p, const double value);
    THELIBGEOM2D_API friend vector operator* (const double value, const vector p);
    THELIBGEOM2D_API friend vector operator/ (const vector p, const double value);
    // скалярное произведение
    THELIBGEOM2D_API friend double operator, (const vector p, const vector q);
    // z-компонента векторного произведения
    THELIBGEOM2D_API friend double operator^ (const vector p, const vector q);
    
    double x{0.0};
    double y{0.0};
  };

}