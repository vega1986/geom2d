#pragma once

#include "point.h"

namespace geom2d
{

  struct vector {

    vector() = default;
    vector(const double xx, const double yy);
    vector(const point p, const point q);
    vector(const point p);

    double length() const;
    void normalize();

    friend vector operator+ (const vector p, const vector q);
    friend vector operator- (const vector p, const vector q);
    friend vector operator* (const vector p, const double value);
    friend vector operator* (const double value, const vector p);
    friend vector operator/ (const vector p, const double value);
    // скалярное произведение
    friend double operator, (const vector p, const vector q);
    // z-компонента векторного произведения
    friend double operator^ (const vector p, const vector q);
    
    double x{0.0};
    double y{0.0};
  };

}