#pragma once

#include "libgeom2d/theLibGeom2d/theLibGeom2d.h"

#include <cmath>

namespace geom2d
{

  struct THELIBGEOM2D_API point
  {
    double x{ 0.0 };
    double y{ 0.0 };
    static double distance(const point& p, const point& q);
    static bool isSame(const point& p, const point& q);
    friend point operator+ (const point& p, const point& q);
    friend point operator- (const point& p, const point& q);
    friend point operator* (const point& p, const double value);
    friend point operator* (const double value, const point& p);
    friend point operator/ (const point& p, const double value);
  };

}