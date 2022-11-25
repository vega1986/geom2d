#pragma once

#include "theLibGeom2d.h"

#include <cmath>

namespace geom2d
{

  struct point
  {
    double x{ 0.0 };
    double y{ 0.0 };
    THELIBGEOM2D_API static double distance(const point& p, const point& q);
    THELIBGEOM2D_API static bool isSame(const point& p, const point& q);
    THELIBGEOM2D_API friend point operator+ (const point& p, const point& q);
    THELIBGEOM2D_API friend point operator- (const point& p, const point& q);
    THELIBGEOM2D_API friend point operator* (const point& p, const double value);
    THELIBGEOM2D_API friend point operator* (const double value, const point& p);
    THELIBGEOM2D_API friend point operator/ (const point& p, const double value);
  };

}