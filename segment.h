#pragma once

#include "point.h"

namespace geom2d
{
  struct segment
  {
    friend double distance(const point p, const segment s);

    point a{0.0, 0.0};
    point b{0.0, 0.0};
  };

}


