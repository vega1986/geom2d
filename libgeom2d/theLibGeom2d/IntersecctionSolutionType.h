#pragma once
#include <tuple>
#include "point.h"

namespace geom2d
{
  using IntersecctionSolutionType = std::tuple<geom2d::point, double, double>;
}