#pragma once
#include <tuple>
#include "point.h"

namespace geom2d
{
  // Тип решения для алгоритма поиска точки пересечения: сама точка пересечения и значения параметров кривых,
  // при которых пересечение наблюдается.
  using IntersecctionSolutionType = std::tuple<geom2d::point, double, double>;
}