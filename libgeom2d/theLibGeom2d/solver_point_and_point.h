#pragma once
#include "point.h"
#include "theLibGeom2d.h"

namespace geom2d
{
  namespace point_and_point
  {
    //  ласс дл€ поиска точки пересечени€ кривых Point.
    // ¬озвращает true если точки пересекаюстс€.
    class THELIBGEOM2D_API solver
    {
    public:
      solver(const point theP, const point theQ):p(theP), q(theQ) {}

      const bool execute() const
      {
        return point::isSame(p, q);
      }

    private:

      const point p;
      const point q;

    };
  }
}
