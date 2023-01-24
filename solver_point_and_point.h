#pragma once
#include "point.h"

namespace geom2d
{
  namespace point_and_point
  {
    class solver
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
