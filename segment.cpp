#include <stdexcept>
#include "segment.h"
#include "vector.h"

double geom2d::distance(const point p, const segment s)
{
  vector ap{ s.a, p };
  vector ab{ s.a, s.b };

  const auto pIsAfterA = (ap, ab) >= 0.0;

  vector bp{ s.b, p };
  vector ba{ s.b, s.a };

  const auto pIsAfterB = (bp, ba) >= 0.0;

  if (pIsAfterA and pIsAfterB)
  {
    vector abDir = ab;
    abDir.normalize();
    const double dist = std::abs(ap ^ abDir);
    return dist;
  }
  else if (pIsAfterA)
  {
    const double dist = bp.length();
    return dist;
  }
  else if (pIsAfterB)
  {
    const double dist = ap.length();
    return dist;
  }
  else
  {
    throw std::logic_error("geom2d::distance impossible");
  }
}
