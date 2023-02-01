#pragma once
#include "baseCurve.h"
#include "point.h"
#include "libgeom2d/theLibGeom2d/theLibGeom2d.h"

namespace geom2d
{

  class pointCurve : public baseCurve
  {
  public:
    THELIBGEOM2D_API pointCurve(const point thePoint);

    THELIBGEOM2D_API virtual point getPoint(double) const;

    THELIBGEOM2D_API virtual point getVelocity(double) const;

    THELIBGEOM2D_API virtual double parameterMin() const;

    THELIBGEOM2D_API virtual double parameterMax() const;
    
  private:
    point m_point;
  };

}
