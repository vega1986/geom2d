#pragma once
#include "baseCurve.h"
#include "point.h"
#include "libgeom2d/theLibGeom2d/theLibGeom2d.h"

namespace geom2d
{

  class THELIBGEOM2D_API pointCurve : public baseCurve
  {
  public:
    pointCurve(const point thePoint);

    virtual point getPoint(double) const;

    virtual point getVelocity(double) const;

    virtual double parameterMin() const;

    virtual double parameterMax() const;

  private:
    point m_point;
  };

}
