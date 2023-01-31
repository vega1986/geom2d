#pragma once
#include "baseCurve.h"

#include "point.h"

namespace geom2d
{

  class pointCurve : public baseCurve
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
