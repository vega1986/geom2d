#pragma once
#include "baseCurve.h"
#include "point.h"
#include "avector.h"
#include "theLibGeom2d.h"

namespace geom2d
{
  // представляет отрезок на плоскости
  // vec{r}(t) = vec{m_startPoint} + t * vec{m_finishPoint - m_startPoint}
  class segmentCurve : public baseCurve
  {
  public:
    THELIBGEOM2D_API segmentCurve(const point first, const point second);

    THELIBGEOM2D_API virtual point getPoint(double t) const;

    THELIBGEOM2D_API virtual point getVelocity(double t) const;

    THELIBGEOM2D_API virtual double parameterMin() const;

    THELIBGEOM2D_API virtual double parameterMax() const;

    // возвращает параметр ближайшей к данной точке точки на отрезке
    THELIBGEOM2D_API double nearestTo(const point p) const;

    // возвращает направляющий вектор (единичный)
    THELIBGEOM2D_API vector getTau() const;

    // возвращает расстояние от отрезка до точки
    THELIBGEOM2D_API double distanceTo(const point pnt) const;

  private:
    point m_startPoint;
    point m_finishPoint;
  };

}