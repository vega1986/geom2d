#pragma once
#include "baseCurve.h"
#include "point.h"
#include "avector.h"
#include "theLibGeom2d.h"

#include <optional>

namespace geom2d
{
  // представляет отрезок на плоскости
  // vec{r}(t) = vec{m_startPoint} + t * vec{m_finishPoint - m_startPoint}
  class segmentCurve : public baseCurve
  {
  public:
    THELIBGEOM2D_API segmentCurve(const point first, const point second);

    THELIBGEOM2D_API
    virtual ~segmentCurve() = default;

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

    THELIBGEOM2D_API const point start() const;

    THELIBGEOM2D_API const point finish() const;

    THELIBGEOM2D_API static
      std::optional<std::pair<double, double>>
      intersectionOf(segmentCurve seg1, segmentCurve seg2);

    THELIBGEOM2D_API static
      std::tuple<double, double, double>
      distanceBetween(segmentCurve seg1, segmentCurve seg2);

    THELIBGEOM2D_API const double length() const;

  private:
    point m_startPoint;
    point m_finishPoint;
  };

}