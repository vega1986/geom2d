#pragma once
#include "baseCurve.h"
#include "point.h"
#include "avector.h"
#include "theLibGeom2d.h"

namespace geom2d
{
  // ������������ ������� �� ���������
  // vec{r}(t) = vec{m_startPoint} + t * vec{m_finishPoint - m_startPoint}
  class segmentCurve : public baseCurve
  {
  public:
    THELIBGEOM2D_API segmentCurve(const point first, const point second);

    THELIBGEOM2D_API virtual point getPoint(double t) const;

    THELIBGEOM2D_API virtual point getVelocity(double t) const;

    THELIBGEOM2D_API virtual double parameterMin() const;

    THELIBGEOM2D_API virtual double parameterMax() const;

    // ���������� �������� ��������� � ������ ����� ����� �� �������
    THELIBGEOM2D_API double nearestTo(const point p) const;

    // ���������� ������������ ������ (���������)
    THELIBGEOM2D_API vector getTau() const;

    // ���������� ���������� �� ������� �� �����
    THELIBGEOM2D_API double distanceTo(const point pnt) const;

  private:
    point m_startPoint;
    point m_finishPoint;
  };

}