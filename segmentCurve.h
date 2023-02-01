#pragma once
#include "baseCurve.h"
#include "point.h"
#include "vector.h"
#include "libgeom2d/theLibGeom2d/theLibGeom2d.h"

namespace geom2d
{
  // ������������ ������� �� ���������
  // vec{r}(t) = vec{m_startPoint} + t * vec{m_finishPoint - m_startPoint}
  class THELIBGEOM2D_API segmentCurve : public baseCurve
  {
  public:
    segmentCurve(const point first, const point second);

    virtual point getPoint(double t) const;

    virtual point getVelocity(double t) const;

    virtual double parameterMin() const;

    virtual double parameterMax() const;

    // ���������� �������� ��������� � ������ ����� ����� �� �������
    double nearestTo(const point p) const;

    // ���������� ������������ ������ (���������)
    vector getTau() const;

    // ���������� ���������� �� ������� �� �����
    double distanceTo(const point pnt) const;

  private:
    point m_startPoint;
    point m_finishPoint;
  };

}