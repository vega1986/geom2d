#pragma once
#include "common.h"
#include "point.h"

#include "theLibGeom2d.h"

namespace geom2d
{
  
  // �������������� ��������� ������� ������ �� ������� ������������
  enum class curveClass : unsigned short int
  {
    Screen = 1, // (dx / dt > 0 && dy / dt < 0) || (dx / dt < 0 && dy / dt > 0)
    Normal = 2, // (dx / dt > 0 && dy / dt > 0) || (dx / dt < 0 && dy / dt < 0)
    PlatoX = 3, // (dx / dt == 0 && (dy / dt>0 || dy / dt < 0)
    PlatoY = 4, // (dy / dt == 0 && (dx / dt>0 || dx / dt < 0)
    Point = 5  // (dx / dt == 0 && dy / dt == 0)
  };

  class baseCurve
  {
  public:
    // �������� ����� ������ �� ���������� �������!
    static
      curveClass
      getCurveClass(const double tmin, const double tmax, const baseCurve& curve);

  protected:

    THELIBGEOM2D_API
    baseCurve () = default;

  public:

    // �������� ����� �� ������ ��� ��������� ��������� t
    THELIBGEOM2D_API
    virtual point getPoint (double t) const = 0;

    // �������� ����������� �� t �� ����� �� ������ ��� ��������� ��������� t
    THELIBGEOM2D_API
    virtual point getVelocity (double t) const = 0;

    // ����������� �������� ��������� t
    THELIBGEOM2D_API
    virtual double parameterMin () const = 0;

    // ������������ �������� ��������� t
    THELIBGEOM2D_API
    virtual double parameterMax () const = 0;

    // �������� �������� t ��� ������� x ��������� �������� theX �� ������� ������������
    THELIBGEOM2D_API
    const double tofX(const double tmin, const double tmax, const double theX) const;

    // �������� �������� t ��� ������� y ��������� �������� theY �� ������� ������������
    THELIBGEOM2D_API
    const double tofY(const double tmin, const double tmax, const double theY) const;

  };

}