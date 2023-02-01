#pragma once
#include "common.h"
#include "point.h"

#include "libgeom2d/theLibGeom2d/theLibGeom2d.h"

namespace geom2d
{

  class baseCurve
  {

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