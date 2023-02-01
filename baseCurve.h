#pragma once
#include "common.h"
#include "point.h"

#include "libgeom2d/theLibGeom2d/theLibGeom2d.h"

namespace geom2d
{

  class THELIBGEOM2D_API baseCurve
  {

  protected:

    baseCurve () = default;

  public:

    // �������� ����� �� ������ ��� ��������� ��������� t
    virtual point getPoint (double t) const = 0;

    // �������� ����������� �� t �� ����� �� ������ ��� ��������� ��������� t
    virtual point getVelocity (double t) const = 0;

    // ����������� �������� ��������� t
    virtual double parameterMin () const = 0;

    // ������������ �������� ��������� t
    virtual double parameterMax () const = 0;

    // �������� �������� t ��� ������� x ��������� �������� theX �� ������� ������������
    const double tofX(const double tmin, const double tmax, const double theX) const;

    // �������� �������� t ��� ������� y ��������� �������� theY �� ������� ������������
    const double tofY(const double tmin, const double tmax, const double theY) const;

  };

}