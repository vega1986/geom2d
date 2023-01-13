#pragma once
#include "common.h"
#include "point.h"

namespace geom2d
{

  class baseCurve
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
      
  };

}