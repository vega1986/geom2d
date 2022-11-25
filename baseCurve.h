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

    // получить точку на кривой для заданного параметра t
    virtual point getPoint (double t) const = 0;

    // получить производные по t от точки на кривой для заданного параметра t
    virtual point getVelocity (double t) const = 0;

    // минимальное значение параметра t
    virtual double parameterMin () const = 0;

    // максимальное значение параметра t
    virtual double parameterMax () const = 0;
      
  };

}