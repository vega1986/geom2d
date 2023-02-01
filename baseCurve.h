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

    // получить точку на кривой для заданного параметра t
    virtual point getPoint (double t) const = 0;

    // получить производные по t от точки на кривой для заданного параметра t
    virtual point getVelocity (double t) const = 0;

    // минимальное значение параметра t
    virtual double parameterMin () const = 0;

    // максимальное значение параметра t
    virtual double parameterMax () const = 0;

    // получаем значение t при котором x принимает значение theX на участке монотонности
    const double tofX(const double tmin, const double tmax, const double theX) const;

    // получаем значение t при котором y принимает значение theY на участке монотонности
    const double tofY(const double tmin, const double tmax, const double theY) const;

  };

}