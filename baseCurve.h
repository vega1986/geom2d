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

    // получить точку на кривой для заданного параметра t
    THELIBGEOM2D_API
    virtual point getPoint (double t) const = 0;

    // получить производные по t от точки на кривой для заданного параметра t
    THELIBGEOM2D_API
    virtual point getVelocity (double t) const = 0;

    // минимальное значение параметра t
    THELIBGEOM2D_API
    virtual double parameterMin () const = 0;

    // максимальное значение параметра t
    THELIBGEOM2D_API
    virtual double parameterMax () const = 0;

    // получаем значение t при котором x принимает значение theX на участке монотонности
    THELIBGEOM2D_API
    const double tofX(const double tmin, const double tmax, const double theX) const;

    // получаем значение t при котором y принимает значение theY на участке монотонности
    THELIBGEOM2D_API
    const double tofY(const double tmin, const double tmax, const double theY) const;

  };

}