#pragma once
#include "common.h"
#include "point.h"

#include "theLibGeom2d.h"

namespace geom2d
{
  
  // классифицируем поведение участка кривой на участке монотонности
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
    // Получаем класс кривой на монотонном участке!
    static
      curveClass
      getCurveClass(const double tmin, const double tmax, const baseCurve& curve);

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