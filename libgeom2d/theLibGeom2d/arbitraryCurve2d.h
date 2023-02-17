#pragma once
#include "common.h"
#include "point.h"
#include "baseCurve.h"

#include "theLibGeom2d.h"

namespace geom2d
{
  /////////////
  //         //
  // CONCEPT //
  //         //
  /////////////
  template <class T>
  concept OneArgFunction = requires(T a, const double t)
  {
    { a(t) } -> std::convertible_to<double>;
  };
  

  template<OneArgFunction XPointf, OneArgFunction YPointf, OneArgFunction XVelocityf, OneArgFunction YVelocityf>
  class arbitraryCurve final : public baseCurve
  {
  private:
    XPointf px;
    YPointf py;

    XVelocityf ux;
    YVelocityf uy;

    const double tmin;
    const double tmax;
  public:

    /*THELIBGEOM2D_API*/
    arbitraryCurve(const XPointf& thePx, const YPointf& thePy, const XVelocityf& theUx, const YVelocityf& theUy, const double theTmin, const double theTmax)
      :px{ thePx }, py{ thePy }, ux{ theUx }, uy{ theUy }, tmin{ theTmin }, tmax{ theTmax } {}

    virtual ~arbitraryCurve() = default;
    
    // получить точку на кривой для заданного параметра t
    /*THELIBGEOM2D_API*/
    virtual point getPoint(double t) const
    {
      return point{ px(t), py(t) };
    }

    // получить производные по t от точки на кривой для заданного параметра t
    /*THELIBGEOM2D_API*/
    virtual point getVelocity(double t) const
    {
      return point{ ux(t), uy(t) };
    }

    // минимальное значение параметра t
    /*THELIBGEOM2D_API*/
    virtual double parameterMin() const
    {
      return tmin;
    }

    // максимальное значение параметра t
    /*THELIBGEOM2D_API*/
    virtual double parameterMax() const
    {
      return tmax;
    }

  };

}