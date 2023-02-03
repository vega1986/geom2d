#include "baseCurve.h"
#include "findFunctionRoots.h"

#include <algorithm>
#include <stdexcept>

//---------------------------------------------------------------------------------------------------------------------

geom2d::curveClass
  geom2d::baseCurve::getCurveClass
  (
    const double tmin,
    const double tmax,
    const baseCurve& curve
  )
{
  const auto PointBegin = curve.getPoint(tmin);
  const auto PointEnd = curve.getPoint(tmax);

  const auto [xbeg, ybeg] = PointBegin;
  const auto [xend, yend] = PointEnd;
  constexpr auto tol2d = math::tolerance::tolPoint;

  const bool xInc = (xend - xbeg) > tol2d;
  const bool xDec = (xbeg - xend) > tol2d;

  const bool yInc = (yend - ybeg) > tol2d;
  const bool yDec = (ybeg - yend) > tol2d;

  if (xInc)
  {
    if (yDec)
    {
      return curveClass::Screen;
    }
    else if (yInc)
    {
      return curveClass::Normal;
    }
    else
    {
      return curveClass::PlatoY;
    }
  }
  else if (xDec)
  {
    if (yDec)
    {
      return curveClass::Normal;
    }
    else if (yInc)
    {
      return curveClass::Screen;
    }
    else
    {
      return curveClass::PlatoY;
    }
  }
  else
  {
    if (yDec)
    {
      return curveClass::PlatoX;
    }
    else if (yInc)
    {
      return curveClass::PlatoX;
    }
    else
    {
      return curveClass::Point;
    }
  }
}

//---------------------------------------------------------------------------------------------------------------------

const double geom2d::baseCurve::tofX(const double tmin, const double tmax, const double theX) const
{
  const auto P = getPoint(tmin);
  const auto Q = getPoint(tmax);

  const auto xmin = std::min(std::initializer_list{ P.x, Q.x });
  const auto xmax = std::max(std::initializer_list{ P.x, Q.x });

  if ((theX < xmin) or (theX > xmax)) throw std::logic_error("bad theX argument in geom2d::baseCurve::tofX");

  auto func = [this, xo = theX](const double t) -> const double
  {
    return getPoint(t).x - xo;
  };

  const auto theT = math::findUniqueFunctionRoot(tmin, tmax, func);
  return theT;
}

//---------------------------------------------------------------------------------------------------------------------

const double geom2d::baseCurve::tofY(const double tmin, const double tmax, const double theY) const
{
  const auto P = getPoint(tmin);
  const auto Q = getPoint(tmax);

  const auto ymin = std::min(std::initializer_list{ P.y, Q.y });
  const auto ymax = std::max(std::initializer_list{ P.y, Q.y });

  if ((theY < ymin) or (theY > ymax)) throw std::logic_error("bad theX argument in geom2d::baseCurve::tofY");

  auto func = [this, yo = theY](const double t) -> const double
  {
    return getPoint(t).y - yo;
  };

  const auto theT = math::findUniqueFunctionRoot(tmin, tmax, func);
  return theT;
}

//---------------------------------------------------------------------------------------------------------------------
