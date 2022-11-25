#include "baseCurve.h"
#include "findFunctionRoots.h"

#include <algorithm>
#include <stdexcept>

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
