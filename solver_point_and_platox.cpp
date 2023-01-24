#include "solver_point_and_platox.h"
#include "StatOfCurvePiece.h"
#include "findFunctionRoots.h"

//---------------------------------------------------------------------------------------------------------------------

std::optional<double> geom2d::point_and_platox::solver::execute() const
{
  StatOfCurvePiece scp{ tmin, curve.getPoint(tmin), tmax, curve.getPoint(tmax) };

  point pntLeft = scp.pointOfxmin();
  point pntRight = scp.pointOfxmax();

  if (P.x >= pntRight.x)
  {
    if (point::isSame(P, pntRight))
    {
      return scp.tOfxmax();
    }
    else
    {
      return std::nullopt;
    }
  }
  else if (P.x <= pntLeft.x)
  {
    if (point::isSame(P, pntLeft))
    {
      return scp.tOfxmin();
    }
    else
    {
      return std::nullopt;
    }
  }
  else
  {
    const auto the_t = curve.tofX(tmin, tmax, P.x);
    const point the_point = curve.getPoint(the_t);
    if (point::isSame(P, the_point))
    {
      return the_t;
    }
    else
    {
      return std::nullopt;
    }
  }
}

//---------------------------------------------------------------------------------------------------------------------
