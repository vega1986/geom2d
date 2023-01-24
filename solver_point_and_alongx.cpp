#include "solver_point_and_alongx.h"
#include "StatOfCurvePiece.h"

//---------------------------------------------------------------------------------------------------------------------

std::optional<double> geom2d::point_and_curve_alongx::solver::execute() const
{
  StatOfCurvePiece scp{ tmin, curve.getPoint(tmin), tmax, curve.getPoint(tmax) };
  if (P.x <= scp.pointOfxmin().x)
  {
    if (point::isSame(P, scp.pointOfxmin()))
    {
      return scp.tOfxmin();
    }
    else
    {
      return std::nullopt;
    }
  }
  else if (P.x >= scp.pointOfxmax().x)
  {
    if (point::isSame(P, scp.pointOfxmax()))
    {
      return scp.tOfxmax();
    }
    else
    {
      return std::nullopt;
    }
  }
  else
  {
    return curve.tofX(tmin, tmax, P.x);
  }
  return std::nullopt;
}

//---------------------------------------------------------------------------------------------------------------------
