#include "solver_point_and_alongy.h"
#include "StatOfCurvePiece.h"

//---------------------------------------------------------------------------------------------------------------------

std::optional<double> geom2d::point_and_curve_alongy::solver::execute() const
{
  StatOfCurvePiece scp{ tmin, curve.getPoint(tmin), tmax, curve.getPoint(tmax) };
  if (P.y <= scp.pointOfymin().y)
  {
    if (point::isSame(P, scp.pointOfymin()))
    {
      return scp.tOfymin();
    }
    else
    {
      return std::nullopt;
    }
  }
  else if (P.y >= scp.pointOfymax().y)
  {
    if (point::isSame(P, scp.pointOfymax()))
    {
      return scp.tOfymax();
    }
    else
    {
      return std::nullopt;
    }
  }
  else
  {
    return curve.tofY(tmin, tmax, P.y);
  }
  return std::nullopt;
}

//---------------------------------------------------------------------------------------------------------------------
