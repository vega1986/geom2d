#include "solver_point_and_platoy.h"
#include "StatOfCurvePiece.h"
#include "findFunctionRoots.h"

//---------------------------------------------------------------------------------------------------------------------

std::optional<double> geom2d::point_and_platoy::solver::execute() const
{
  StatOfCurvePiece scp{ tmin, curve.getPoint(tmin), tmax, curve.getPoint(tmax) };

  point pntUp = scp.pointOfymax();
  point pntDown = scp.pointOfymin();

  if (P.y >= pntUp.y)
  {
    if (point::isSame(P, pntUp))
    {
      return scp.tOfymax();
    }
    else
    {
      return std::nullopt;
    }
  }
  else if (P.y <= pntDown.y)
  {
    if (point::isSame(P, pntDown))
    {
      return scp.tOfymin();
    }
    else
    {
      return std::nullopt;
    }
  }
  else
  {
    const auto the_t = curve.tofY(tmin, tmax, P.y);
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
