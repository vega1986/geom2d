#pragma once
#include "point.h"
#include "baseCurve.h"

#include <optional>

namespace geom2d
{
  namespace point_and_platox
  {
    class solver
    {
    public:
      solver(const point thePoint, const double theTmin, const double theTmax, const baseCurve& theCurve)
        : P(thePoint), tmin(theTmin), tmax(theTmax), curve(theCurve) {}

      std::optional<double> execute() const;

    private:
      const point P;
      const double tmin;
      const double tmax;
      const baseCurve& curve;
    };
  }
}