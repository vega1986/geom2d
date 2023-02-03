#pragma once
#include "baseCurve.h"

namespace geom2d
{
  namespace extrema
  {
    namespace point_and_curve
    {

      class nearestPoint
      {
      public:
        nearestPoint(const point theP, const double theTmin, const double theTmax, const baseCurve& theCurve)
          : P{ theP }, tmin(theTmin), tmax{ theTmax }, curve{ theCurve } {}

        std::pair<double, double>
          execute() const;

      private:
        const point P;
        const double tmin;
        const double tmax;
        const baseCurve& curve;
      };
    }
  }
}