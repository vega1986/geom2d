#pragma once
#include "baseCurve.h"

namespace geom2d
{
  namespace extrema
  {
    namespace point_and_curve
    {
      //  ласс дл€ расчЄта рассто€ни€ между точкой и произвольной кривой.
      // Ѕлижайша€ точка можетбыть либо на концах кривой либо внутри ќƒ«,
      // но тогда отрезок, соедин€ющий точку и точку на привой должен быть перпендикул€рен кривой
      // (касательной к кривой).
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