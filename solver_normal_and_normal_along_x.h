#pragma once
#include "baseCurve.h"
#include "IntersecctionSolutionType.h"

#include <vector>
#include <optional>

namespace geom2d
{
  namespace solver_normal_and_normal_along_x
  {
    class solver
    {
    public:
      solver(const double theTofmin1, const double theTofmax1, const baseCurve& theCurve1,
             const double theTofmin2, const double theTofmax2, const baseCurve& theCurve2)
        :
        tofmin1(theTofmin1), tofmax1(theTofmax1), curve1(theCurve1),
        tofmin2(theTofmin2), tofmax2(theTofmax2), curve2(theCurve2){}

      std::vector<IntersecctionSolutionType>
        execute() const;

      static
        std::optional<std::pair<double, double>>
        shiftForDiverge
        (
          const double tofmin1,
          const double tofmax1,
          const baseCurve& curve1,
          const double tofmin2,
          const double tofmax2,
          const baseCurve& curve2
        );

      static
        std::optional<std::tuple<double, double, bool>>
        shiftForConverge
        (
          const double tofminUpper,
          const double tofmaxUpper,
          const baseCurve& curveUpper,
          const double tofminLower,
          const double tofmaxLower,
          const baseCurve& curveLower
        );

    private:
      const double tofmin1;
      const double tofmax1;
      const baseCurve& curve1;

      const double tofmin2;
      const double tofmax2;
      const baseCurve& curve2;
    };
  }
}