#pragma once
#include "point.h"
#include "baseCurve.h"
#include "CommonRangeHelper.h"

#include <optional>
#include <vector>


namespace geom2d
{
  namespace points_and_curve_alongaxis
  {
    // Класс для поиска точек пересечения точек с кривой, не являющейся точкой на заданном участке монотонности.
    template <CurveDataGetter dataGetter>
    class solver
    {
    public:
      template<class anyPointsContainer>
      solver(const anyPointsContainer& points, const double theTmin, const double theTmax, const baseCurve& theCurve)
        :Ps(begin(points), end(points)), tmin{theTmin}, tmax{theTmax}, curve{theCurve} {}

      std::optional<std::pair<size_t, double>>
        execute() const
      {
        dataGetter theGetter{ tmin, tmax, curve };
        //for (const auto P : Ps)
        for (size_t j = 0; j < Ps.size(); ++j)
        {
          const auto P = Ps[j];
          if (dataGetter::abscissaOf(P) <= theGetter.getMinCoord())
          {
            if (point::isSame(P, theGetter.getPointOfMin()))
            {
              return std::pair{ j, theGetter.getTofMin() };
            }
            else
            {
              continue;
            }
          }
          else if (dataGetter::abscissaOf(P) >= theGetter.getMaxCoord())
          {
            if (point::isSame(P, theGetter.getPointOfMax()))
            {
              return std::pair{ j, theGetter.getTofMax() };
            }
            else
            {
              continue;
            }
          }
          else
          {
            const auto tofCurve = theGetter.getTofCoord(dataGetter::abscissaOf(P));
            const auto pofCurve = curve.getPoint(tofCurve);
            if (point::isSame(P, pofCurve))
            {
              return std::pair{ j, tofCurve };
            }
            else
            {
              continue;
            }
          }
        }
        return std::nullopt;
      }

    private:
      const std::vector<point> Ps;
      const double tmin;
      const double tmax;
      const baseCurve& curve;
    };
  }
}