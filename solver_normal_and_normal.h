#pragma once
#include "baseCurve.h"
#include "IntersecctionSolutionType.h"
#include "CommonRangeHelper.h"

#include <vector>
#include <optional>
#include <algorithm>

namespace geom2d
{
  namespace solver_normal_and_normal
  {
    // здесь используется паттерн проектирования builder (строитель)
    // так как алгоритмически не важно вдоль оси X|Y мы идём
    template<CurveDataGetter dataGetter>
    class solver
    {
    public:
      solver(const double theTofmin1, const double theTofmax1, const baseCurve& theCurve1,
        const double theTofmin2, const double theTofmax2, const baseCurve& theCurve2)
        :
        tofmin1(theTofmin1), tofmax1(theTofmax1), curve1(theCurve1),
        tofmin2(theTofmin2), tofmax2(theTofmax2), curve2(theCurve2) {}

      std::vector<IntersecctionSolutionType>
        execute() const
      {
        std::vector<geom2d::IntersecctionSolutionType> result;

        auto t1 = tofmin1;
        auto p1 = curve1.getPoint(tofmin1);
        const auto pointOfMax1 = curve1.getPoint(tofmax1);

        auto t2 = tofmin2;
        auto p2 = curve2.getPoint(tofmin2);
        const auto pointOfMax2 = curve2.getPoint(tofmax2);

        bool shiftFacedWithEnd = false;

        while (true)
        {
          if (point::isSame(p1, p2))
          {
            // в этой ветке утилизирован случай, когда ordinate(p1) == ordinate(p2)
            // так как если это так, то dist(p1, p2) <= 0.1 * tolerance(точки)
            // (в соответствии с точностью нахождения корня уравнения, которая в 10 раз меньше толеранса точки)
            result.push_back(IntersecctionSolutionType{ 0.5 * (p1 + p2), t1, t2 });
            if (shiftFacedWithEnd) break;
            const auto theResultOfShift = shiftForDiverge(t1, tofmax1, curve1, t2, tofmax2, curve2);
            if (not theResultOfShift) break;

            const auto [new_t1, new_t2] = theResultOfShift.value();
            // update t1 & t2, p1 & p2
            t1 = new_t1;
            t2 = new_t2;

            p1 = curve1.getPoint(t1);
            p2 = curve2.getPoint(t2);
          }
          else
          {
            if (dataGetter::ordinateOf(p2) > dataGetter::ordinateOf(p1))
            {
              // 1 - lower curve
              // 2 - upper curve
              const auto theResultOfShift = shiftForConverge(t2, tofmax2, curve2, t1, tofmax1, curve1);
              if (not theResultOfShift) break;

              auto [theT2, theT1, isLast] =
                theResultOfShift.value();

              // update t1 & t2, p1 & p2
              t1 = theT1;
              t2 = theT2;

              p1 = curve1.getPoint(t1);
              p2 = curve2.getPoint(t2);

              shiftFacedWithEnd = isLast;
              continue;
            }
            else if (dataGetter::ordinateOf(p1) > dataGetter::ordinateOf(p2))
            {
              // 1 - upper curve
              // 2 - lower curve
              const auto theResultOfShift = shiftForConverge(t1, tofmax1, curve1, t2, tofmax2, curve2);
              if (not theResultOfShift) break;

              auto [theT1, theT2, isLast] =
                theResultOfShift.value();

              // update t1 & t2, p1 & p2
              t1 = theT1;
              t2 = theT2;

              p1 = curve1.getPoint(t1);
              p2 = curve2.getPoint(t2);

              shiftFacedWithEnd = isLast;
              continue;
            }
            else
            {
              throw std::logic_error("impossible situation in function geom2d::solver_normal_and_normal::solver::execute");
            }
          }
        }
        return result;
      }

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
        )
      {
        const auto tmin1 = std::min(tofmin1, tofmax1);
        const auto tmax1 = std::max(tofmin1, tofmax1);
        dataGetter getter1{ tmin1 , tmax1 , curve1 };

        const auto tmin2 = std::min(tofmin2, tofmax2);
        const auto tmax2 = std::max(tofmin2, tofmax2);
        dataGetter getter2{ tmin2 , tmax2 , curve2 };

        const auto pointMore1 = curve1.getPoint(tofmax1);
        const auto pointMore2 = curve2.getPoint(tofmax2);

        auto p1 = curve1.getPoint(tofmin1);
        auto p2 = curve2.getPoint(tofmin2);

        auto t1 = tofmin1;
        auto t2 = tofmin2;

        auto theShift = math::tolerance::tolPoint;
        const auto initial_vline = std::max(dataGetter::abscissaOf(p1), dataGetter::abscissaOf(p2));

        while (true)
        {
          if (not point::isSame(p1, p2))
          {
            return std::pair{ t1, t2 };
          }
          const double vline = initial_vline + theShift;
          theShift *= 2;
          if ((vline > dataGetter::abscissaOf(pointMore1)) or (vline > dataGetter::abscissaOf(pointMore2)))
          {
            return std::nullopt;
          }
          
          // *******************************************
          // update t1
          if (vline < dataGetter::abscissaOf(pointMore1))
          {
            t1 = getter1.getTofAbscissa(vline);
            // instead of curve1.tofX(tmin1, tmax1, vline);
          }
          else
          {
            t1 = tofmax1;
          }
          // update t2
          if (vline < dataGetter::abscissaOf(pointMore2))
          {
            t2 = getter2.getTofAbscissa(vline);
            // instead of curve2.tofX(tmin2, tmax2, vline);
          }
          else
          {
            t2 = tofmax2;
          }
          // update p1 & p2
          p1 = curve1.getPoint(t1);
          p2 = curve2.getPoint(t2);
        }
      }

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
        )
      {
        const auto tminUpper = std::min(tofminUpper, tofmaxUpper);
        const auto tmaxUpper = std::max(tofminUpper, tofmaxUpper);
        dataGetter getterUpper{ tminUpper, tmaxUpper, curveUpper };

        const auto tminLower = std::min(tofminLower, tofmaxLower);
        const auto tmaxLower = std::max(tofminLower, tofmaxLower);
        dataGetter getterLower{tminLower, tmaxLower, curveLower};

        const auto pointMoreUpper = curveUpper.getPoint(tofmaxUpper);
        const auto pointMoreLower = curveLower.getPoint(tofmaxLower);

        auto pu = curveUpper.getPoint(tofminUpper);
        auto pl = curveLower.getPoint(tofminLower);

        auto tu = tofminUpper;
        auto tl = tofminLower;

        while (true)
        {
          if (point::isSame(pu, pl))
          {
            return std::tuple{ tu, tl, false };
          }
          else
          {
            const auto hline = dataGetter::ordinateOf(pu);

            if (hline >= dataGetter::ordinateOf(pointMoreLower))
            {
              const auto vline = dataGetter::abscissaOf(pointMoreLower);
              const auto abscissaOfMoreUpper = dataGetter::abscissaOf(pointMoreUpper);
              if (vline < abscissaOfMoreUpper)
              {
                const auto tonUpper = getterUpper.getTofAbscissa(vline);
                // instead of curveUpper.tofX(tminUpper, tmaxUpper, vline);
                const auto ponUpper = getterUpper.getPoint(tonUpper);
                // instead of curveUpper.getPoint(tonUpper);
                if (point::isSame(ponUpper, pointMoreLower))
                {
                  return std::tuple{ tonUpper, tofmaxLower , true };
                }
                else
                {
                  return std::nullopt;
                }
              }
              else if (vline > abscissaOfMoreUpper)
              {
                const auto tonLower = getterLower.getTofAbscissa(abscissaOfMoreUpper);
                // instead of curveLower.tofX(tminLower, tmaxLower, pointMoreUpper.x);
                const auto ponLower = getterLower.getPoint(tonLower);
                // instead of curveLower.getPoint(tonLower);
                if (point::isSame(pointMoreUpper, ponLower))
                {
                  return std::tuple{ tofmaxUpper, tonLower, true };
                }
                else
                {
                  return std::nullopt;
                }
              }
              else
              {
                if (point::isSame(pointMoreUpper, pointMoreLower))
                {
                  return std::tuple{ tofmaxUpper, tofmaxLower, true };
                }
                else
                {
                  return std::nullopt;
                }
              }
            }
            else
            {
              const auto tonLower = getterLower.getTofOrdinate(hline);
              // instead of curveLower.tofY(tminLower, tmaxLower, hline);
              const auto ponLower = getterLower.getPoint(tonLower);
              // instead of curveLower.getPoint(tonLower);

              const auto vline = dataGetter::abscissaOf(ponLower);
              // instead of ponLower.x;
              const auto abscissaOfMoreUpper = dataGetter::abscissaOf(pointMoreUpper);
              if (vline < abscissaOfMoreUpper)
              {
                const auto tonUpper = getterUpper.getTofAbscissa(vline);
                // instead of curveUpper.tofX(tminUpper, tmaxUpper, vline);
                const auto ponUpper = getterUpper.getPoint(tonUpper);
                // instead of curveUpper.getPoint(tonUpper);

                pu = ponUpper;
                pl = ponLower;

                tu = tonUpper;
                tl = tonLower;

                continue;
              }
              else if (vline > abscissaOfMoreUpper)
              {
                const auto tonLowerBack = getterLower.getTofAbscissa(abscissaOfMoreUpper);
                // instead of curveLower.tofX(tminLower, tmaxLower, pointMoreUpper.x);
                const auto ponLowerBack = getterLower.getPoint(tonLowerBack);
                // instead of curveLower.getPoint(tonLowerBack);

                if (point::isSame(pointMoreUpper, ponLowerBack))
                {
                  return std::tuple{ tofmaxUpper, tonLowerBack, true };
                }
                else
                {
                  return std::nullopt;
                }
              }
              else
              {
                if (point::isSame(pointMoreUpper, ponLower))
                {
                  return std::tuple{ tofmaxUpper, tonLower, true };
                }
                else
                {
                  return std::nullopt;
                }
              }
            }
          }
        }

        return std::nullopt;
      }

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