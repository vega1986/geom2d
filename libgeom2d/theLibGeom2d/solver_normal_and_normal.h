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
    // Класс для поиска точек пересечения кривых типа Normal на участках монотонности.
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

        const auto tmin1 = std::min(tofmin1, tofmax1);
        const auto tmax1 = std::max(tofmin1, tofmax1);
        const dataGetter getter1{tmin1, tmax1, curve1};

        const auto tmin2 = std::min(tofmin2, tofmax2);
        const auto tmax2 = std::max(tofmin2, tofmax2);
        const dataGetter getter2{ tmin2, tmax2, curve2 };

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
            const auto theResultOfShift = shiftForDiverge(t1, getter1, t2, getter2);
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
              const auto theResultOfShift = shiftForConverge(t2, getter2, t1, getter1);
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
              const auto theResultOfShift = shiftForConverge(t1, getter1, t2, getter2);
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
          const double tofLeft1,
          const dataGetter& getter1,
          const double tofLeft2,
          const dataGetter& getter2
        )
      {
        const auto pointMore1 = getter1.getPointOfMax(); // instead of curve1.getPoint(tofmax1);
        const auto pointMore2 = getter2.getPointOfMax(); // instead of curve2.getPoint(tofmax2);

        auto p1 = getter1.getPoint(tofLeft1);
        auto p2 = getter2.getPoint(tofLeft2);

        auto t1 = tofLeft1;
        auto t2 = tofLeft2;

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
            t1 = getter1.getTofMax(); // instead of tofmax1;
          }
          // update t2
          if (vline < dataGetter::abscissaOf(pointMore2))
          {
            t2 = getter2.getTofAbscissa(vline);
            // instead of curve2.tofX(tmin2, tmax2, vline);
          }
          else
          {
            t2 = getter2.getTofMax(); // instead of tofmax2;
          }
          // update p1 & p2
          p1 = getter1.getPoint(t1);
          p2 = getter2.getPoint(t2);
        }
      }

      static
        std::optional<std::tuple<double, double, bool>>
        shiftForConverge
        (
          const double tofLeftUpper,
          const dataGetter& getterUpper,
          const double tofLeftLower,
          const dataGetter& getterLower
        )
      {
        const auto pointMoreUpper = getterUpper.getPointOfMax();
        const auto pointMoreLower = getterLower.getPointOfMax();

        const auto tofmaxUpper = getterUpper.getTofMax();
        const auto tofmaxLower = getterLower.getTofMax();

        auto pu = getterUpper.getPoint(tofLeftUpper);
        auto pl = getterLower.getPoint(tofLeftLower);

        auto tu = tofLeftUpper;
        auto tl = tofLeftLower;

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
                const auto ponUpper = getterUpper.getPoint(tonUpper);
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
                const auto ponLower = getterLower.getPoint(tonLower);
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
              const auto ponLower = getterLower.getPoint(tonLower);

              const auto vline = dataGetter::abscissaOf(ponLower);
              const auto abscissaOfMoreUpper = dataGetter::abscissaOf(pointMoreUpper);
              if (vline < abscissaOfMoreUpper)
              {
                const auto tonUpper = getterUpper.getTofAbscissa(vline);
                const auto ponUpper = getterUpper.getPoint(tonUpper);

                pu = ponUpper;
                pl = ponLower;

                tu = tonUpper;
                tl = tonLower;

                continue;
              }
              else if (vline > abscissaOfMoreUpper)
              {
                const auto tonLowerBack = getterLower.getTofAbscissa(abscissaOfMoreUpper);
                const auto ponLowerBack = getterLower.getPoint(tonLowerBack);

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