#include "solver_normal_and_normal_along_y.h"

#include <algorithm>

//*********************************************************************************************************************

std::vector<geom2d::IntersecctionSolutionType>
geom2d::solver_normal_and_normal_along_y::solver::execute() const
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
      // в этой ветке утилизирован случай, когда p1.x == p2.x
      // так как если это так, то dist(p1, p2) <= 0.1 * tolerance точки
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
      if (p2.x > p1.x)
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
      else if (p1.x > p2.x)
      {
        // 1 - upper curve
        // 2 - lower curve
        const auto theResultOfShift = shiftForConverge(t1, tofmax1, curve1, t2, tofmax2, curve2);
        if (not theResultOfShift)
        {
          return result;
        }
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
        throw std::logic_error("impossible situation in function geom2d::solver_normal_and_normal_along_y::solver::execute");
      }
    }
  }
  return result;
}

//*********************************************************************************************************************

std::optional<std::pair<double, double>>
geom2d::solver_normal_and_normal_along_y::solver::shiftForDiverge
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

  const auto tmin2 = std::min(tofmin2, tofmax2);
  const auto tmax2 = std::max(tofmin2, tofmax2);

  const auto pointMore1 = curve1.getPoint(tofmax1);
  const auto pointMore2 = curve2.getPoint(tofmax2);

  auto p1 = curve1.getPoint(tofmin1);
  auto p2 = curve2.getPoint(tofmin2);

  auto t1 = tofmin1;
  auto t2 = tofmin2;

  auto theShift = math::tolerance::tolPoint;
  const auto initial_vline = std::max(p1.y, p2.y);

  while (true)
  {
    if (not point::isSame(p1, p2))
    {
      return std::pair{ t1, t2 };
    }
    const double vline = initial_vline + theShift;
    theShift *= 2;
    if ((vline > pointMore1.y) or (vline > pointMore2.y))
    {
      return std::nullopt;
    }
    if (vline < pointMore1.y) t1 = curve1.tofY(tmin1, tmax1, vline);
    else                       t1 = tofmax1;

    if (vline < pointMore2.y) t2 = curve2.tofY(tmin2, tmax2, vline);
    else                       t2 = tofmax2;

    // update p1 & p2
    p1 = curve1.getPoint(t1);
    p2 = curve2.getPoint(t2);
  }
}

//*********************************************************************************************************************

std::optional<std::tuple<double, double, bool>>
geom2d::solver_normal_and_normal_along_y::solver::shiftForConverge
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

  const auto tminLower = std::min(tofminLower, tofmaxLower);
  const auto tmaxLower = std::max(tofminLower, tofmaxLower);

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
      const auto hline = pu.x;

      if (hline >= pointMoreLower.x)
      {
        const auto vline = pointMoreLower.y;
        if (vline < pointMoreUpper.y)
        {
          const auto tonUpper = curveUpper.tofY(tminUpper, tmaxUpper, vline);
          const auto ponUpper = curveUpper.getPoint(tonUpper);
          if (point::isSame(ponUpper, pointMoreLower))
          {
            return std::tuple{ tonUpper, tofmaxLower , true };
          }
          else
          {
            return std::nullopt;
          }
        }
        else if (vline > pointMoreUpper.y)
        {
          const auto tonLower = curveLower.tofY(tminLower, tmaxLower, pointMoreUpper.y);
          const auto ponLower = curveLower.getPoint(tonLower);
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
        const auto tonLower = curveLower.tofX(tminLower, tmaxLower, hline);
        const auto ponLower = curveLower.getPoint(tonLower);

        const auto vline = ponLower.y;

        if (vline < pointMoreUpper.y)
        {
          const auto tonUpper = curveUpper.tofY(tminUpper, tmaxUpper, vline);
          const auto ponUpper = curveUpper.getPoint(tonUpper);

          pu = ponUpper;
          pl = ponLower;

          tu = tonUpper;
          tl = tonLower;

          continue;
        }
        else if (vline > pointMoreUpper.y)
        {
          const auto tonLowerBack = curveLower.tofY(tminLower, tmaxLower, pointMoreUpper.y);
          const auto ponLowerBack = curveLower.getPoint(tonLowerBack);

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
}

//*********************************************************************************************************************
