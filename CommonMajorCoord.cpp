#include "CommonMajorCoord.h"
#include "findFunctionRoots.h"

geom2d::CommonMajorCoord::CommonMajorCoord
(
  const double tlow1,
  const double thigh1,
  const baseCurve& crv1,
  const double tlow2,
  const double thigh2,
  const baseCurve& crv2
)
  :
  tmin1(tlow1),
  tmax1(thigh1),
  curve1(crv1),
  tmin2(tlow2),
  tmax2(thigh2),
  curve2(crv2)
{
}

const
  std::tuple<double, double, geom2d::point, double, geom2d::point>
    geom2d::CommonMajorCoord::getMinCoords(const CurveMajorHelper& getter1, const CurveMajorHelper& getter2)
{
  if (getter2.getCoordMin() < getter1.getCoordMin())
  {
    const auto commonMin = getter1.getCoordMin();
    const auto tofmin1 = getter1.getTofCoordMin();
    const auto pntofmin1 = getter1.getPointMin();
    auto func = [&getter2, varo = commonMin](const double t) -> double
    {
      return getter2.getCoord(t) - varo;
    };
    const auto tofmin2 = math::findUniqueFunctionRoot(tmin2, tmax2, func);
    const auto pntofmin2 = curve2.getPoint(tofmin2);
    return std::tuple{commonMin, tofmin1, pntofmin1, tofmin2, pntofmin2};
  }
  else if (getter1.getCoordMin() < getter2.getCoordMin())
  {
    const auto commonMin = getter2.getCoordMin();
    const auto tofmin2 = getter2.getTofCoordMin();
    const auto pntofmin2 = getter2.getPointMin();
    auto func = [&getter1, varo = commonMin](const double t) -> double
    {
      return getter1.getCoord(t) - varo;
    };
    const auto tofmin1 = math::findUniqueFunctionRoot(tmin1, tmax1, func);
    const auto pntofmin1 = curve1.getPoint(tofmin1);
    return std::tuple{ commonMin, tofmin1, pntofmin1, tofmin2, pntofmin2 };
  }
  else
  {
    const auto commonMin = getter1.getCoordMin();
    const auto tofmin1 = getter1.getTofCoordMin();
    const auto pntofmin1 = getter1.getPointMin();
    const auto tofmin2 = getter2.getTofCoordMin();
    const auto pntofmin2 = getter2.getPointMin();
    return std::tuple{ commonMin, tofmin1, pntofmin1, tofmin2, pntofmin2 };
  }
}

const
  std::tuple<double, double, geom2d::point, double, geom2d::point>
    geom2d::CommonMajorCoord::getMaxCoords(const CurveMajorHelper& getter1, const CurveMajorHelper& getter2)
{
  if (getter2.getCoordMax() < getter1.getCoordMax())
  {
    const auto commonMax = getter2.getCoordMax();
    const auto tofmax2 = getter2.getTofCoordMax();
    const auto pntofmax2 = getter2.getPointMax();
    auto func = [&getter1, varo = commonMax](const double t) -> double
    {
      return getter1.getCoord(t) - varo;
    };
    const auto tofmax1 = math::findUniqueFunctionRoot(tmin1, tmax1, func);
    const auto pntofmax1 = curve1.getPoint(tofmax1);
    return std::tuple{ commonMax, tofmax1, pntofmax1, tofmax2, pntofmax2 };
  }
  else if (getter1.getCoordMax() < getter2.getCoordMax())
  {
    const auto commonMax = getter1.getCoordMax();
    const auto tofmax1 = getter1.getTofCoordMax();
    const auto pntofmax1 = getter1.getPointMax();
    auto func = [&getter2, varo = commonMax](const double t) -> double
    {
      return getter2.getCoord(t) - varo;
    };
    const auto tofmax2 = math::findUniqueFunctionRoot(tmin2, tmax2, func);
    const auto pntofmax2 = curve2.getPoint(tofmax2);
    return std::tuple{ commonMax, tofmax1, pntofmax1, tofmax2, pntofmax2 };
  }
  else
  {
    const auto commonMax = getter1.getCoordMax();
    const auto tofmax1 = getter1.getTofCoordMax();
    const auto pntofmax1 = getter1.getPointMax();
    const auto tofmax2 = getter2.getTofCoordMax();
    const auto pntofmax2 = getter2.getPointMax();
    return std::tuple{ commonMax, tofmax1, pntofmax1, tofmax2, pntofmax2};
  }
}

void geom2d::XCommonMajorCoord::setXmin()
{
  const XofCurvePointGetter getter1{ tmin1, tmax1, curve1 };
  const XofCurvePointGetter getter2{ tmin2, tmax2, curve2 };

  const auto [aXmin, aTofxmin1, pnt1, aTofxmin2, pnt2] = CommonMajorCoord::getMinCoords(getter1, getter2);
  xmin = aXmin;
  tofxmin1 = aTofxmin1;
  pntofxmin1 = pnt1;
  tofxmin2 = aTofxmin2;
  pntofxmin2 = pnt2;
}

void geom2d::XCommonMajorCoord::setXmax()
{
  const XofCurvePointGetter getter1{ tmin1, tmax1, curve1 };
  const XofCurvePointGetter getter2{ tmin2, tmax2, curve2 };

  const auto [aXmax, aTofxmax1, pnt1, aTofxmax2, pnt2] = CommonMajorCoord::getMaxCoords(getter1, getter2);
  xmax = aXmax;
  tofxmax1 = aTofxmax1;
  pntofxmax1 = pnt1;
  tofxmax2 = aTofxmax2;
  pntofxmax2 = pnt2;
}

void geom2d::YCommonMajorCoord::setYmin()
{
  const YofCurvePointGetter getter1{ tmin1, tmax1, curve1 };
  const YofCurvePointGetter getter2{ tmin2, tmax2, curve2 };

  const auto [aYmin, aTofymin1, pnt1, aTofymin2, pnt2] = CommonMajorCoord::getMinCoords(getter1, getter2);
  ymin = aYmin;
  tofymin1 = aTofymin1;
  pntofymin1 = pnt1;
  tofymin2 = aTofymin2;
  pntofymin2 = pnt2;
}

void geom2d::YCommonMajorCoord::setYmax()
{
  const XofCurvePointGetter getter1{ tmin1, tmax1, curve1 };
  const XofCurvePointGetter getter2{ tmin2, tmax2, curve2 };

  const auto [aYmax, aTofymax1, pnt1, aTofymax2, pnt2] = CommonMajorCoord::getMaxCoords(getter1, getter2);
  ymax = aYmax;
  tofymax1 = aTofymax1;
  pntofymax1 = pnt1;
  tofymax2 = aTofymax2;
  pntofymax2 = pnt2;
}
