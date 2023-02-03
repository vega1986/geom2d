#include "curveIntersector.h"
#include "findFunctionRoots.h"
#include "avector.h"
#include "segmentCurve.h"
#include "amath.h"
#include "Axis.h"
#include "StatOfCurvePiece.h"
#include "CommonRangeHelper.h"
#include "solver_unique_intersection.h"
#include "solver_point_and_alongaxis.h"
#include "solver_points_and_alongaxis.h"
#include "solver_point_and_point.h"
#include "solver_normal_and_normal.h"
#include "solver_screen_and_screen.h"

#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>

namespace geom2d
{
  // Исследуем решение, если первая и вторая кривая - PlatoX или PlatoY
  template <CurveDataGetter dataGetter>
  std::optional<IntersecctionSolutionType>
    execParallelPlatos
    (
      const double tmin1,
      const double tmax1,
      const baseCurve& curvePlato1,
      const double tmin2,
      const double tmax2,
      const baseCurve& curvePlato2
    )
  {
    double currentDistance = math::infinite::distance;
    double tof1 = 0.0;
    point pof1;
    double tof2 = 0.0;
    point pof2;
    bool theMinDistanceFound = false;

    // рассматриваем точки 1-ой кривой и их положение по отношению ко кривой 2
    std::initializer_list ts1{ tmin1 , tmax1 };
    for (const auto at : ts1)
    {
      using namespace point_and_curve_alongaxis;
      const point P = curvePlato1.getPoint(at);

      solver<dataGetter> solv{ P , tmin2 , tmax2 , curvePlato2 };
      const auto result = solv.execute();

      if (not result) continue;
      const auto theT = result.value();
      const auto thePoint = curvePlato2.getPoint(theT);
      const auto dist = point::distance(P, thePoint);
      if (dist < currentDistance)
      {
        theMinDistanceFound = true;
        currentDistance = dist;
        tof1 = at;
        pof1 = P;
        tof2 = theT;
        pof2 = thePoint;
      }
    }

    // рассматриваем точки 2-ой кривой и их положение по отношению ко кривой 1
    std::initializer_list ts2{ tmin2 , tmax2 };
    for (const auto at : ts2)
    {
      using namespace point_and_curve_alongaxis;
      const point P = curvePlato2.getPoint(at);

      solver<dataGetter> solv{ P , tmin1 , tmax1 , curvePlato1 };
      const auto result = solv.execute();

      if (not result) continue;
      const auto theT = result.value();
      const auto thePoint = curvePlato1.getPoint(theT);
      const auto dist = point::distance(P, thePoint);
      if (dist < currentDistance)
      {
        theMinDistanceFound = true;
        currentDistance = dist;
        tof2 = at;
        pof2 = P;
        tof1 = theT;
        pof1 = thePoint;
      }
    }

    if (theMinDistanceFound and point::isSame(pof1, pof2))
    {
      return IntersecctionSolutionType{ 0.5 * (pof1 + pof2), tof1 , tof2 };
    }
    return std::nullopt;
  }
}

//-----------------------------------------------------------------------------

void geom2d::curveIntersector::fulfill()
{
  curveAnalizerBase::perform();
}

//-----------------------------------------------------------------------------

void geom2d::curveIntersector::postProcessing()
{
  excludeDuplicatesFromSolution();
}

//-----------------------------------------------------------------------------

void geom2d::curveIntersector::excludeDuplicatesFromSolution()
{
  auto& pnts = solutionPoints;
  auto& t1s = solutionParameterOnCurve1;
  auto& t2s = solutionParameterOnCurve2;

  if (pnts.size() < 2) return;

  std::vector<bool> forExclude;
  forExclude.assign(pnts.size(), false);

  auto need2exclude = [&pnts, &t1s, &t2s](const size_t j) -> bool
  {
    const point current{ pnts[j] };
    for (size_t i = j + 1; i < pnts.size(); ++i)
    {
      if (point::isSame(current, pnts[i])) return true;
    }
    return false;
  };

  for (size_t k = 0; k < pnts.size() - 1; ++k)
  {
    forExclude[k] = need2exclude(k);
  }

  std::remove_const_t<std::remove_reference_t<decltype(pnts)>> newPnts;
  std::remove_const_t<std::remove_reference_t<decltype(t1s)>> newT1s;
  std::remove_const_t<std::remove_reference_t<decltype(t2s)>> newT2s;

  for (size_t k = 0; k < pnts.size(); ++k)
  {
    if (forExclude[k]) continue;

    newPnts.push_back(pnts[k]);
    newT1s.push_back(t1s[k]);
    newT2s.push_back(t2s[k]);
  }

  pnts = std::move(newPnts);
  t1s = std::move(newT1s);
  t2s = std::move(newT2s);
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performScreen1Screen2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // кривая 1 - Screen
  // кривая 2 - Screen
  const auto result =
    execScreenAndScreen(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (not result.empty())
  {
    for (const auto [intPoint, t1, t2] : result)
    {
      solutionPoints.push_back(intPoint);
      solutionParameterOnCurve1.push_back(t1);
      solutionParameterOnCurve2.push_back(t2);
    }
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performScreen1Normal2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // кривая 1 - Screen
  // кривая 2 - Normal
  const auto result =
    execNormalAndScreen(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (result)
  {
    const auto [intPoint, t2, t1] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performScreen1PlatoX2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // первая кривая - Any (Screen)
  // вторая кривая - PlatoX
  const auto result =
    execPlatoXAndAny(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (result)
  {
    const auto [intPoint, t2, t1] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performScreen1PlatoY2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // первая кривая - Any (Screen)
  // вторая кривая - PlatoY
  const auto result =
    execPlatoYAndAny(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (result)
  {
    const auto [intPoint, t2, t1] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performScreen1Point2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto result =
    execPointAndAny(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (result)
  {
    const auto [intPoint, t2, t1] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void 
  geom2d::curveIntersector::performNormal1Screen2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // кривая 1 - Normal
  // кривая 2 - Screen
  const auto result =
    execNormalAndScreen(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (result)
  {
    const auto [intPoint, t1, t2] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performNormal1Normal2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // кривая 1 - Normal
        // кривая 2 - Normal
  const auto result =
    execNormalAndNormal(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (not result.empty())
  {
    for (const auto [intPoint, t1, t2] : result)
    {
      solutionPoints.push_back(intPoint);
      solutionParameterOnCurve1.push_back(t1);
      solutionParameterOnCurve2.push_back(t2);
    }
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performNormal1PlatoX2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // первая кривая - Any (Normal)
  // вторая кривая - PlatoX
  const auto result =
    execPlatoXAndAny(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (result)
  {
    const auto [intPoint, t2, t1] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performNormal1PlatoY2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // первая кривая - Any (Normal)
  // вторая кривая - PlatoY
  const auto result =
    execPlatoYAndAny(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (result)
  {
    const auto [intPoint, t2, t1] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performNormal1Point2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // кривая 1 - Normal
  // кривая 2 - Point
  const auto result =
    execPointAndAny(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (result)
  {
    const auto [intPoint, t2, t1] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPlatoX1Any2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // первая кривая - PlatoX
  // вторая кривая - Any (Screen or Normal)
  const auto result =
    execPlatoXAndAny(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (result)
  {
    const auto [intPoint, t1, t2] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPlatoX1PlatoX2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // обе кривые - PlatoX
  const auto result =
    execPlatoXAndPlatoX(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (result)
  {
    const auto [intPoint, t1, t2] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPlatoX1PlatoY2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // первая кривая - PlatoX
      // вторая кривая - PlatoY
  const auto result =
    execPlatoXAndPlatoY(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (result)
  {
    const auto [intPoint, t1, t2] = result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPlatoX1Point2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // кривая 1 - PlatoX
  // кривая 2 - Point
  const auto result =
    execPointAndPlatoX(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (result)
  {
    const auto [intPoint, t2, t1] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPlatoY1Any2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // первая кривая - PlatoY
  // вторая кривая - Any (Screen or Normal)
  const auto result =
    execPlatoYAndAny(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (result)
  {
    const auto [intPoint, t1, t2] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPlatoY1PlatoX2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // кривая 1 - PlatoY
  // кривая 2 - PlatoX
  const auto result =
    execPlatoXAndPlatoY(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (result)
  {
    const auto [intPoint, t2, t1] = result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPlatoY1PlatoY2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // обе кривые - PlatoY
  const auto result =
    execPlatoYAndPlatoY(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (result)
  {
    const auto [intPoint, t1, t2] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPlatoY1Point2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // первая кривая - PlatoY
  // вторая кривая - Point
  const auto result =
    execPointAndPlatoY(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (result)
  {
    const auto [intPoint, t2, t1] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPoint1Any2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto result = execPointAndAny(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (result)
  {
    const auto [intPoint, t1, t2] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPoint1PlatoX2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // первая кривая - Point
  // вторая кривая - PlatoX
  const auto result =
    execPointAndPlatoX(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (result)
  {
    const auto [intPoint, t1, t2] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPoint1PlatoY2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // первая кривая - Point
  // вторая кривая - PlatoY
  const auto result =
    execPointAndPlatoY(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (result)
  {
    const auto [intPoint, t1, t2] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

void
  geom2d::curveIntersector::performPoint1Point2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto result =
    execPointAndPoint(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (result)
  {
    const auto [intPoint, t1, t2] =
      result.value();
    solutionPoints.push_back(intPoint);
    solutionParameterOnCurve1.push_back(t1);
    solutionParameterOnCurve2.push_back(t2);
  }
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execPointAndPoint
  (
    const double tminOfPoint1,
    const double tmaxOfPoint1,
    const baseCurve& curvePoint1,
    const double tminOfPoint2,
    const double tmaxOfPoint2,
    const baseCurve& curvePoint2
  )
{
  // рассматриваем всевозможные пары точек
  std::initializer_list<std::pair<double, double>> listOfPairsOfT
  {
    {tminOfPoint1, tminOfPoint2},
    {tminOfPoint1, tmaxOfPoint2},
    {tmaxOfPoint1, tminOfPoint2},
    {tmaxOfPoint1, tmaxOfPoint2}
  };

  double currentDistance = math::infinite::distance;
  double tofExtrema1 = 0.0;
  point pofExtrema1;
  double tofExtrema2 = 0.0;
  point pofExtrema2;
  bool theMinimalFound = false;
  for (const auto [t1, t2] : listOfPairsOfT)
  {
    const auto p1 = curvePoint1.getPoint(t1);
    const auto p2 = curvePoint2.getPoint(t2);
    const auto dist12 = point::distance(p1, p2);
    if (dist12 < currentDistance)
    {
      theMinimalFound = true;
      currentDistance = dist12;
      tofExtrema1 = t1;
      pofExtrema1 = p1;
      tofExtrema2 = t2;
      pofExtrema2 = p2;
    }
  }

  if (not theMinimalFound) return std::nullopt;

  using namespace geom2d::point_and_point;
  solver solv{ pofExtrema1 , pofExtrema2 };
  if (solv.execute())
  {
    return IntersecctionSolutionType{ 0.5 * (pofExtrema1 + pofExtrema2), tofExtrema1 , tofExtrema2 };
  }
  return std::nullopt;
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
geom2d::curveIntersector::execPointAndPlatoX(
  const double tminOfPoint,
  const double tmaxOfPoint,
  const baseCurve& curvePoint,
  const double tminOfPlatoX,
  const double tmaxOfPlatoX,
  const baseCurve& curvePlatoX)
{
  std::initializer_list tsofPoint{ tminOfPoint , tmaxOfPoint };

  double currentDistance = math::infinite::distance;
  double tofPoint = 0.0;
  point pofPoint;
  double tofPlato = 0.0;
  point pofPlato;
  bool theMinDistanceFound = false;
  for (const auto at : tsofPoint)
  {
    using namespace point_and_curve_alongaxis;
    const point P = curvePoint.getPoint(at);
    
    solver<DataGetterOfY> solv{ P , tminOfPlatoX , tmaxOfPlatoX , curvePlatoX };
    const auto result = solv.execute();

    if (not result) continue;
    const auto theT = result.value();
    const auto thePoint = curvePlatoX.getPoint(theT);
    const auto dist = point::distance(P, thePoint);
    if (dist < currentDistance)
    {
      theMinDistanceFound = true;
      currentDistance = dist;
      tofPoint = at;
      pofPoint = P;
      tofPlato = theT;
      pofPlato = thePoint;
    }
  }
  
  if (theMinDistanceFound and point::isSame(pofPoint, pofPlato))
  {
    return IntersecctionSolutionType{0.5 * (pofPoint + pofPlato), tofPoint , tofPlato };
  }
  return std::nullopt;
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execPointAndPlatoY
  (
    const double tminOfPoint,
    const double tmaxOfPoint,
    const baseCurve& curvePoint,
    const double tminOfPlatoY,
    const double tmaxOfPlatoY,
    const baseCurve& curvePlatoY
  )
{
  std::initializer_list tsofPoint{ tminOfPoint , tmaxOfPoint };

  double currentDistance = math::infinite::distance;
  double tofPoint = 0.0;
  point pofPoint;
  double tofPlato = 0.0;
  point pofPlato;
  bool theMinDistanceFound = false;
  for (const auto at : tsofPoint)
  {
    using namespace point_and_curve_alongaxis;
    const point P = curvePoint.getPoint(at);

    solver<DataGetterOfX> solv{ P , tminOfPlatoY , tmaxOfPlatoY , curvePlatoY };
    const auto result = solv.execute();

    if (not result) continue;
    const auto theT = result.value();
    const auto thePoint = curvePlatoY.getPoint(theT);
    const auto dist = point::distance(P, thePoint);
    if (dist < currentDistance)
    {
      theMinDistanceFound = true;
      currentDistance = dist;
      tofPoint = at;
      pofPoint = P;
      tofPlato = theT;
      pofPlato = thePoint;
    }
  }

  if (theMinDistanceFound and point::isSame(pofPoint, pofPlato))
  {
    return IntersecctionSolutionType{ 0.5 * (pofPoint + pofPlato), tofPoint , tofPlato };
  }
  return std::nullopt;
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execPointAndAny
  (
    const double tminOfPoint,
    const double tmaxOfPoint,
    const baseCurve& curvePoint,
    const double tminOfAny,
    const double tmaxOfAny,
    const baseCurve& curveAny
  )
{
  std::initializer_list listofTofPoint{ tminOfPoint , tmaxOfPoint };

  double currentDistance = math::infinite::distance;

  double tpointResult = 0.0;
  point ppointResult;
  
  double tanyResult = 0.0;
  point panyResult;
  
  bool solutionFound = false;
  // вдоль оси X
  for (const auto t : listofTofPoint)
  {
    const auto p = curvePoint.getPoint(t);
    using namespace point_and_curve_alongaxis;
    solver<DataGetterOfX> solv{ p, tminOfAny, tmaxOfAny, curveAny };
    const auto result = solv.execute();

    if (not result) continue;

    const auto tonAny = result.value();
    const auto ponAny = curveAny.getPoint(tonAny);

    const auto dist = point::distance(p, ponAny);

    if (dist < currentDistance)
    {
      solutionFound = true;
      currentDistance = dist;
      tpointResult = t;
      ppointResult = p;
      tanyResult = tonAny;
      panyResult = ponAny;
    }
  }
  // вдоль оси Y
  for (const auto t : listofTofPoint)
  {
    const auto p = curvePoint.getPoint(t);
    using namespace point_and_curve_alongaxis;
    solver<DataGetterOfY> solv{ p, tminOfAny, tmaxOfAny, curveAny };
    const auto result = solv.execute();

    if (not result) continue;

    const auto tonAny = result.value();
    const auto ponAny = curveAny.getPoint(tonAny);

    const auto dist = point::distance(p, ponAny);

    if (dist < currentDistance)
    {
      solutionFound = true;
      currentDistance = dist;
      tpointResult = t;
      ppointResult = p;
      tanyResult = tonAny;
      panyResult = ponAny;
    }
  }
  if (solutionFound)
  {
    return IntersecctionSolutionType{ 0.5 * (ppointResult + panyResult), tpointResult , tanyResult };
  }
  return std::nullopt;  
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execPlatoYAndPlatoY
  (
    const double tmin1,
    const double tmax1,
    const baseCurve& curvePlatoY1,
    const double tmin2,
    const double tmax2,
    const baseCurve& curvePlatoY2
  )
{
  return geom2d::execParallelPlatos<DataGetterOfX>(tmin1, tmax1, curvePlatoY1, tmin2, tmax2, curvePlatoY2);
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execPlatoXAndPlatoX
  (
    const double tmin1,
    const double tmax1,
    const baseCurve& curvePlatoX1,
    const double tmin2,
    const double tmax2,
    const baseCurve& curvePlatoX2
  )
{
  return geom2d::execParallelPlatos<DataGetterOfY>(tmin1, tmax1, curvePlatoX1, tmin2, tmax2, curvePlatoX2);
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execPlatoXAndPlatoY
  (
    const double tminOfPlatox,
    const double tmaxOfPlatox,
    const baseCurve& curvePlatoX,
    const double tminOfPlatoy,
    const double tmaxOfPlatoy,
    const baseCurve& curvePlatoY)
{
  return exec_Any_and_Any_Unique_Intersection(tminOfPlatox, tmaxOfPlatox, curvePlatoX, tminOfPlatoy, tmaxOfPlatoy, curvePlatoY);
}

//-------------------------------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
geom2d::curveIntersector::findUniqueIntersection
(
  const double tmin1,
  const double tmax1,
  const baseCurve& curve1,
  const double tmin2,
  const double tmax2,
  const baseCurve& curve2
)
{
  const auto p1 = curve1.getPoint(tmin1);
  const auto q1 = curve1.getPoint(tmax1);
  const aabb aabb1{ std::initializer_list{p1, q1} };

  const auto p2 = curve2.getPoint(tmin2);
  const auto q2 = curve2.getPoint(tmax2);
  const aabb aabb2{ std::initializer_list{p2, q2} };

  using StatCurveType = std::tuple<int, Axis, double>;
  StatCurveType curve1AlongX{ int{1}, Axis::X, aabb1.lengthx() };
  StatCurveType curve1AlongY{ int{1}, Axis::Y, aabb1.lengthy() };
  StatCurveType curve2AlongX{ int{2}, Axis::X, aabb2.lengthx() };
  StatCurveType curve2AlongY{ int{2}, Axis::Y, aabb2.lengthy() };
  auto statCurveLess = [](const StatCurveType& sc1, const StatCurveType& sc2)
  {
    const auto [i1, axis1, len1] = sc1;
    const auto [i2, axis2, len2] = sc2;
    return len1 < len2;
  };
  const auto [maxCurveId, maxAxis, maxLength] = 
    std::max(std::initializer_list{ curve1AlongX , curve1AlongY , curve2AlongX , curve2AlongY }, statCurveLess);
  if (maxCurveId == 1)
  {
    const auto result = (maxAxis == Axis::X) ?
      findUniqueIntersectionRefAlongX(tmin1, tmax1, curve1, tmin2, tmax2, curve2)
      :
      findUniqueIntersectionRefAlongY(tmin1, tmax1, curve1, tmin2, tmax2, curve2);
    return result;
  }
  else
  {
    const auto result = (maxAxis == Axis::X) ?
      findUniqueIntersectionRefAlongX(tmin2, tmax2, curve2, tmin1, tmax1, curve1)
      :
      findUniqueIntersectionRefAlongY(tmin2, tmax2, curve2, tmin1, tmax1, curve1);
    if (result)
    {
      const auto [dist, t2, t1] = result.value();
      return std::tuple{ dist, t1, t2 };
    }
    else
    {
      return std::nullopt;
    }
  }
  return std::nullopt;
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::findUniqueIntersectionRefAlongX
  (
    const double trefmin,
    const double trefmax,
    const baseCurve& referenceCurve,
    const double tothmin,
    const double tothmax,
    const baseCurve& otherCurve
  )
{
  using namespace geom2d::uniqueIntersection;
  solver<DataGetterOfX> the_solver{trefmin, trefmax, referenceCurve, tothmin, tothmax, otherCurve };

  return the_solver.execute();
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::findUniqueIntersectionRefAlongY
  (
    const double trefmin,
    const double trefmax,
    const baseCurve& referenceCurve,
    const double tothmin,
    const double tothmax,
    const baseCurve& otherCurve
  )
{
  using namespace geom2d::uniqueIntersection;
  solver<DataGetterOfY> the_solver{ trefmin, trefmax, referenceCurve, tothmin, tothmax, otherCurve };

  return the_solver.execute();
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execPlatoXAndAny
  (
    const double tminOfPlatox,
    const double tmaxOfPlatox,
    const baseCurve& curvePlatoX,
    const double tminOfAny,
    const double tmaxOfAny,
    const baseCurve& curveAny
  )
{
  return exec_Any_and_Any_Unique_Intersection(tminOfPlatox, tmaxOfPlatox, curvePlatoX, tminOfAny, tmaxOfAny, curveAny);
#if 0
  std::initializer_list listOfPlatoXt{ tminOfPlatox , tmaxOfPlatox };
  // Сначала тестируем конечные точки кривой PlatoX - лежат ли они на Any?
    // trace along X axis of curveAny
  for (const auto t : listOfPlatoXt)
  {
    const auto p = curvePlatoX.getPoint(t);
    using namespace geom2d::point_and_curve_alongaxis;

    solver<DataGetterOfX> solv{ p , tminOfAny , tmaxOfAny , curveAny };
    const auto result = solv.execute();
    if (not result) continue;

    const auto tofAny = result.value();
    const auto pofAny = curveAny.getPoint(tofAny);
    return IntersecctionSolutionType{ 0.5 * (p + pofAny), t, tofAny };
  }
    // trace along Y axis of curveAny
  for (const auto t : listOfPlatoXt)
  {
    const auto p = curvePlatoX.getPoint(t);
    using namespace geom2d::point_and_curve_alongaxis;

    solver<DataGetterOfY> solv{ p, tminOfAny , tmaxOfAny , curveAny };
    const auto result = solv.execute();
    if (not result) continue;

    const auto tofAny = result.value();
    const auto pofAny = curveAny.getPoint(tofAny);
    return IntersecctionSolutionType{ 0.5 * (p + pofAny), t, tofAny };
  }


  // А теперь тестируем конечные точки кривой Any - лежат ли они на PlatoX?
  std::initializer_list listOfAnyt{ tminOfAny , tmaxOfAny };
  for (const auto t : listOfAnyt)
  {
    const auto p = curveAny.getPoint(t);
    using namespace geom2d::point_and_curve_alongaxis;

    solver<DataGetterOfY> solv{ p, tminOfPlatox, tmaxOfPlatox, curvePlatoX };
    const auto result = solv.execute();
    if (not result) continue;

    const auto tofPlatox = result.value();
    const auto pofPlatox = curvePlatoX.getPoint(tofPlatox);
    return IntersecctionSolutionType{ 0.5 * (p + pofPlatox), tofPlatox, t };
  }

  // далее проверяем полноценное пересечение или его отсутствие (граничные случаи утилизированы выше)
  return findUniqueIntersection(tminOfPlatox, tmaxOfPlatox, curvePlatoX, tminOfAny, tmaxOfAny, curveAny);
#endif
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execPlatoYAndAny
  (
    const double tminOfPlatoy,
    const double tmaxOfPlatoy,
    const baseCurve& curvePlatoY,
    const double tminOfAny,
    const double tmaxOfAny,
    const baseCurve& curveAny
  )
{
  return exec_Any_and_Any_Unique_Intersection(tminOfPlatoy, tmaxOfPlatoy, curvePlatoY, tminOfAny, tmaxOfAny, curveAny);
}

//-----------------------------------------------------------------------------

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execNormalAndScreen
  (
    const double tminOfNormal,
    const double tmaxOfNormal,
    const baseCurve& curveNormal,
    const double tminOfScreen,
    const double tmaxOfScreen,
    const baseCurve& curveScreen
  )
{
  return exec_Any_and_Any_Unique_Intersection(tminOfNormal, tmaxOfNormal, curveNormal, tminOfScreen, tmaxOfScreen, curveScreen);
}

//-----------------------------------------------------------------------------

std::vector<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execNormalAndNormal
  (
    const double tmin1,
    const double tmax1,
    const baseCurve& curve1,
    const double tmin2,
    const double tmax2,
    const baseCurve& curve2)
{
  std::vector<geom2d::IntersecctionSolutionType> result;
  // начинаем с нижних точек обеих кривых и идём к точке пересечения методом ступенек
  const auto P1 = curve1.getPoint(tmin1);
  const auto Q1 = curve1.getPoint(tmax1);
  StatOfCurvePiece scp1{ tmin1, P1, tmax1, Q1 };

  const auto P2 = curve2.getPoint(tmin2);
  const auto Q2 = curve2.getPoint(tmax2);
  StatOfCurvePiece scp2{ tmin2, P2, tmax2, Q2 };

  // сначала утилизируем тривиальные случаи
  if (scp1.pointOfxmin().x >= scp2.pointOfxmax().x)
  {
    // проверяем эти точки на совпадение
    const auto t1 = scp1.tOfxmin();
    const auto p1 = scp1.pointOfxmin();

    const auto t2 = scp2.tOfxmax();
    const auto p2 = scp2.pointOfxmax();

    if (point::isSame(p1, p2))
    {
      result.push_back(IntersecctionSolutionType{ 0.5 * (p1 + p2), t1, t2 });
    }
    return result;
  }
  else if (scp2.pointOfxmin().x >= scp1.pointOfxmax().x)
  {
    // проверяем эти точки на совпадение
    const auto t1 = scp1.tOfxmax();
    const auto p1 = scp1.pointOfxmax();

    const auto t2 = scp2.tOfxmin();
    const auto p2 = scp2.pointOfxmin();

    if (point::isSame(p1, p2))
    {
      result.push_back(IntersecctionSolutionType{ 0.5 * (p1 + p2), t1, t2 });
    }
    return result;
  }
  else if (scp1.pointOfymin().y >= scp2.pointOfymax().y)
  {
    const auto t1 = scp1.tOfymin();
    const auto p1 = scp1.pointOfymin();

    const auto t2 = scp2.tOfymax();
    const auto p2 = scp2.pointOfymax();

    if (point::isSame(p1, p2))
    {
      result.push_back(IntersecctionSolutionType{0.5 * (p1 + p2), t1, t2});
    }
    return result;
  }
  else if (scp2.pointOfymin().y >= scp1.pointOfymax().y)
  {
    const auto t1 = scp1.tOfymax();
    const auto p1 = scp1.pointOfymax();

    const auto t2 = scp2.tOfymin();
    const auto p2 = scp2.pointOfymin();

    if (point::isSame(p1, p2))
    {
      result.push_back(IntersecctionSolutionType{ 0.5 * (p1 + p2), t1, t2 });
    }
    return result;
  }
  else
  {
    // если мы пришли сюда, значит есть пересечение ОДЗ двух кривых как по X, так и по Y
    // находим общее ОДЗ по оси X
    DataGetterOfX getter1alongx{ tmin1, tmax1, curve1 };
    DataGetterOfX getter2alongx{ tmin2, tmax2, curve2 };

    ///////////
    //       //
    // MIN X //
    //       //
    ///////////
    auto [
      commonXmin,
      tofCommonXmin1,
      pointofCommonXmin1,
      tofCommonXmin2,
      pointofCommonXmin2
    ] = CommonRangeHelper::ofLowest(getter1alongx, getter2alongx);
    ///////////
    //       //
    // MAX X //
    //       //
    ///////////
    auto [
      commonXmax,
      tofCommonXmax1,
      pointofCommonXmax1,
      tofCommonXmax2,
      pointofCommonXmax2
    ] = CommonRangeHelper::ofHighest(getter1alongx, getter2alongx);
    
    const auto commonRangeX = commonXmax - commonXmin;
    // находим общее ОДЗ по оси Y
    DataGetterOfY getter1alongy{ tmin1, tmax1, curve1 };
    DataGetterOfY getter2alongy{ tmin2, tmax2, curve2 };

    ///////////
    //       //
    // MIN Y //
    //       //
    ///////////
    auto [
      commonYmin,
      tofCommonYmin1,
      pointofCommonYmin1,
      tofCommonYmin2,
      pointofCommonYmin2
    ] = CommonRangeHelper::ofLowest(getter1alongy, getter2alongy);
    ///////////
    //       //
    // MAX Y //
    //       //
    ///////////
    auto [
      commonYmax,
      tofCommonYmax1,
      pointofCommonYmax1,
      tofCommonYmax2,
      pointofCommonYmax2
    ] = CommonRangeHelper::ofHighest(getter1alongy, getter2alongy);
    
    const auto commonRangeY = commonYmax - commonYmin;
    if (commonRangeX > commonRangeY)
    {
      using namespace geom2d::solver_normal_and_normal;

      solver<DataGetterOfX> theSolver{ tofCommonXmin1, tofCommonXmax1, curve1, tofCommonXmin2, tofCommonXmax2, curve2 };
      result = theSolver.execute();
    }
    else
    {
      using namespace geom2d::solver_normal_and_normal;

      solver<DataGetterOfY> theSolver{ tofCommonYmin1, tofCommonYmax1, curve1, tofCommonYmin2, tofCommonYmax2, curve2 };
      result = theSolver.execute();
    }
  }
  return result;
}

std::vector<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::execScreenAndScreen
  (
    const double tmin1,
    const double tmax1,
    const baseCurve& curve1,
    const double tmin2,
    const double tmax2,
    const baseCurve& curve2
  )
{
  std::vector<geom2d::IntersecctionSolutionType> result;
  // начинаем с нижних точек обеих кривых и идём к точке пересечения методом ступенек
  const auto P1 = curve1.getPoint(tmin1);
  const auto Q1 = curve1.getPoint(tmax1);
  StatOfCurvePiece scp1{ tmin1, P1, tmax1, Q1 };

  const auto P2 = curve2.getPoint(tmin2);
  const auto Q2 = curve2.getPoint(tmax2);
  StatOfCurvePiece scp2{ tmin2, P2, tmax2, Q2 };

  // сначала утилизируем тривиальные случаи
  if (scp1.pointOfxmin().x >= scp2.pointOfxmax().x)
  {
    // проверяем эти точки на совпадение
    const auto t1 = scp1.tOfxmin();
    const auto p1 = scp1.pointOfxmin();

    const auto t2 = scp2.tOfxmax();
    const auto p2 = scp2.pointOfxmax();

    if (point::isSame(p1, p2))
    {
      result.push_back(IntersecctionSolutionType{ 0.5 * (p1 + p2), t1, t2 });
    }
    return result;
  }
  else if (scp2.pointOfxmin().x >= scp1.pointOfxmax().x)
  {
    // проверяем эти точки на совпадение
    const auto t1 = scp1.tOfxmax();
    const auto p1 = scp1.pointOfxmax();

    const auto t2 = scp2.tOfxmin();
    const auto p2 = scp2.pointOfxmin();

    if (point::isSame(p1, p2))
    {
      result.push_back(IntersecctionSolutionType{ 0.5 * (p1 + p2), t1, t2 });
    }
    return result;
  }
  else if (scp1.pointOfymin().y >= scp2.pointOfymax().y)
  {
    const auto t1 = scp1.tOfymin();
    const auto p1 = scp1.pointOfymin();

    const auto t2 = scp2.tOfymax();
    const auto p2 = scp2.pointOfymax();

    if (point::isSame(p1, p2))
    {
      result.push_back(IntersecctionSolutionType{ 0.5 * (p1 + p2), t1, t2 });
    }
    return result;
  }
  else if (scp2.pointOfymin().y >= scp1.pointOfymax().y)
  {
    const auto t1 = scp1.tOfymax();
    const auto p1 = scp1.pointOfymax();

    const auto t2 = scp2.tOfymin();
    const auto p2 = scp2.pointOfymin();

    if (point::isSame(p1, p2))
    {
      result.push_back(IntersecctionSolutionType{ 0.5 * (p1 + p2), t1, t2 });
    }
    return result;
  }
  else
  {
    // если мы пришли сюда, значит есть пересечение ОДЗ двух кривых как по X, так и по Y
    // находим общее ОДЗ по оси X
    DataGetterOfX getter1alongx{ tmin1, tmax1, curve1 };
    DataGetterOfX getter2alongx{ tmin2, tmax2, curve2 };

    ///////////
    //       //
    // MIN X //
    //       //
    ///////////
    auto [
      commonXmin,
      tofCommonXmin1,
      pointofCommonXmin1,
      tofCommonXmin2,
      pointofCommonXmin2
    ] = CommonRangeHelper::ofLowest(getter1alongx, getter2alongx);
    ///////////
    //       //
    // MAX X //
    //       //
    ///////////
    auto [
      commonXmax,
      tofCommonXmax1,
      pointofCommonXmax1,
      tofCommonXmax2,
      pointofCommonXmax2
    ] = CommonRangeHelper::ofHighest(getter1alongx, getter2alongx);
    
    const auto commonRangeX = commonXmax - commonXmin;
    // находим общее ОДЗ по оси Y
    DataGetterOfY getter1alongy{ tmin1, tmax1, curve1 };
    DataGetterOfY getter2alongy{ tmin2, tmax2, curve2 };

    ///////////
    //       //
    // MIN Y //
    //       //
    ///////////
    auto [
      commonYmin,
      tofCommonYmin1,
      pointofCommonYmin1,
      tofCommonYmin2,
      pointofCommonYmin2
    ] = CommonRangeHelper::ofLowest(getter1alongy, getter2alongy);
    ///////////
    //       //
    // MAX Y //
    //       //
    ///////////
    auto [
      commonYmax,
      tofCommonYmax1,
      pointofCommonYmax1,
      tofCommonYmax2,
      pointofCommonYmax2
    ] = CommonRangeHelper::ofHighest(getter1alongy, getter2alongy);

    const auto commonRangeY = commonYmax - commonYmin;
    if (commonRangeX > commonRangeY)
    {
      using namespace geom2d::solver_screen_and_screen;

      solver<DataGetterOfX> theSolver{ tofCommonXmin1, tofCommonXmax1, curve1, tofCommonXmin2, tofCommonXmax2, curve2 };
      result = theSolver.execute();
    }
    else
    {
      using namespace geom2d::solver_screen_and_screen;

      solver<DataGetterOfY> theSolver{ tofCommonYmin1, tofCommonYmax1, curve1, tofCommonYmin2, tofCommonYmax2, curve2 };
      result = theSolver.execute();
    }
  }
  return result;
}

void geom2d::curveIntersector::dumpIntersections(std::ostream& ost) const
{
  const auto n = solutionPoints.size();
  if ((solutionParameterOnCurve1.size() != n) or (solutionParameterOnCurve2.size() != n))
  {
    throw std::logic_error("impossible situation in geom2d::curveIntersector::dumpIntersections");
  }

  ost << "intersection count: " << n << std::endl;
  for (size_t j = 0; j < n; ++j)
  {
    const auto pnt = solutionPoints[j];
    
    const auto t1 = solutionParameterOnCurve1[j];
    const auto t2 = solutionParameterOnCurve2[j];


    ost << "point: {" << pnt.x << ", " << pnt.y << "}" << " | ";
    ost << "t on curve No. 1: " << t1 << " | ";
    ost << "t on curve No. 2: " << t2 << std::endl;

    // verification
    const auto p1 = m_curve1.getPoint(t1);
    const auto p2 = m_curve2.getPoint(t2);

    if (not point::isSame(p1, p2))
    {
      throw std::logic_error("intersection point of two curves is incorrect!");
    }
  }
}

std::optional<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::exec_Any_and_Any_Unique_Intersection
  (
    const double tmin1,
    const double tmax1,
    const baseCurve& curve1,
    const double tmin2,
    const double tmax2,
    const baseCurve& curve2
  )
{
  const auto cc1 = baseCurve::getCurveClass(tmin1, tmax1, curve1);
  const auto cc2 = baseCurve::getCurveClass(tmin2, tmax2, curve2);

  if ((cc1 == curveClass::Point) or (cc2 == curveClass::Point))
  {
    return std::nullopt;
  }

  if (cc1 == cc2)
  {
    return std::nullopt;
  }
  // утилизация случаев, когда конец или начало одной кривой лежат на второй и наоборот
  {
    // Концы кривой Normal лежат ли на кривой Screen?
    std::initializer_list tof1{ tmin1, tmax1 };
    if (cc2 != curveClass::PlatoX)
    {
      // вдоль оси X
      for (const auto tof1 : std::initializer_list{ tmin1, tmax1 })
      {
        const auto pof1 = curve1.getPoint(tof1);

        using namespace geom2d::point_and_curve_alongaxis;

        solver<DataGetterOfX> solv{ pof1, tmin2, tmax2, curve2 };
        const auto result = solv.execute();
        if (not result) continue;

        const auto ton2 = result.value();
        const auto pon2 = curve2.getPoint(ton2);
        return IntersecctionSolutionType{ 0.5 * (pof1 + pon2), tof1, ton2 };
      }
    }
    if (cc2 != curveClass::PlatoY)
    {
      // вдоль оси Y
      for (const auto tof1 : std::initializer_list{ tmin1, tmax1 })
      {
        const auto pon1 = curve1.getPoint(tof1);

        using namespace geom2d::point_and_curve_alongaxis;

        solver<DataGetterOfY> solv{ pon1, tmin2, tmax2, curve2 };
        const auto result = solv.execute();
        if (not result) continue;

        const auto ton2 = result.value();
        const auto pon2 = curve2.getPoint(ton2);
        return IntersecctionSolutionType{ 0.5 * (pon1 + pon2), tof1, ton2 };
      }
    }
  }

  {
    // Концы кривой Screen лежат ли на кривой Normal?
    std::initializer_list tsofScreen{ tmin2, tmax2 };
    if (cc1 != curveClass::PlatoX)
    {
      // вдоль оси X
      for (const auto tof2 : tsofScreen)
      {
        const auto pof2 = curve2.getPoint(tof2);

        using namespace geom2d::point_and_curve_alongaxis;

        solver<DataGetterOfX> solv{ pof2, tmin1, tmax1, curve1 };
        const auto result = solv.execute();
        if (not result) continue;

        const auto ton1 = result.value();
        const auto pon1 = curve1.getPoint(ton1);
        return IntersecctionSolutionType{ 0.5 * (pof2 + pon1), ton1, tof2 };
      }
    }
    if (cc1 != curveClass::PlatoY)
    {
      // вдоль оси Y
      for (const auto tof2 : tsofScreen)
      {
        const auto pof2 = curve2.getPoint(tof2);

        using namespace geom2d::point_and_curve_alongaxis;

        solver<DataGetterOfY> solv{ pof2, tmin1, tmax1, curve1 };
        const auto result = solv.execute();
        if (not result) continue;

        const auto ton1 = result.value();
        const auto pon1 = curve1.getPoint(ton1);
        return IntersecctionSolutionType{ 0.5 * (pof2 + pon1), ton1, tof2 };
      }
    }
  }

  // далее проверяем полноценное пересечение или его отсутствие (граничные случаи утилизированы выше)
  return findUniqueIntersection(tmin1, tmax1, curve1, tmin2, tmax2, curve2);
}

//-----------------------------------------------------------------------------
