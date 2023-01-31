#include "curveIntersector.h"
#include "findFunctionRoots.h"
#include "vector.h"
#include "segmentCurve.h"
#include "math.h"
#include "Axis.h"
#include "StatOfCurvePiece.h"
#include "CommonRangeHelper.h"
#include "solver_unique_intersection.h"
#include "solver_point_and_alongaxis.h"
#include "solver_points_and_alongaxis.h"
#include "solver_point_and_point.h"
#include "solver_normal_and_normal_along_x.h"
#include "solver_normal_and_normal.h"

#include <vector>
#include <array>
#include <algorithm>

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

void geom2d::curveIntersector::perform()
{
  // Разбиваем кривую 1 на точки монотонности
  const auto tmanifold1 = rootsOfCurveVelocity(m_curve1);
  
  // Тоже для кривой 2
  const auto tmanifold2 = rootsOfCurveVelocity(m_curve2);

  // точки монотонности не включают начало и конец
  // для каждой пары значений интервала монотонности и каждой кривой ищем пересечение,
  // зная что кривые монотонны на выбранном участке
  // tmanifold - отсортированный массив t
  auto assembleVectorOft = [](const double tmin, const double tmax, const auto tmanifold)
  {
    std::vector<double> result;
    result.push_back(tmin);
    for (const auto t : tmanifold)
    {
      result.push_back(t);
    }
    result.push_back(tmax);
    return result;
  };
  const auto tman1 = assembleVectorOft(m_curve1.parameterMin(), m_curve1.parameterMax(), tmanifold1);
  const auto tman2 = assembleVectorOft(m_curve2.parameterMin(), m_curve2.parameterMax(), tmanifold2);
  for (size_t i = 0; i < tman1.size() - 1; ++i)
  {
    const double tmin1 = tman1[i];
    const double tmax1 = tman1[i + 1];
    for (size_t j = 0; j < tman2.size() - 1; ++j)
    {
      const double tmin2 = tman2[j];
      const double tmax2 = tman2[j + 1];
      perform(tmin1, tmax1, tmin2, tmax2);
    }
  }
}

//-----------------------------------------------------------------------------

void geom2d::curveIntersector::perform(
  const double tmin1,
  const double tmax1,
  const double tmin2,
  const double tmax2)
{
  const auto pBegOfCurve1 = m_curve1.getPoint(tmin1);
  const auto pEndOfCurve1 = m_curve1.getPoint(tmax1);

  const auto pBegOfCurve2 = m_curve1.getPoint(tmin2);
  const auto pEndOfCurve2 = m_curve1.getPoint(tmax2);

  // возвращаем класс кривой
  auto getCurveClass = [](const point pBegOfCurve, const point pEndOfCurve) -> curveClass
  {
    const auto [xbeg, ybeg] = pBegOfCurve;
    const auto [xend, yend] = pEndOfCurve;
    constexpr auto tol2d = math::tolerance::tolPoint;
    //constexpr auto halfTol2d = tol2d / 2.0;

    const bool xInc = (xend - xbeg) > tol2d;
    const bool xDec = (xbeg - xend) > tol2d;

    const bool yInc = (yend - ybeg) > tol2d;
    const bool yDec = (ybeg - yend) > tol2d;

    if (xInc)
    {
      if (yDec)
      {
        return curveClass::Screen;
      }
      else if (yInc)
      {
        return curveClass::Normal;
      }
      else
      {
        return curveClass::PlatoY;
      }
    }
    else if (xDec)
    {
      if (yDec)
      {
        return curveClass::Normal;
      }
      else if (yInc)
      {
        return curveClass::Screen;
      }
      else
      {
        return curveClass::PlatoY;
      }
    }
    else
    {
      if (yDec)
      {
        return curveClass::PlatoX;
      }
      else if (yInc)
      {
        return curveClass::PlatoX;
      }
      else
      {
        return curveClass::Point;
      }
    }
  };

  const auto classOfCurve1 = getCurveClass(pBegOfCurve1, pEndOfCurve1);
  const auto classOfCurve2 = getCurveClass(pBegOfCurve2, pEndOfCurve2);

  // решаем задачу для каждого возможного сочетания классов отдельно
  switch (classOfCurve1)
  {
  case curveClass::Screen:
    switch (classOfCurve2)
    {
    case curveClass::Screen:
    break;
    case curveClass::Normal:
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
    break;
    case curveClass::PlatoX:
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
    break;
    case curveClass::PlatoY:
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
    break;
    case curveClass::Point:
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
    break;
  }
  break;
  case curveClass::Normal:
    switch (classOfCurve2)
    {
      case curveClass::Screen:
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
      break;
      case curveClass::Normal:
      break;
      case curveClass::PlatoX:
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
      break;
      case curveClass::PlatoY:
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
      break;
      case curveClass::Point:
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
      break;
    }
    break;
  case curveClass::PlatoX:
    switch (classOfCurve2)
    {
      case curveClass::Screen:
      case curveClass::Normal:
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
      break;
      case curveClass::PlatoX:
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
      break;
      case curveClass::PlatoY:
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
      break;
      case curveClass::Point:
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
      break;
    }
    break;
  case curveClass::PlatoY:
    switch (classOfCurve2)
    {
      case curveClass::Screen:
      case curveClass::Normal:
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
      break;
      case curveClass::PlatoX:
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
      break;
      case curveClass::PlatoY:
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
      break;
      case curveClass::Point:
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
      break;
    }
    break;
  case curveClass::Point:
    switch (classOfCurve2)
    {
      case curveClass::Screen:
      case curveClass::Normal:
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
      break;
      case curveClass::PlatoX:
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
      break;
      case curveClass::PlatoY:
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
      break;
      case curveClass::Point:
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
      break;
    }
    break;
  }
}

//-----------------------------------------------------------------------------

inline std::set<double>
geom2d::curveIntersector::rootsOfCurveVelocity(const baseCurve & curve)
{
  const double tmin1 = curve.parameterMin();
  const double tmax1 = curve.parameterMax();

  // используем свой компаратор для исключения повторяющихся значений параметра
  auto curveParameterLess = [](const double lhs, const double rhs) -> bool
  {
    return (rhs - lhs) > math::tolerance::tolNumeric;
  };

  using parameterContainer = std::set<double, decltype(curveParameterLess)>;

  // I - ищем множество значений параметра, в которых кривая по x меняет знак монотонности
  // то есть, либо меняет уменьшение на увеличение, либо наоборот
  auto curve1u = [&curve](const double t) -> double
  {
    return curve.getVelocity(t).x;
  };
  const auto rootsOfCurve1u = math::findFunctionRoots(tmin1, tmax1, curve1u);

  // II - ищем множество значений параметра, в которых кривая по y меняет знак монотонности
  // то есть, либо меняет уменьшение на увеличение, либо наоборот
  auto curve1v = [&curve](const double t) -> double
  {
    return curve.getVelocity(t).y;
  };
  const auto rootsOfCurve1v = math::findFunctionRoots(tmin1, tmax1, curve1v);

  parameterContainer allRoots(curveParameterLess);
  // собираем корни по u
  for (const auto t : rootsOfCurve1u)
  {
    allRoots.insert(t);
  }
  // собираем корни по v
  for (const auto t : rootsOfCurve1v)
  {
    allRoots.insert(t);
  }
  // таким образом исключили все "дупликаты"
  // теперь собираем всё в обычный std::set
  std::set<double> result;
  for (const auto t : allRoots)
  {
    result.insert(t);
  }
  return result;
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
  // проверяем сначала концы кривой curvePlatoX
  {
    std::initializer_list<double> listoft{ tminOfPlatox, tmaxOfPlatox };
    for (const auto tofPlatoX : listoft)
    {
      using namespace geom2d::point_and_curve_alongaxis;
      const auto pointofPlatoX = curvePlatoX.getPoint(tofPlatoX);
      solver<DataGetterOfX> solv{ pointofPlatoX, tminOfPlatoy , tmaxOfPlatoy , curvePlatoY };
      const auto result = solv.execute();
      if (result)
      {
        const auto tofPlatoY = result.value();
        const auto pointofPlatoY = curvePlatoY.getPoint(tofPlatoY);
        return IntersecctionSolutionType{ 0.5 * (pointofPlatoY + pointofPlatoX), tofPlatoX, tofPlatoY };
      }
    }
  }

  // теперь проверяем концы кривой curvePlatoY
  {
    std::initializer_list<double> listoft{ tminOfPlatoy, tmaxOfPlatoy };
    for (const auto tofPlatoY : listoft)
    {
      using namespace geom2d::point_and_curve_alongaxis;
      const auto pointofPlatoY = curvePlatoY.getPoint(tofPlatoY);
      solver<DataGetterOfY> solv{ pointofPlatoY, tminOfPlatox , tmaxOfPlatox , curvePlatoX };
      const auto result = solv.execute();
      if (result)
      {
        const auto tofPlatoX = result.value();
        const auto pointofPlatoX = curvePlatoX.getPoint(tofPlatoX);
        return IntersecctionSolutionType{ 0.5 * (pointofPlatoY + pointofPlatoX), tofPlatoX, tofPlatoY };
      }
    }
  }
  // далее проверяем полноценное пересечение или его отсутствие (граничные случаи утилизированы выше)
  return findUniqueIntersection(tminOfPlatox, tmaxOfPlatox, curvePlatoX, tminOfPlatoy, tmaxOfPlatoy, curvePlatoY);
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

  
  // Все граничные случаи рассмотрены, значит кривые либо четко пересекаются, либо нет
  return findUniqueIntersection(tminOfPlatox, tmaxOfPlatox, curvePlatoX, tminOfAny, tmaxOfAny, curveAny);
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
  std::initializer_list listOfPlatoYt{ tminOfPlatoy , tmaxOfPlatoy };
  // Сначала тестируем конечные точки кривой PlatoY - лежат ли они на Any?
    // trace along X axis of curveAny
  for (const auto t : listOfPlatoYt)
  {
    const auto p = curvePlatoY.getPoint(t);
    using namespace geom2d::point_and_curve_alongaxis;

    solver<DataGetterOfX> solv{ p , tminOfAny , tmaxOfAny , curveAny };
    const auto result = solv.execute();
    if (not result) continue;

    const auto tofAny = result.value();
    const auto pofAny = curveAny.getPoint(tofAny);
    return IntersecctionSolutionType{ 0.5 * (p + pofAny), t, tofAny };
  }
  // trace along Y axis of curveAny
  for (const auto t : listOfPlatoYt)
  {
    const auto p = curvePlatoY.getPoint(t);
    using namespace geom2d::point_and_curve_alongaxis;

    solver<DataGetterOfY> solv{ p, tminOfAny , tmaxOfAny , curveAny };
    const auto result = solv.execute();
    if (not result) continue;

    const auto tofAny = result.value();
    const auto pofAny = curveAny.getPoint(tofAny);
    return IntersecctionSolutionType{ 0.5 * (p + pofAny), t, tofAny };
  }


  // А теперь тестируем конечные точки кривой Any - лежат ли они на PlatoY?
  // проходимся только вдоль оси X так как у кривой PlatoY нет диапазона по Y
  std::initializer_list listOfAnyt{ tminOfAny , tmaxOfAny };
  for (const auto t : listOfAnyt)
  {
    const auto p = curveAny.getPoint(t);
    using namespace geom2d::point_and_curve_alongaxis;

    solver<DataGetterOfX> solv{ p, tminOfPlatoy, tmaxOfPlatoy, curvePlatoY };
    const auto result = solv.execute();
    if (not result) continue;

    const auto tofPlatoy = result.value();
    const auto pofPlatoy = curvePlatoY.getPoint(tofPlatoy);
    return IntersecctionSolutionType{ 0.5 * (p + pofPlatoy), tofPlatoy, t };
  }


  // Все граничные случаи рассмотрены, значит кривые либо четко пересекаются, либо отстоят друг от друга на достаточном
  // расстоянии, чтобы мы могли считать их непересекающимися
  return findUniqueIntersection(tminOfPlatoy, tmaxOfPlatoy, curvePlatoY, tminOfAny, tmaxOfAny, curveAny);
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
  // утилизация случаев, когда конец или начало одной кривой лежат на второй и наоборот
  {
    // Концы кривой Normal лежат ли на кривой Screen?
    std::initializer_list tsofNormal{ tminOfNormal, tmaxOfNormal };
    // вдоль оси X
    for (const auto tofNormal : tsofNormal)
    {
      const auto pofNormal = curveNormal.getPoint(tofNormal);

      using namespace geom2d::point_and_curve_alongaxis;

      solver<DataGetterOfX> solv{ pofNormal, tminOfScreen, tmaxOfScreen, curveScreen };
      const auto result = solv.execute();
      if (not result) continue;

      const auto tofScreen = result.value();
      const auto pofScreen = curveScreen.getPoint(tofScreen);
      return IntersecctionSolutionType{ 0.5 * (pofNormal + pofScreen), tofNormal, tofScreen };
    }
    // вдоль оси Y
    for (const auto tofNormal : tsofNormal)
    {
      const auto pofNormal = curveNormal.getPoint(tofNormal);

      using namespace geom2d::point_and_curve_alongaxis;

      solver<DataGetterOfY> solv{ pofNormal, tminOfScreen, tmaxOfScreen, curveScreen };
      const auto result = solv.execute();
      if (not result) continue;

      const auto tofScreen = result.value();
      const auto pofScreen = curveScreen.getPoint(tofScreen);
      return IntersecctionSolutionType{ 0.5 * (pofNormal + pofScreen), tofNormal, tofScreen };
    }
  }

  {
    // Концы кривой Screen лежат ли на кривой Normal?
    std::initializer_list tsofScreen{ tminOfScreen, tmaxOfScreen };
    // вдоль оси X
    for (const auto tofScreen : tsofScreen)
    {
      const auto pofScreen = curveScreen.getPoint(tofScreen);

      using namespace geom2d::point_and_curve_alongaxis;

      solver<DataGetterOfX> solv{ pofScreen, tminOfNormal, tmaxOfNormal, curveNormal };
      const auto result = solv.execute();
      if (not result) continue;

      const auto tofNormal = result.value();
      const auto pofNormal = curveNormal.getPoint(tofNormal);
      return IntersecctionSolutionType{ 0.5 * (pofScreen + pofNormal), tofNormal, tofScreen };
    }
    // вдоль оси Y
    for (const auto tofScreen : tsofScreen)
    {
      const auto pofScreen = curveScreen.getPoint(tofScreen);

      using namespace geom2d::point_and_curve_alongaxis;

      solver<DataGetterOfY> solv{ pofScreen, tminOfNormal, tmaxOfNormal, curveNormal };
      const auto result = solv.execute();
      if (not result) continue;

      const auto tofNormal = result.value();
      const auto pofNormal = curveNormal.getPoint(tofNormal);
      return IntersecctionSolutionType{ 0.5 * (pofScreen + pofNormal), tofNormal, tofScreen };
    }
  }
  return std::nullopt;
}

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

      solver<DataGetterOfY> theSolver{ tofCommonXmin1, tofCommonXmax1, curve1, tofCommonXmin2, tofCommonXmax2, curve2 };
      result = theSolver.execute();
    }
  }
  return result;
}

std::vector<geom2d::IntersecctionSolutionType>
  geom2d::curveIntersector::traceTwoNormalsThroughXaxis
  (
    const double tofmin1,
    const double tofmax1,
    const baseCurve& curve1,
    const double tofmin2,
    const double tofmax2,
    const baseCurve& curve2
  )
{
  // в данной функции ОДЗ по X для двух кривых совпадает
  const auto pointOfLess1 = curve1.getPoint(tofmin1);
  const auto pointOfMore1 = curve1.getPoint(tofmax1);

  const auto pointOfLess2 = curve2.getPoint(tofmin2);
  const auto pointOfMore2 = curve2.getPoint(tofmax2);

  auto p1 = pointOfLess1;
  auto p2 = pointOfLess2;

  auto t1 = tofmin1;
  auto t2 = tofmin2;

  while (true)
  {
    if (point::isSame(p1, p2))
    {
      // здесь мы добавляем решение в контейнер и смещаемся вправо на столько, чтобы точки не совпадали
    }
    else if (p2.y > p1.y)
    {

    }
    else // p1.y > p2.y
    {

    }
  }

}

//-----------------------------------------------------------------------------

