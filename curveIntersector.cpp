#include "curveIntersector.h"
#include "findFunctionRoots.h"
#include "vector.h"
#include "segmentCurve.h"
#include "math.h"
#include "Axis.h"
#include "StatOfCurvePiece.h"
#include "CommonRangeHelper.h"
#include "solver_unique_intersection.h"
#include "solver_point_and_platox.h"
#include "solver_point_and_platoy.h"
#include "solver_point_and_point.h"

#include <vector>
#include <array>
#include <algorithm>

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
        break;
        case curveClass::PlatoX:
        break;
        case curveClass::PlatoY:
        break;
        case curveClass::Point:
        {
          const auto result =
            execPointAndAny(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
          if (result)
          {
            const auto [intPoint, tOnPoint, tOnAny] =
              result.value();
            solutionPoints.push_back(intPoint);
            solutionParameterOnCurve1.push_back(tOnAny);
            solutionParameterOnCurve2.push_back(tOnPoint);
          }
        }
        break;
      }
    break;
  case curveClass::Normal:
    switch (classOfCurve2)
    {
      case curveClass::Screen:
        break;
      case curveClass::Normal:
        break;
      case curveClass::PlatoX:
        break;
      case curveClass::PlatoY:
        break;
      case curveClass::Point:
      {
        const auto result =
          execPointAndAny(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
        if (result)
        {
          const auto [intPoint, tOnPoint, tOnAny] =
            result.value();
          solutionPoints.push_back(intPoint);
          solutionParameterOnCurve1.push_back(tOnAny);
          solutionParameterOnCurve2.push_back(tOnPoint);
        }
      }
      break;
    }
    break;
  case curveClass::PlatoX:
    switch (classOfCurve2)
    {
      case curveClass::Screen:
        break;
      case curveClass::Normal:
        break;
      case curveClass::PlatoX:
      {
        // обе кривые - PlatoX
        const auto result =
          execPlatoXAndPlatoX(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
        if (result)
        {
          const auto [intPoint, tOnPlatoX1, tOnPlatoX2] =
            result.value();
          solutionPoints.push_back(intPoint);
          solutionParameterOnCurve1.push_back(tOnPlatoX1);
          solutionParameterOnCurve2.push_back(tOnPlatoX2);
        }
      }
      break;
      case curveClass::PlatoY:
        break;
      case curveClass::Point:
      {
        // первая кривая - PlatoX
        // вторая кривая - Point
        const auto result =
          execPointAndPlatoX(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
        if (result)
        {
          const auto [intPoint, tOnPoint, tOnPlatoX] =
            result.value();
          solutionPoints.push_back(intPoint);
          solutionParameterOnCurve1.push_back(tOnPlatoX);
          solutionParameterOnCurve2.push_back(tOnPoint);
        }
      }
      break;
    }
    break;
  case curveClass::PlatoY:
    switch (classOfCurve2)
    {
      case curveClass::Screen:
        break;
      case curveClass::Normal:
        break;
      case curveClass::PlatoX:
      {
        // вторая кривая - PlatoX
        // первая кривая - PlatoY
        const auto result =
          execPlatoXAndPlatoY(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
        if (result)
        {
          const auto [intPoint, tonCurve2, tonCurve1] = result.value();
          solutionPoints.push_back(intPoint);
          solutionParameterOnCurve1.push_back(tonCurve1);
          solutionParameterOnCurve2.push_back(tonCurve2);
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
          const auto [intPoint, tOnPlatoY1, tOnPlatoY2] =
            result.value();
          solutionPoints.push_back(intPoint);
          solutionParameterOnCurve1.push_back(tOnPlatoY1);
          solutionParameterOnCurve2.push_back(tOnPlatoY2);
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
          const auto [intPoint, tOnPoint, tOnPlatoY] =
            result.value();
          solutionPoints.push_back(intPoint);
          solutionParameterOnCurve1.push_back(tOnPlatoY);
          solutionParameterOnCurve2.push_back(tOnPoint);
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
          const auto [intPoint, tOnPoint, tOnAny] =
            result.value();
          solutionPoints.push_back(intPoint);
          solutionParameterOnCurve1.push_back(tOnPoint);
          solutionParameterOnCurve2.push_back(tOnAny);
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
          const auto [intPoint, tOnPoint, tOnPlatoX] =
            result.value();
          solutionPoints.push_back(intPoint);
          solutionParameterOnCurve1.push_back(tOnPoint);
          solutionParameterOnCurve2.push_back(tOnPlatoX);
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
          const auto [intPoint, tOnPoint, tOnPlatoY] =
            result.value();
          solutionPoints.push_back(intPoint);
          solutionParameterOnCurve1.push_back(tOnPoint);
          solutionParameterOnCurve2.push_back(tOnPlatoY);
        }
      }
      break;
      case curveClass::Point:
      {
        const auto result =
          execPointAndPoint(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
        if (result)
        {
          const auto [intPoint, tOnPoint1, tOnPoint2] =
            result.value();
          solutionPoints.push_back(intPoint);
          solutionParameterOnCurve1.push_back(tOnPoint1);
          solutionParameterOnCurve2.push_back(tOnPoint2);
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

std::optional<double> geom2d::curveIntersector::performPointVSAnyAlongY(
  const point pnt,
  const double tmin,
  const double tmax,
  const baseCurve& curve)
{
  std::optional<double> result;
  const auto firstPoint = curve.getPoint(tmin);
  const auto secondPoint = curve.getPoint(tmax);

  const auto minPoint = (firstPoint.y < secondPoint.y) ? firstPoint : secondPoint;
  const auto maxPoint = (firstPoint.y < secondPoint.y) ? secondPoint : firstPoint;

  const auto tOfMinPoint = (firstPoint.y < secondPoint.y) ? tmin : tmax;
  const auto tOfMaxPoint = (firstPoint.y < secondPoint.y) ? tmax : tmin;

  if (pnt.y <= minPoint.y)
  {
    if (point::isSame(pnt, minPoint))
    {
      return tOfMinPoint;
    }
    else
    {
      return std::nullopt;
    }
  }
  else if (pnt.y >= maxPoint.y)
  {
    if (point::isSame(pnt, maxPoint))
    {
      return tOfMaxPoint;
    }
    else
    {
      return std::nullopt;
    }
  }
  else
  {
    auto func = [&curve, yo = pnt.y](double t) -> double
    {
      return curve.getPoint(t).y - yo;
    };
    const double tRoot = math::findUniqueFunctionRoot(tmin, tmax, func);
    return tRoot;
  }
  return std::nullopt;
}

//-----------------------------------------------------------------------------

std::optional<double> geom2d::curveIntersector::performPointVSAnyAlongX(
  const point pnt,
  const double tmin,
  const double tmax,
  const baseCurve& curve)
{
  std::optional<double> result;
  const auto firstPoint = curve.getPoint(tmin);
  const auto secondPoint = curve.getPoint(tmax);

  const auto minPoint = (firstPoint.x < secondPoint.x) ? firstPoint : secondPoint;
  const auto maxPoint = (firstPoint.x < secondPoint.x) ? secondPoint : firstPoint;

  const auto tOfMinPoint = (firstPoint.x < secondPoint.x) ? tmin : tmax;
  const auto tOfMaxPoint = (firstPoint.x < secondPoint.x) ? tmax : tmin;

  if (pnt.x <= minPoint.x)
  {
    if (point::isSame(pnt, minPoint))
    {
      return tOfMinPoint;
    }
    else
    {
      return std::nullopt;
    }
  }
  else if (pnt.x >= maxPoint.x)
  {
    if (point::isSame(pnt, maxPoint))
    {
      return tOfMaxPoint;
    }
    else
    {
      return std::nullopt;
    }
  }
  else
  {
    auto func = [&curve, xo = pnt.x](double t) -> double
    {
      return curve.getPoint(t).x - xo;
    };
    const double tRoot = math::findUniqueFunctionRoot(tmin, tmax, func);
    return tRoot;
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
  // кривая 1 - точка
  // кривая 2 - плато 'X = const'
  
  bool intersectionFound{ false };
  point intersectionPoint;
  double tOfIntersectionOnCurvePoint{ 0.0 };
  double tOfIntersectionOnCurvePlato{ 0.0 };
  double currentDistance = math::infinite::distance;

  auto handlePointAndPlatoX =
    [
      &curvePoint,
      &curvePlatoX,
      &intersectionFound,
      &intersectionPoint,
      &tOfIntersectionOnCurvePoint,
      &tOfIntersectionOnCurvePlato,
      &currentDistance,
      tminOfPlatoX,
      tmaxOfPlatoX
    ](const double tOfPoint)
  {
    const point pnt = curvePoint.getPoint(tOfPoint);
    const auto tOfPlatoX = performPointVSAnyAlongY(pnt, tminOfPlatoX, tmaxOfPlatoX, curvePlatoX);
    if (tOfPlatoX)
    {
      const auto pointOnPlatoX = curvePlatoX.getPoint(tOfPlatoX.value());
      const auto dist = point::distance(pnt, pointOnPlatoX);

      if (dist < currentDistance)
      {
        intersectionFound = true;
        currentDistance = dist;
        intersectionPoint = 0.5 * (pnt + pointOnPlatoX);
        tOfIntersectionOnCurvePoint = tOfPoint;
        tOfIntersectionOnCurvePlato = tOfPlatoX.value();
      }
    }
  };

  handlePointAndPlatoX(tminOfPoint);
  handlePointAndPlatoX(tmaxOfPoint);

  if (intersectionFound)
  {
    return
      std::tuple
        {
          intersectionPoint,
          tOfIntersectionOnCurvePoint,
          tOfIntersectionOnCurvePlato
        };
  }
  else
  {
    return std::nullopt;
  }
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
  // кривая 1 - точка
  // кривая 2 - плато 'Y = const'

  bool intersectionFound{ false };
  point intersectionPoint;
  double tOfIntersectionOnCurvePoint{ 0.0 };
  double tOfIntersectionOnCurvePlato{ 0.0 };
  double currentDistance = math::infinite::distance;

  auto handlePointAndPlatoY =
    [
      &curvePoint,
      &curvePlatoY,
      &intersectionFound,
      &intersectionPoint,
      &tOfIntersectionOnCurvePoint,
      &tOfIntersectionOnCurvePlato,
      &currentDistance,
      tminOfPlatoY,
      tmaxOfPlatoY
    ](const double tOfPoint)
  {
    const point pnt = curvePoint.getPoint(tOfPoint);
    const auto tOfPlatoY =
      performPointVSAnyAlongX(pnt, tminOfPlatoY, tmaxOfPlatoY, curvePlatoY);
    if (tOfPlatoY)
    {
      const auto pointOnPlatoY = curvePlatoY.getPoint(tOfPlatoY.value());
      const auto dist = point::distance(pnt, pointOnPlatoY);

      if (dist < currentDistance)
      {
        intersectionFound = true;
        currentDistance = dist;
        intersectionPoint = 0.5 * (pnt + pointOnPlatoY);
        tOfIntersectionOnCurvePoint = tOfPoint;
        tOfIntersectionOnCurvePlato = tOfPlatoY.value();
      }
    }
  };

  handlePointAndPlatoY(tminOfPoint);
  handlePointAndPlatoY(tmaxOfPoint);

  if (intersectionFound)
  {
    return
      std::tuple
        {
          intersectionPoint,
          tOfIntersectionOnCurvePoint,
          tOfIntersectionOnCurvePlato
        };
  }
  else
  {
    return std::nullopt;
  }
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
  // кривая 1 - точка
  // кривая 2 - не платно и не точка

  bool intersectionFound{ false };
  point intersectionPoint;
  double tOfIntersectionOnPoint{ 0.0 };
  double tOfIntersectionOnAny{ 0.0 };
  double currentDistance = math::infinite::distance;

  // исследуем вторую кривую вдоль X или Y
  auto handleParsingAlongXorY = 
  [
    &curvePoint,
    &curveAny,
    tminOfAny,
    tmaxOfAny,
    &currentDistance,
    &intersectionFound,
    &intersectionPoint,
    &tOfIntersectionOnPoint,
    &tOfIntersectionOnAny
  ] (const double tOfAny, const parseAxis parseThrough)
  {
    const auto pnt = curvePoint.getPoint(tOfAny);
    const auto tOfPnt = (parseThrough == parseAxis::X) ?
      performPointVSAnyAlongX(pnt, tminOfAny, tmaxOfAny, curveAny) :
      performPointVSAnyAlongY(pnt, tminOfAny, tmaxOfAny, curveAny);
    if (tOfPnt)
    {
      const auto pntOnCurve = curveAny.getPoint(tOfPnt.value());
      const auto dist = point::distance(pnt, pntOnCurve);

      if (dist < currentDistance)
      {
        currentDistance = dist;
        intersectionFound = true;
        intersectionPoint = 0.5 * (pnt + pntOnCurve);
        tOfIntersectionOnPoint = tOfAny;
        tOfIntersectionOnAny = tOfPnt.value();
      }
    }
  };

  handleParsingAlongXorY(tminOfPoint, parseAxis::X);
  handleParsingAlongXorY(tmaxOfPoint, parseAxis::X);

  handleParsingAlongXorY(tminOfPoint, parseAxis::Y);
  handleParsingAlongXorY(tmaxOfPoint, parseAxis::Y);

  if (intersectionFound)
  {
    return
      std::tuple
      {
        intersectionPoint,
        tOfIntersectionOnPoint,
        tOfIntersectionOnAny
      };
  }
  else
  {
    return std::nullopt;
  }
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
  // про кривую №1
  const auto P1 = curvePlatoY1.getPoint(tmin1);
  const auto Q1 = curvePlatoY1.getPoint(tmax1);

  const auto less1 = (P1.x < Q1.x) ? P1 : Q1;
  const auto more1 = (P1.x < Q1.x) ? Q1 : P1;

  const auto tofLess1 = (P1.x < Q1.x) ? tmin1 : tmax1;
  const auto tofMore1 = (P1.x < Q1.x) ? tmax1 : tmin1;

  aabb aabbOfPlatoY1(std::initializer_list{ P1, Q1 });

  // про кривую №2
  const auto P2 = curvePlatoY2.getPoint(tmin2);
  const auto Q2 = curvePlatoY2.getPoint(tmax2);

  const auto less2 = (P2.x < Q2.x) ? P2 : Q2;
  const auto more2 = (P2.x < Q2.x) ? Q2 : P2;

  const auto tofLess2 = (P2.x < Q2.x) ? tmin2 : tmax2;
  const auto tofMore2 = (P2.x < Q2.x) ? tmax2 : tmin2;

  aabb aabbOfPlatoY2(std::initializer_list{ P2, Q2 });
  
  const bool doNotIntersected =
    ((aabbOfPlatoY1.ymin() - aabbOfPlatoY2.ymax()) > math::tolerance::tolPoint)
      or
    ((aabbOfPlatoY2.ymin() - aabbOfPlatoY1.ymax()) > math::tolerance::tolPoint);

  if (doNotIntersected)
  {
    return std::nullopt;
  }

  if (more1.x <= less2.x)
  {
    if (point::isSame(more1, less2))
    {
      return std::tuple{0.5 * (more1 + less2), tofMore1, tofLess2};
    }
    else
    {
      return std::nullopt;
    }
  }
  else if (more2.x <= less1.x)
  {
    if (point::isSame(more2, less1))
    {
      return std::tuple{0.5 * (more2 + less1), tofMore2, tofLess1};
    }
    else
    {
      return std::nullopt;
    }
  }
  else
  {
    // есть пересечение интервалов
    double commonXmin = 0.0;
    double commonXmax = 0.0;
    double t1OfCommonXmin = 0.0;
    double t2OfCommonXmin = 0.0;
    double t1OfCommonXmax = 0.0;
    double t2OfCommonXmax = 0.0;
    // рассматриваем начало общего интервала
    if (less1.x < less2.x)
    {
      commonXmin = less2.x;
      t2OfCommonXmin = tofLess2;
      auto func = [&curvePlatoY1, xo = commonXmin](const double t) -> double
      {
        return curvePlatoY1.getPoint(t).x - xo;
      };
      t1OfCommonXmin = math::findUniqueFunctionRoot(tmin1, tmax1, func);

    }
    else if (less2.x < less1.x)
    {
      commonXmin = less1.x;
      t1OfCommonXmin = tofLess1;
      auto func = [&curvePlatoY2, xo = commonXmin](const double t) -> double
      {
        return curvePlatoY2.getPoint(t).x - xo;
      };
      t2OfCommonXmin = math::findUniqueFunctionRoot(tmin2, tmax2, func);
    }
    else
    {
      commonXmin = less1.x;
      t1OfCommonXmin = tofLess1;
      t2OfCommonXmin = tofLess2;
    }
    
    // рассматриваем конец общего интервала
    if (more1.x < more2.x)
    {
      commonXmax = more1.x;
      t1OfCommonXmax = tofMore1;
      auto func = [&curvePlatoY2, xo = commonXmax](const double t) -> double
      {
        return curvePlatoY2.getPoint(t).x - xo;
      };
      t2OfCommonXmax = math::findUniqueFunctionRoot(tmin2, tmax2, func);
    }
    else if (more2.x < more1.x)
    {
      commonXmax = more2.x;
      t2OfCommonXmax = tofMore2;
      auto func = [&curvePlatoY1, xo = commonXmax](const double t) -> double
      {
        return curvePlatoY1.getPoint(t).x - xo;
      };
      t1OfCommonXmax = math::findUniqueFunctionRoot(tmin1, tmax1, func);
    }
    else
    {
      commonXmax = more1.x;
      t1OfCommonXmax = tofMore1;
      t2OfCommonXmax = tofMore2;
    }
    // теперь проверяем три точки:
    // - начало
    // - конец
    // общего отразка
    // Если хотя бы на в одной точке кривые совпадают, возвращаем точку пересечения
    // Середину рассматривать не буду - лень
    const auto pointOn1OfLower = curvePlatoY1.getPoint(t1OfCommonXmin);
    const auto pointOn2OfLower = curvePlatoY2.getPoint(t2OfCommonXmin);
    if (point::isSame(pointOn1OfLower, pointOn2OfLower))
    {
      return std::tuple{ 0.5 * (pointOn1OfLower + pointOn2OfLower), t1OfCommonXmin, t2OfCommonXmin };
    }
    const auto pointOn1OfUpper = curvePlatoY1.getPoint(t1OfCommonXmax);
    const auto pointOn2OfUpper = curvePlatoY2.getPoint(t2OfCommonXmax);
    if (point::isSame(pointOn1OfUpper, pointOn2OfUpper))
    {
      return std::tuple{ 0.5 * (pointOn1OfUpper + pointOn2OfUpper), t1OfCommonXmax, t2OfCommonXmax };
    }
  }
  return std::nullopt;
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
  // про кривую №1
  const auto P1 = curvePlatoX1.getPoint(tmin1);
  const auto Q1 = curvePlatoX1.getPoint(tmax1);

  const auto less1 = (P1.y < Q1.y) ? P1 : Q1;
  const auto more1 = (P1.y < Q1.y) ? Q1 : P1;

  const auto tofLess1 = (P1.y < Q1.y) ? tmin1 : tmax1;
  const auto tofMore1 = (P1.y < Q1.y) ? tmax1 : tmin1;

  aabb aabbOfPlatoY1(std::initializer_list{ P1, Q1 });

  // про кривую №2
  const auto P2 = curvePlatoX2.getPoint(tmin2);
  const auto Q2 = curvePlatoX2.getPoint(tmax2);

  const auto less2 = (P2.y < Q2.y) ? P2 : Q2;
  const auto more2 = (P2.y < Q2.y) ? Q2 : P2;

  const auto tofLess2 = (P2.y < Q2.y) ? tmin2 : tmax2;
  const auto tofMore2 = (P2.y < Q2.y) ? tmax2 : tmin2;

  aabb aabbOfPlatoY2(std::initializer_list{ P2, Q2 });

  const bool doNotIntersected =
    ((aabbOfPlatoY1.xmin() - aabbOfPlatoY2.xmax()) > math::tolerance::tolPoint)
    or
    ((aabbOfPlatoY2.xmin() - aabbOfPlatoY1.xmax()) > math::tolerance::tolPoint);

  if (doNotIntersected)
  {
    return std::nullopt;
  }

  if (more1.y <= less2.y)
  {
    if (point::isSame(more1, less2))
    {
      return std::tuple{ 0.5 * (more1 + less2), tofMore1, tofLess2 };
    }
    else
    {
      return std::nullopt;
    }
  }
  else if (more2.y <= less1.y)
  {
    if (point::isSame(more2, less1))
    {
      return std::tuple{ 0.5 * (more2 + less1), tofMore2, tofLess1 };
    }
    else
    {
      return std::nullopt;
    }
  }
  else
  {
    // есть пересечение интервалов
    double commonYmin = 0.0;
    double commonYmax = 0.0;
    double t1OfCommonYmin = 0.0;
    double t2OfCommonYmin = 0.0;
    double t1OfCommonYmax = 0.0;
    double t2OfCommonYmax = 0.0;
    // рассматриваем начало общего интервала
    if (less1.y < less2.y)
    {
      commonYmin = less2.y;
      t2OfCommonYmin = tofLess2;
      auto func = [&curvePlatoX1, yo = commonYmin](const double t) -> double
      {
        return curvePlatoX1.getPoint(t).y - yo;
      };
      t1OfCommonYmin = math::findUniqueFunctionRoot(tmin1, tmax1, func);
    }
    else if (less2.y < less1.y)
    {
      commonYmin = less1.y;
      t1OfCommonYmin = tofLess1;
      auto func = [&curvePlatoX2, yo = commonYmin](const double t) -> double
      {
        return curvePlatoX2.getPoint(t).y - yo;
      };
      t2OfCommonYmin = math::findUniqueFunctionRoot(tmin2, tmax2, func);
    }
    else
    {
      commonYmin = less1.y;
      t1OfCommonYmin = tofLess1;
      t2OfCommonYmin = tofLess2;
    }

    // рассматриваем конец общего интервала
    if (more1.y < more2.y)
    {
      commonYmax = more1.y;
      t1OfCommonYmax = tofMore1;
      auto func = [&curvePlatoX2, yo = commonYmax](const double t) -> double
      {
        return curvePlatoX2.getPoint(t).y - yo;
      };
      t2OfCommonYmax = math::findUniqueFunctionRoot(tmin2, tmax2, func);
    }
    else if (more2.y < more1.y)
    {
      commonYmax = more2.y;
      t2OfCommonYmax = tofMore2;
      auto func = [&curvePlatoX1, yo = commonYmax](const double t) -> double
      {
        return curvePlatoX1.getPoint(t).y - yo;
      };
      t1OfCommonYmax = math::findUniqueFunctionRoot(tmin1, tmax1, func);
    }
    else
    {
      commonYmax = more1.y;
      t1OfCommonYmax = tofMore1;
      t2OfCommonYmax = tofMore2;
    }
    // теперь проверяем две точки:
    // - начало
    // - конец
    // общего отразка
    // Если хотя бы на в одной точке кривые совпадают, возвращаем точку пересечения
    // Середину рассматривать не буду - лень
    
    const auto pointOn1OfLower = curvePlatoX1.getPoint(t1OfCommonYmin);
    const auto pointOn2OfLower = curvePlatoX2.getPoint(t2OfCommonYmin);
    if (point::isSame(pointOn1OfLower, pointOn2OfLower))
    {
      return std::tuple{ 0.5 * (pointOn1OfLower + pointOn2OfLower), t1OfCommonYmin, t2OfCommonYmin };
    }

    const auto pointOn1OfUpper = curvePlatoX1.getPoint(t1OfCommonYmax);
    const auto pointOn2OfUpper = curvePlatoX2.getPoint(t2OfCommonYmax);
    if (point::isSame(pointOn1OfUpper, pointOn2OfUpper))
    {
      return std::tuple{ 0.5 * (pointOn1OfUpper + pointOn2OfUpper), t1OfCommonYmax, t2OfCommonYmax };
    }
  }
  return std::nullopt;
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
      using namespace geom2d::point_and_platoy;
      const auto pointofPlatoX = curvePlatoX.getPoint(tofPlatoX);
      solver solv{ pointofPlatoX, tminOfPlatoy , tmaxOfPlatoy , curvePlatoY };
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
      using namespace geom2d::point_and_platox;
      const auto pointofPlatoY = curvePlatoY.getPoint(tofPlatoY);
      solver solv{ pointofPlatoY, tminOfPlatox , tmaxOfPlatox , curvePlatoX };
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
