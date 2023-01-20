#include "curveIntersector.h"
#include "findFunctionRoots.h"
#include "vector.h"
#include "segmentCurve.h"
#include "math.h"
#include "Axis.h"
#include "StatOfCurvePiece.h"
#include "CommonRangeHelper.h"

#include <vector>
#include <array>
#include <algorithm>

//-----------------------------------------------------------------------------

void geom2d::curveIntersector::perform()
{
  // ��������� ������ 1 �� ����� ������������
  const auto tmanifold1 = rootsOfCurveVelocity(m_curve1);
  
  // ���� ��� ������ 2
  const auto tmanifold2 = rootsOfCurveVelocity(m_curve2);

  // ����� ������������ �� �������� ������ � �����
  // ��� ������ ���� �������� ��������� ������������ � ������ ������ ���� �����������,
  // ���� ��� ������ ��������� �� ��������� �������
  // tmanifold - ��������������� ������ t
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

  // ���������� ����� ������
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

  // ������ ������ ��� ������� ���������� ��������� ������� ��������
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
        // ��� ������ - PlatoX
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
        // ������ ������ - PlatoX
        // ������ ������ - Point
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
        break;
      case curveClass::PlatoY:
      {
        // ��� ������ - PlatoY
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
        // ������ ������ - PlatoY
        // ������ ������ - Point
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
        // ������ ������ - Point
        // ������ ������ - PlatoX
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
        // ������ ������ - Point
        // ������ ������ - PlatoY
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

  // ���������� ���� ���������� ��� ���������� ������������� �������� ���������
  auto curveParameterLess = [](const double lhs, const double rhs) -> bool
  {
    return (rhs - lhs) > math::tolerance::tolNumeric;
  };

  using parameterContainer = std::set<double, decltype(curveParameterLess)>;

  // I - ���� ��������� �������� ���������, � ������� ������ �� x ������ ���� ������������
  // �� ����, ���� ������ ���������� �� ����������, ���� ��������
  auto curve1u = [&curve](const double t) -> double
  {
    return curve.getVelocity(t).x;
  };
  const auto rootsOfCurve1u = math::findFunctionRoots(tmin1, tmax1, curve1u);

  // II - ���� ��������� �������� ���������, � ������� ������ �� y ������ ���� ������������
  // �� ����, ���� ������ ���������� �� ����������, ���� ��������
  auto curve1v = [&curve](const double t) -> double
  {
    return curve.getVelocity(t).y;
  };
  const auto rootsOfCurve1v = math::findFunctionRoots(tmin1, tmax1, curve1v);

  parameterContainer allRoots(curveParameterLess);
  // �������� ����� �� u
  for (const auto t : rootsOfCurve1u)
  {
    allRoots.insert(t);
  }
  // �������� ����� �� v
  for (const auto t : rootsOfCurve1v)
  {
    allRoots.insert(t);
  }
  // ����� ������� ��������� ��� "���������"
  // ������ �������� �� � ������� std::set
  std::set<double> result;
  for (const auto t : allRoots)
  {
    result.insert(t);
  }
  return result;
}

//-----------------------------------------------------------------------------

bool geom2d::curveIntersector::performPointXPoint(const point p, const point q)
{
  return point::isSame(p, q);
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
  bool intersectionFound{ false };
  point intersectionPoint;
  double tOnCurve1OfIntersection{ 0.0 };
  double tOnCurve2OfIntersection{ 0.0 };

  double currentDistance = 1.0e+10;
  auto handlePairOfPoints =
    [
      &curvePoint1,
      &curvePoint2,
      &intersectionFound,
      &intersectionPoint,
      &tOnCurve1OfIntersection,
      &tOnCurve2OfIntersection,
      &currentDistance
    ]
  (const double tOnCurve1, const double tOnCurve2)
  {
    const auto pnt1 = curvePoint1.getPoint(tOnCurve1);
    const auto pnt2 = curvePoint2.getPoint(tOnCurve2);

    const auto result = performPointXPoint(pnt1, pnt2);

    if (result)
    {
      const auto dist = point::distance(pnt1, pnt2);
      if (dist < currentDistance)
      {
        currentDistance = dist;
        intersectionFound = true;
        intersectionPoint = (pnt1 + pnt2) / 2.0;
        tOnCurve1OfIntersection = tOnCurve1;
        tOnCurve2OfIntersection = tOnCurve2;
      }
    }
  };
  
  handlePairOfPoints(tminOfPoint1, tminOfPoint2);
  handlePairOfPoints(tminOfPoint1, tmaxOfPoint2);
  handlePairOfPoints(tmaxOfPoint1, tminOfPoint2);
  handlePairOfPoints(tmaxOfPoint1, tmaxOfPoint2);

  if (intersectionFound)
  {
    return
      std::tuple
      {
        intersectionPoint,
        tOnCurve1OfIntersection,
        tOnCurve2OfIntersection
      };
  }
  else
  {
    return std::nullopt;
  }
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
  // ������ 1 - �����
  // ������ 2 - ����� 'X = const'
  
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

  handlePointAndPlatoX(tmaxOfPoint);
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
  // ������ 1 - �����
  // ������ 2 - ����� 'Y = const'

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
  // ������ 1 - �����
  // ������ 2 - �� ������ � �� �����

  bool intersectionFound{ false };
  point intersectionPoint;
  double tOfIntersectionOnPoint{ 0.0 };
  double tOfIntersectionOnAny{ 0.0 };
  double currentDistance = math::infinite::distance;

  // ��������� ������ ������ ����� X ��� Y
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
  // ��� ������ �1
  const auto P1 = curvePlatoY1.getPoint(tmin1);
  const auto Q1 = curvePlatoY1.getPoint(tmax1);

  const auto less1 = (P1.x < Q1.x) ? P1 : Q1;
  const auto more1 = (P1.x < Q1.x) ? Q1 : P1;

  const auto tofLess1 = (P1.x < Q1.x) ? tmin1 : tmax1;
  const auto tofMore1 = (P1.x < Q1.x) ? tmax1 : tmin1;

  aabb aabbOfPlatoY1(std::initializer_list{ P1, Q1 });

  // ��� ������ �2
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
    // ���� ����������� ����������
    double commonXmin = 0.0;
    double commonXmax = 0.0;
    double t1OfCommonXmin = 0.0;
    double t2OfCommonXmin = 0.0;
    double t1OfCommonXmax = 0.0;
    double t2OfCommonXmax = 0.0;
    // ������������� ������ ������ ���������
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
    
    // ������������� ����� ������ ���������
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
    // ������ ��������� ��� �����:
    // - ������
    // - �����
    // ������ �������
    // ���� ���� �� �� � ����� ����� ������ ���������, ���������� ����� �����������
    // �������� ������������� �� ���� - ����
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
  // ��� ������ �1
  const auto P1 = curvePlatoX1.getPoint(tmin1);
  const auto Q1 = curvePlatoX1.getPoint(tmax1);

  const auto less1 = (P1.y < Q1.y) ? P1 : Q1;
  const auto more1 = (P1.y < Q1.y) ? Q1 : P1;

  const auto tofLess1 = (P1.y < Q1.y) ? tmin1 : tmax1;
  const auto tofMore1 = (P1.y < Q1.y) ? tmax1 : tmin1;

  aabb aabbOfPlatoY1(std::initializer_list{ P1, Q1 });

  // ��� ������ �2
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
    // ���� ����������� ����������
    double commonYmin = 0.0;
    double commonYmax = 0.0;
    double t1OfCommonYmin = 0.0;
    double t2OfCommonYmin = 0.0;
    double t1OfCommonYmax = 0.0;
    double t2OfCommonYmax = 0.0;
    // ������������� ������ ������ ���������
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

    // ������������� ����� ������ ���������
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
    // ������ ��������� ��� �����:
    // - ������
    // - �����
    // ������ �������
    // ���� ���� �� �� � ����� ����� ������ ���������, ���������� ����� �����������
    // �������� ������������� �� ���� - ����
    
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
  /////////
  //     //
  //  *  // Platox curve
  //     //
  /////////
  // ������ ��������� ������ �����-������������ ��� OY (PlatoX)
  const auto PofPlatox = curvePlatoX.getPoint(tminOfPlatox);
  const auto QofPlatox = curvePlatoX.getPoint(tmaxOfPlatox);
  // ***
  const auto pointOfYmin = (PofPlatox.y < QofPlatox.y) ? PofPlatox : QofPlatox;
  const auto tOfYmin = (PofPlatox.y < QofPlatox.y) ? tminOfPlatox : tmaxOfPlatox;
  // ***
  const auto pointOfYmax = (PofPlatox.y < QofPlatox.y) ? QofPlatox : PofPlatox;
  const auto tOfYmax = (PofPlatox.y < QofPlatox.y) ? tmaxOfPlatox : tminOfPlatox;
  // ***
  const aabb pyAABB{ std::initializer_list{PofPlatox, QofPlatox} };
  /////////
  //     //
  //  *  // Platoy curve
  //     //
  /////////
  // ������ ��������� ������ �����-������������ ��� OX (PlatoY)
  const auto PofPlatoy = curvePlatoY.getPoint(tminOfPlatoy);
  const auto QofPlatoy = curvePlatoY.getPoint(tmaxOfPlatoy);
  // ***
  const auto pointOfXmin = (PofPlatoy.x < QofPlatoy.x) ? PofPlatoy : QofPlatoy;
  const auto tOfXmin = (PofPlatoy.x < QofPlatoy.x) ? tminOfPlatoy : tmaxOfPlatoy;
  // ***
  const auto pointOfXmax = (PofPlatoy.x < QofPlatoy.x) ? QofPlatoy : PofPlatoy;
  const auto tOfXmax = (PofPlatoy.x < QofPlatoy.x) ? tmaxOfPlatoy : tminOfPlatoy;
  // ***
  const aabb pxAABB{ std::initializer_list{PofPlatoy, QofPlatoy} };

  // ������������� ��������� ����� ������ PlatoX
  auto marginalIntersectorOfPointWithPlatoY =
    [
      &curvePlatoX,
      &curvePlatoY,
      pointOfXmin,
      tOfXmin,
      tminOfPlatoy,
      pointOfXmax,
      tOfXmax,
      tmaxOfPlatoy
    ]
  (const double tVertexPlatoX) -> std::optional<geom2d::IntersecctionSolutionType>
  {
    const auto vertexOfPlatoX = curvePlatoX.getPoint(tVertexPlatoX);
    if (vertexOfPlatoX.x <= pointOfXmin.x)
    {
      if (point::isSame(vertexOfPlatoX, pointOfXmin))
      {
        return std::tuple{ 0.5 * (vertexOfPlatoX + pointOfXmin), tVertexPlatoX , tOfXmin };
      }
      else
      {
        return std::nullopt;
      }
    }
    else if (vertexOfPlatoX.x >= pointOfXmax.x)
    {
      if (point::isSame(vertexOfPlatoX, pointOfXmax))
      {
        return std::tuple{ 0.5 * (vertexOfPlatoX + pointOfXmax), tVertexPlatoX , tOfXmax };
      }
      else
      {
        return std::nullopt;
      }
    }
    else
    {
      auto func = [&curvePlatoY, xo = vertexOfPlatoX.x](const double t) -> double
      {
        return curvePlatoY.getPoint(t).x - xo;
      };
      const auto tOnPlatoY = math::findUniqueFunctionRoot(tminOfPlatoy, tmaxOfPlatoy, func);
      const auto pointOnPlatoY = curvePlatoY.getPoint(tOnPlatoY);
      if (point::isSame(vertexOfPlatoX, pointOnPlatoY))
      {
        return std::tuple{ 0.5 * (vertexOfPlatoX + pointOnPlatoY), tVertexPlatoX, tOnPlatoY };
      }
      else
      {
        std::nullopt;
      }
    }
  };
  // ��������� ���������� ��������� ����� ������ PlatoX
  const auto result1 = marginalIntersectorOfPointWithPlatoY(tminOfPlatox);
  if (result1) return result1;

  const auto result2 = marginalIntersectorOfPointWithPlatoY(tmaxOfPlatox);
  if (result2) return result2;

  // ������������� ��������� ����� ������ PlatoY
  auto marginalIntersectorOfPointWithPlatoX =
    [
      &curvePlatoY,
      &curvePlatoX,
      pointOfYmin,
      tOfYmin,
      tminOfPlatox,
      pointOfYmax,
      tOfYmax,
      tmaxOfPlatox
    ]
  (const double tVertexPlatoY) -> std::optional<geom2d::IntersecctionSolutionType>
  {
    const auto vertexOfPlatoY = curvePlatoY.getPoint(tVertexPlatoY);
    if (vertexOfPlatoY.y <= pointOfYmin.y)
    {
      if (point::isSame(vertexOfPlatoY, pointOfYmin))
      {
        return std::tuple{ 0.5 * (vertexOfPlatoY + pointOfYmin) , tOfYmin, tVertexPlatoY };
      }
      else
      {
        return std::nullopt;
      }
    }
    else if (vertexOfPlatoY.y >= pointOfYmax.y)
    {
      if (point::isSame(vertexOfPlatoY, pointOfYmax))
      {
        return std::tuple{ 0.5 * (vertexOfPlatoY + pointOfYmax), tOfYmax, tVertexPlatoY };
    }
      else
      {
        return std::nullopt;
      }
  }
    else
    {
      auto func = [&curvePlatoX, yo = vertexOfPlatoY.y](const double t) -> double
      {
        return curvePlatoX.getPoint(t).y - yo;
      };
      const auto tOnPlatoX = math::findUniqueFunctionRoot(tminOfPlatox, tmaxOfPlatox, func);
      const auto pointOnPlatoX = curvePlatoX.getPoint(tOnPlatoX);
      if (point::isSame(vertexOfPlatoY, pointOnPlatoX))
      {
        return std::tuple{ 0.5 * (vertexOfPlatoY + pointOnPlatoX), tOnPlatoX, tVertexPlatoY };
      }
      else
      {
        std::nullopt;
      }
    }
  };
  // ��������� ���������� ��������� ����� ������ PlatoY
  const auto result3 = marginalIntersectorOfPointWithPlatoX(tminOfPlatoy);
  if (result3) return result3;

  const auto result4 = marginalIntersectorOfPointWithPlatoY(tmaxOfPlatoy);
  if (result4) return result4;

  // ����� ��������� ����������� �����������
  
  
#if 0
  const bool doNotIntersected =
    ((xAABB.xmin() - yAABB.xmax()) > math::tolerance::tolPoint)
      or
    ((yAABB.xmin() - xAABB.xmax()) > math::tolerance::tolPoint)
      or
    ((xAABB.ymin() - yAABB.ymax()) > math::tolerance::tolPoint)
      or
    ((yAABB.ymin() - xAABB.ymax()) > math::tolerance::tolPoint);

  if (doNotIntersected)
  {
    return std::nullopt;
  }

  // ��������� ������ PlatoX �� ��������� � ������ PlatoY
  XMutualPos mutPosByX {XMutualPos::inside};
  if      (yAABB.xmax() <= xAABB.xmin()) mutPosByX = XMutualPos::right;
  else if (yAABB.xmax() <= xAABB.xmax()) mutPosByX = XMutualPos::rightPartially;
  else if (yAABB.xmin() <  xAABB.xmin()) mutPosByX = XMutualPos::inside;
  else if (yAABB.xmin() <  xAABB.xmax()) mutPosByX = XMutualPos::leftPartially;
  else if (xAABB.xmax() <= yAABB.xmin()) mutPosByX = XMutualPos::left;
  else throw std::logic_error("impossible situ in intersection of PlatoX & PlatoY");

  // ��������� ������ PlatoY �� ��������� � ������ PlatoX


  if      ((yAABB.xmax() <= xAABB.xmin())) mutPosByX = XMutualPos::right;
  // *********************************************************************************************************************
  else if ((xAABB.xmin() <  yAABB.xmax())
            and
           (yAABB.xmax() <= xAABB.xmax())) mutPosByX = XMutualPos::rightPartially;
  // *********************************************************************************************************************
  else if ((yAABB.xmin() <  xAABB.xmin())
            and
           (xAABB.xmax() <  yAABB.xmax())) mutPosByX = XMutualPos::inside;
  // *********************************************************************************************************************
  else if ((yAABB.xmin() <  xAABB.xmax())
            and
           (xAABB.xmin() <= yAABB.xmin())) mutPosByX = XMutualPos::leftPartially;
  // *********************************************************************************************************************
  else if ((xAABB.xmax() <= yAABB.xmin())) mutPosByX = XMutualPos::left;
  // *********************************************************************************************************************
  else                                                   throw std::logic_error("impossible mutual position of PlatoX & PlatoY by x axis");

  // ��������� ������

#endif
  return std::nullopt;
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


  //const double tmin[] = { tmin1, tmin2 };
  //const double tmax[] = { tmax1, tmax2 };
  //const baseCurve* curves[] = { &curve1 , &curve2 };

  //const point p[2] =
  //{
  //  curves[0]->getPoint(tmin[0]),
  //  curves[1]->getPoint(tmin[1])
  //};
  //const point q[2] =
  //{
  //  curves[0]->getPoint(tmax[0]),
  //  curves[1]->getPoint(tmax[1])
  //};

  //const aabb aabb[2] =
  //{
  //  {std::initializer_list{p[0], q[0]}},
  //  {std::initializer_list{p[1], q[1]}}
  //};
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
    if (maxAxis == Axis::X)
    {
      return findUniqueIntersectionRefAlongX(tmin1, tmax1, curve1, tmin2, tmax2, curve2);
    }
    else
    {
      return findUniqueIntersectionRefAlongY(tmin1, tmax1, curve1, tmin2, tmax2, curve2);
    }
  }
  else
  {
    const auto result = (maxAxis == Axis::X) ?
      findUniqueIntersectionRefAlongX(tmin2, tmax2, curve2, tmin1, tmax1, curve1)
      :
      findUniqueIntersectionRefAlongX(tmin2, tmax2, curve2, tmin1, tmax1, curve1);
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

#if 0


  int theReferenceCurve = 1; // ������ � ���������� �������� ����� ����� �� ������������ ����
  double largestLength = 0.0;
  Axis theLargestAxis = Axis::X;
  if (aabb1.lengthx() > largestLength)
  {
    theReferenceCurve = 1;
    largestLength = aabb[0].lengthx();
    theLargestAxis = Axis::X;
  }
  if (aabb[0].lengthy() > largestLength)
  {
    theLargestCurve = 0;
    largestLength = aabb[0].lengthy();
    theLargestAxis = Axis::Y;
  }
  if (aabb[1].lengthx() > largestLength)
  {
    theLargestCurve = 1;
    largestLength = aabb[1].lengthx();
    theLargestAxis = Axis::X;
  }
  if (aabb[1].lengthy() > largestLength)
  {
    theLargestCurve = 1;
    largestLength = aabb[1].lengthy();
    theLargestAxis = Axis::Y;
  }
  const baseCurve& referenceCurve = *(curves[theLargestCurve]);
  
  const double trefmin = tmin[theLargestCurve];
  const point pref = referenceCurve.getPoint(trefmin);

  const double trefmax = tmax[theLargestCurve];
  const point qref = referenceCurve.getPoint(trefmax);
  
  StatOfCurvePiece scpRef{ trefmin , pref , trefmax , qref };

  const baseCurve& otherCurve = *(curves[1 - theLargestCurve]);

  const double tothmin = tmin[theLargestCurve];
  const point poth = otherCurve.getPoint(tothmin);
  
  const double tothmax = tmax[theLargestCurve];
  const point qoth = otherCurve.getPoint(tothmax);
  StatOfCurvePiece scpOth{ tothmin , poth , tothmax , qoth };

  if (theLargestAxis == Axis::X) return findUniqueIntersectionRefAlongX(trefmin, trefmax, referenceCurve, tothmin, tothmax, otherCurve);
  else                           return findUniqueIntersectionRefAlongX(trefmin, trefmax, referenceCurve, tothmin, tothmax, otherCurve);

#endif
#if 0
  if (theLargestAxis == Axis::X)
  {
    if (scpOth.pointOfxmax().x <= scpRef.pointOfxmin().x)
    {
      if (point::isSame(scpOth.pointOfxmax(), scpRef.pointOfxmin()))
      {
        if (theLargestCurve == 0)
        {
          return IntersecctionSolutionType{ 0.5 * (scpOth.pointOfxmax() + scpRef.pointOfxmin()), scpRef.tOfxmin(), scpOth.tOfxmax() };
        }
        else
        {
          return IntersecctionSolutionType{ 0.5 * (scpOth.pointOfxmax() + scpRef.pointOfxmin()), scpOth.tOfxmax(), scpRef.tOfxmin() };
        }
      }
      else
      {
        return std::nullopt;
      }
    }
    else if (scpOth.pointOfxmin().x >= scpRef.pointOfxmax().x)
    {
      if (point::isSame(scpOth.pointOfxmin(), scpRef.pointOfxmax()))
      {
        if (theLargestCurve == 0)
        {
          return IntersecctionSolutionType{ 0.5 * (scpOth.pointOfxmin() + scpRef.pointOfxmax()), scpRef.tOfxmax(), scpOth.tOfxmin() };
        }
        else
        {
          return IntersecctionSolutionType{ 0.5 * (scpOth.pointOfxmin() + scpRef.pointOfxmax()), scpOth.tOfxmin(), scpRef.tOfxmax() };
        }
      }
      else
      {
        return std::nullopt;
      }
    }
    else
    {
      // MIN X
      double commonXmin = 0.0;
      double trefOfCommonXmin = 0.0;
      point pointRefOfCommonXmin;
      double tothOfCommonXmin = 0.0;
      point pointOthOfCommonXmin;
      if (scpRef.pointOfxmin().x > scpOth.pointOfxmin().x)
      {
        commonXmin = scpRef.pointOfxmin().x;
        trefOfCommonXmin = scpRef.tOfxmin();
        auto func = [&otherCurve, xo = commonXmin](const double t) -> double
        {
          return otherCurve.getPoint(t).x - xo;
        };
        tothOfCommonXmin = math::findUniqueFunctionRoot(tothmin, tothmax, func);
        pointOthOfCommonXmin = otherCurve.getPoint(tothOfCommonXmin);
      }
      else if (scpRef.pointOfxmin().x < scpOth.pointOfxmin().x)
      {
        commonXmin = scpOth.pointOfxmin().x;
        tothOfCommonXmin = scpOth.tOfxmin();
        auto func = [&referenceCurve, xo = commonXmin](const double t) -> double
        {
          return referenceCurve.getPoint(t).x - xo;
        };
        trefOfCommonXmin = math::findUniqueFunctionRoot(trefmin, trefmax, func);
        pointRefOfCommonXmin = referenceCurve.getPoint(trefOfCommonXmin);
      }
      else
      {
        commonXmin = scpRef.pointOfxmin().x;
        trefOfCommonXmin = scpRef.tOfxmin();
        pointRefOfCommonXmin = scpRef.pointOfxmin();
        tothOfCommonXmin = scpOth.tOfxmin();
        pointOthOfCommonXmin = scpOth.pointOfxmin();
      }
      // MAX X
      double commonXmax = 0.0;
      double trefOfCommonXmax = 0.0;
      point pointRefOfCommonXmax;
      double tothOfCommonXmax = 0.0;
      point pointOthOfCommonXmax;
      if (scpRef.pointOfxmax().x < scpOth.pointOfxmax().x)
      {
        commonXmax = scpRef.pointOfxmax().x;
        trefOfCommonXmax = scpRef.tOfxmax();
        auto func = [&otherCurve, xo = commonXmax](const double t) -> double
        {
          return otherCurve.getPoint(t).x - xo;
        };
        tothOfCommonXmax = math::findUniqueFunctionRoot(tothmin, tothmax, func);
        pointOthOfCommonXmax = otherCurve.getPoint(tothOfCommonXmax);
      }
      else if (scpRef.pointOfxmax().x > scpOth.pointOfxmax().x)
      {
        commonXmax = scpOth.pointOfxmax().x;
        tothOfCommonXmax = scpOth.tOfxmax();
        auto func = [&referenceCurve, xo = commonXmax](const double t) -> double
        {
          return referenceCurve.getPoint(t).x - xo;
        };
        trefOfCommonXmax = math::findUniqueFunctionRoot(trefmin, trefmax, func);
        pointRefOfCommonXmax = referenceCurve.getPoint(trefOfCommonXmax);
      }
      else
      {
        commonXmax = scpRef.pointOfxmax().x;
        trefOfCommonXmax = scpRef.tOfxmax();
        pointRefOfCommonXmax = scpRef.pointOfxmax();
        tothOfCommonXmax = scpOth.tOfxmax();
        pointOthOfCommonXmax = scpOth.pointOfxmax();
      }
      // ��������� ����� ��� �� X
      //  (*) ��������� ����� �����        
      if (point::isSame(pointOthOfCommonXmin, pointRefOfCommonXmin))
      {
        if (theLargestCurve == 0)
        {
          return IntersecctionSolutionType{ 0.5 * (pointOthOfCommonXmin + pointRefOfCommonXmin), trefOfCommonXmin, tothOfCommonXmin };
        }
        else
        {
          return IntersecctionSolutionType{ 0.5 * (pointOthOfCommonXmin + pointRefOfCommonXmin), tothOfCommonXmin, trefOfCommonXmin };
        }
      }
      //  (*) ��������� ������ �����
      if (point::isSame(pointOthOfCommonXmax, pointRefOfCommonXmax))
      {
        if (theLargestCurve == 0)
        {
          return IntersecctionSolutionType{ 0.5 * (pointOthOfCommonXmax + pointRefOfCommonXmax), trefOfCommonXmax, tothOfCommonXmax };
        }
        else
        {
          return IntersecctionSolutionType{ 0.5 * (pointOthOfCommonXmax + pointRefOfCommonXmax), tothOfCommonXmax, trefOfCommonXmax };
        }
      }
      // ����� �� ���� ������ ��������� ���� ��� ���� ������
      double fleft  = pointOthOfCommonXmin.y - pointRefOfCommonXmin.y;
      double fright = pointOthOfCommonXmax.y - pointRefOfCommonXmax.y;
      const bool othIsAboveRef = (fleft > 0.0) and (fright > 0.0);
      const bool refIsUnderRef = (fleft < 0.0) and (fright < 0.0);
      if (othIsAboveRef or refIsUnderRef)
      {
        return std::nullopt;
      }
      while (true)
      {
        double tothMiddle = 0.5 * (tothOfCommonXmin + tothOfCommonXmax);
        point pointOthOnMiddle = otherCurve.getPoint(tothMiddle);
        
        double trefMiddle = 0.0;
        point pointRefOnMiddle;
        if (pointOthOnMiddle.x <= pointRefOfCommonXmin.x)
        {
          trefMiddle = trefOfCommonXmin;
          pointRefOnMiddle = pointRefOfCommonXmin;
        }
        else if (pointOthOnMiddle.x >= pointRefOfCommonXmax.x)
        {
          trefMiddle = trefOfCommonXmax;
          pointRefOnMiddle = pointRefOfCommonXmax;
        }
        else
        {
          auto func = [&referenceCurve, xo = pointOthOnMiddle.x](const double t) -> double
          {
            return referenceCurve.getPoint(t).x - xo;
          };
          const double trefMinCurrent = std::min(trefOfCommonXmin, trefOfCommonXmax);
          const double trefMaxCurrent = std::max(trefOfCommonXmin, trefOfCommonXmax);
          trefMiddle = math::findUniqueFunctionRoot(trefMinCurrent, trefMaxCurrent, func);
          pointRefOnMiddle = referenceCurve.getPoint(trefMiddle);
        }
        // ��������� ��������
        if (point::isSame(pointOthOnMiddle, pointRefOnMiddle))
        {
          if (theLargestCurve == 0)
          {
            return IntersecctionSolutionType{ 0.5 * (pointOthOnMiddle + pointRefOnMiddle), trefMiddle, tothMiddle};
          }
          else
          {
            return IntersecctionSolutionType{ 0.5 * (pointOthOnMiddle + pointRefOnMiddle), tothMiddle, trefMiddle };
          }
        }
        // �������� �� �������� ��������, ����� ������� ���� �����, ���� ������
        fleft = pointOthOfCommonXmin.y - pointRefOfCommonXmin.y;
        fright = pointOthOfCommonXmax.y - pointRefOfCommonXmax.y;
        const double fmiddle = pointOthOnMiddle.y - pointRefOnMiddle.y;
        if (((fleft > 0.0) and (fmiddle < 0.0)) or ((fleft < 0.0) and (fmiddle > 0.0)))
        {
          // �������� ���������� ������ ������, ��� ��� ����������� ��������� �����
          pointOthOfCommonXmax = pointOthOnMiddle;
          tothOfCommonXmax     = tothMiddle;
          pointRefOfCommonXmax = pointRefOnMiddle;
          trefOfCommonXmax     = trefMiddle;
        }
        else if (((fright > 0.0) and (fmiddle < 0.0)) or ((fright < 0.0) and (fmiddle > 0.0)))
        {
          // �������� ���������� ����� ������, ��� ��� ����������� ��������� ������
          pointOthOfCommonXmin = pointOthOnMiddle;
          tothOfCommonXmin     = tothMiddle;
          pointRefOfCommonXmin = pointRefOnMiddle;
          trefOfCommonXmin     = trefMiddle;
        }
        else
        {
          throw std::logic_error("impossible situation");
        }
      }
    }
  }
  else
  {

  }
#endif
  return std::nullopt;
}

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
  const auto Pref = referenceCurve.getPoint(trefmin);
  const auto Qref = referenceCurve.getPoint(trefmax);
  StatOfCurvePiece scpRef{ trefmin, Pref, trefmax, Qref };

  const auto Poth = otherCurve.getPoint(trefmin);
  const auto Qoth = otherCurve.getPoint(trefmax);
  StatOfCurvePiece scpOth{ tothmin, Poth, tothmax, Qoth };

  if (scpOth.pointOfxmax().x <= scpRef.pointOfxmin().x)
  {
    if (point::isSame(scpOth.pointOfxmax(), scpRef.pointOfxmin()))
    {
      return IntersecctionSolutionType{ 0.5 * (scpOth.pointOfxmax() + scpRef.pointOfxmin()), scpRef.tOfxmin(), scpOth.tOfxmax() };
    }
    else
    {
      return std::nullopt;
    }
  }
  else if (scpOth.pointOfxmin().x >= scpRef.pointOfxmax().x)
  {
    if (point::isSame(scpOth.pointOfxmin(), scpRef.pointOfxmax()))
    {
      return IntersecctionSolutionType{ 0.5 * (scpOth.pointOfxmin() + scpRef.pointOfxmax()), scpRef.tOfxmax(), scpOth.tOfxmin() };
    }
    else
    {
      return std::nullopt;
    }
  }
  else
  {
    // MIN X
    const DataGetterOfX referenceGetter{ trefmin, trefmax, referenceCurve };
    const DataGetterOfX otherGetter{ tothmin, tothmax, otherCurve };
    point pointRefOfCommonXmin;
    auto [commonXmin, trefOfCommonXmin, pointRefOfCommonXmin, tothOfCommonXmin, pointOthOfCommonXmin] =
      CommonRangeHelper::ofLowest(referenceGetter, otherGetter);

    double commonXmin = 0.0;
    double trefOfCommonXmin = 0.0;
    point pointRefOfCommonXmin;
    double tothOfCommonXmin = 0.0;
    point pointOthOfCommonXmin;
    if (scpRef.pointOfxmin().x > scpOth.pointOfxmin().x)
    {
      commonXmin = scpRef.pointOfxmin().x;
      trefOfCommonXmin = scpRef.tOfxmin();
      auto func = [&otherCurve, xo = commonXmin](const double t) -> double
      {
        return otherCurve.getPoint(t).x - xo;
      };
      tothOfCommonXmin = math::findUniqueFunctionRoot(tothmin, tothmax, func);
      pointOthOfCommonXmin = otherCurve.getPoint(tothOfCommonXmin);
    }
    else if (scpRef.pointOfxmin().x < scpOth.pointOfxmin().x)
    {
      commonXmin = scpOth.pointOfxmin().x;
      tothOfCommonXmin = scpOth.tOfxmin();
      auto func = [&referenceCurve, xo = commonXmin](const double t) -> double
      {
        return referenceCurve.getPoint(t).x - xo;
      };
      trefOfCommonXmin = math::findUniqueFunctionRoot(trefmin, trefmax, func);
      pointRefOfCommonXmin = referenceCurve.getPoint(trefOfCommonXmin);
    }
    else
    {
      commonXmin = scpRef.pointOfxmin().x;
      trefOfCommonXmin = scpRef.tOfxmin();
      pointRefOfCommonXmin = scpRef.pointOfxmin();
      tothOfCommonXmin = scpOth.tOfxmin();
      pointOthOfCommonXmin = scpOth.pointOfxmin();
    }
    // MAX X
    double commonXmax = 0.0;
    double trefOfCommonXmax = 0.0;
    point pointRefOfCommonXmax;
    double tothOfCommonXmax = 0.0;
    point pointOthOfCommonXmax;
    if (scpRef.pointOfxmax().x < scpOth.pointOfxmax().x)
    {
      commonXmax = scpRef.pointOfxmax().x;
      trefOfCommonXmax = scpRef.tOfxmax();
      auto func = [&otherCurve, xo = commonXmax](const double t) -> double
      {
        return otherCurve.getPoint(t).x - xo;
      };
      tothOfCommonXmax = math::findUniqueFunctionRoot(tothmin, tothmax, func);
      pointOthOfCommonXmax = otherCurve.getPoint(tothOfCommonXmax);
    }
    else if (scpRef.pointOfxmax().x > scpOth.pointOfxmax().x)
    {
      commonXmax = scpOth.pointOfxmax().x;
      tothOfCommonXmax = scpOth.tOfxmax();
      auto func = [&referenceCurve, xo = commonXmax](const double t) -> double
      {
        return referenceCurve.getPoint(t).x - xo;
      };
      trefOfCommonXmax = math::findUniqueFunctionRoot(trefmin, trefmax, func);
      pointRefOfCommonXmax = referenceCurve.getPoint(trefOfCommonXmax);
    }
    else
    {
      commonXmax = scpRef.pointOfxmax().x;
      trefOfCommonXmax = scpRef.tOfxmax();
      pointRefOfCommonXmax = scpRef.pointOfxmax();
      tothOfCommonXmax = scpOth.tOfxmax();
      pointOthOfCommonXmax = scpOth.pointOfxmax();
    }
    // ��������� ����� ��� �� X
    //  (*) ��������� ����� �����        
    if (point::isSame(pointOthOfCommonXmin, pointRefOfCommonXmin))
    {
      return IntersecctionSolutionType{ 0.5 * (pointOthOfCommonXmin + pointRefOfCommonXmin), trefOfCommonXmin, tothOfCommonXmin };
    }
    //  (*) ��������� ������ �����
    if (point::isSame(pointOthOfCommonXmax, pointRefOfCommonXmax))
    {
      return IntersecctionSolutionType{ 0.5 * (pointOthOfCommonXmax + pointRefOfCommonXmax), trefOfCommonXmax, tothOfCommonXmax };
    }
    // ����� �� ���� ������ ��������� ���� ��� ���� ������
    double fleft = pointOthOfCommonXmin.y - pointRefOfCommonXmin.y;
    double fright = pointOthOfCommonXmax.y - pointRefOfCommonXmax.y;
    const bool othIsAboveRef = (fleft > 0.0) and (fright > 0.0);
    const bool refIsUnderRef = (fleft < 0.0) and (fright < 0.0);
    if (othIsAboveRef or refIsUnderRef)
    {
      return std::nullopt;
    }
    while (true)
    {
      double tothMiddle = 0.5 * (tothOfCommonXmin + tothOfCommonXmax);
      point pointOthOnMiddle = otherCurve.getPoint(tothMiddle);

      double trefMiddle = 0.0;
      point pointRefOnMiddle;
      if (pointOthOnMiddle.x <= pointRefOfCommonXmin.x)
      {
        trefMiddle = trefOfCommonXmin;
        pointRefOnMiddle = pointRefOfCommonXmin;
      }
      else if (pointOthOnMiddle.x >= pointRefOfCommonXmax.x)
      {
        trefMiddle = trefOfCommonXmax;
        pointRefOnMiddle = pointRefOfCommonXmax;
      }
      else
      {
        auto func = [&referenceCurve, xo = pointOthOnMiddle.x](const double t) -> double
        {
          return referenceCurve.getPoint(t).x - xo;
        };
        const double trefMinCurrent = std::min(trefOfCommonXmin, trefOfCommonXmax);
        const double trefMaxCurrent = std::max(trefOfCommonXmin, trefOfCommonXmax);
        trefMiddle = math::findUniqueFunctionRoot(trefMinCurrent, trefMaxCurrent, func);
        pointRefOnMiddle = referenceCurve.getPoint(trefMiddle);
      }
      // ��������� ��������
      if (point::isSame(pointOthOnMiddle, pointRefOnMiddle))
      {
        return IntersecctionSolutionType{ 0.5 * (pointOthOnMiddle + pointRefOnMiddle), trefMiddle, tothMiddle };
      }
      // �������� �� �������� ��������, ����� ������� ���� �����, ���� ������
      fleft = pointOthOfCommonXmin.y - pointRefOfCommonXmin.y;
      fright = pointOthOfCommonXmax.y - pointRefOfCommonXmax.y;
      const double fmiddle = pointOthOnMiddle.y - pointRefOnMiddle.y;
      if (((fleft > 0.0) and (fmiddle < 0.0)) or ((fleft < 0.0) and (fmiddle > 0.0)))
      {
        // �������� ���������� ������ ������, ��� ��� ����������� ��������� �����
        pointOthOfCommonXmax = pointOthOnMiddle;
        tothOfCommonXmax = tothMiddle;
        pointRefOfCommonXmax = pointRefOnMiddle;
        trefOfCommonXmax = trefMiddle;
      }
      else if (((fright > 0.0) and (fmiddle < 0.0)) or ((fright < 0.0) and (fmiddle > 0.0)))
      {
        // �������� ���������� ����� ������, ��� ��� ����������� ��������� ������
        pointOthOfCommonXmin = pointOthOnMiddle;
        tothOfCommonXmin = tothMiddle;
        pointRefOfCommonXmin = pointRefOnMiddle;
        trefOfCommonXmin = trefMiddle;
      }
      else
      {
        throw std::logic_error("impossible situation in functin findUniqueIntersectionRefAlongX");
      }
    }
  }
  return std::nullopt;
}

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

}

std::tuple<double, double, double>
  geom2d::curveIntersector::findLargestCommonXmin
  (
    const double tmin1,
    const double tmax1,
    const baseCurve& curve1,
    const double tmin2,
    const double tmax2,
    const baseCurve& curve2
  )
{
  const point p1 = curve1.getPoint(tmin1);
  const point q1 = curve1.getPoint(tmax1);
  StatOfCurvePiece scp1{ tmin1, p1, tmax1, q1 };

  const point p2 = curve2.getPoint(tmin1);
  const point q2 = curve2.getPoint(tmax1);
  StatOfCurvePiece scp2{tmin2, p2, tmax2, q2};

  if (scp1.pointOfxmin().x > scp2.pointOfxmin().x)
  {
    const double commonXmin = scp1.pointOfxmin().x;
    const double tleft1 = scp1.tOfxmin();
    auto func = [&curve2, xo = commonXmin](const double t) -> double
    {
      return curve2.getPoint(t).x - xo;
    };
    const double tleft2 = math::findUniqueFunctionRoot(tmin2, tmax2, func);
    return std::tuple{ commonXmin, tleft1, tleft2 };
  }
  else if (scp1.pointOfxmin().x < scp2.pointOfxmin().x)
  {
    const double commonXmin = scp2.pointOfxmin().x;
    const double tleft2 = scp2.tOfxmin();
    auto func = [&curve1, xo = commonXmin](const double t) -> double
    {
      return curve1.getPoint(t).x - xo;
    };
    const double tleft1 = math::findUniqueFunctionRoot(tmin1, tmax1, func);
    return std::tuple{ commonXmin, tleft1, tleft2 };
  }
  else
  {
    const double commonXmin = scp1.pointOfxmin().x;
    const double tleft1 = scp1.tOfxmin();
    const double tleft2 = scp2.tOfxmin();
    return std::tuple{ commonXmin, tleft1, tleft2 };
  }
}

//-----------------------------------------------------------------------------
