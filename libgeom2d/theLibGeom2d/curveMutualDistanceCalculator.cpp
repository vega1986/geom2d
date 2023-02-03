#include "curveMutualDistanceCalculator.h"
#include "StatOfCurvePiece.h"
#include "distanceCalculator_point_and_plato.h"
#include "distanceCalculator_two_parallel_plato.h"
#include "point_and_curve_extrema.h"
#include "segmentCurve.h"
#include "curveIntersector.h"

//*********************************************************************************************************************

void geom2d::curveMutualDistanceCalculator::fulfill()
{
  curveAnalizerBase::perform();
}

//*********************************************************************************************************************

const std::tuple<double, double, double> geom2d::curveMutualDistanceCalculator::getExtrema() const
{
  return std::tuple{m_dist, m_t1, m_t2};
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPlatoY1PlatoX2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT2, theT1] = execPlatoX_and_PlatoY(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);

  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPlatoY1PlatoY2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // points of curve 1
  for (const auto t1 : std::initializer_list{ tmin1, tmax1 })
  {
    const auto pnt1 = m_curve1.getPoint(t1);

    using namespace geom2d::extrema::point_and_plato;
    distanceCalculator<DataGetterOfX> calculator{ pnt1, tmin2, tmax2, m_curve2 };

    const auto [dist, ton2] = calculator.execute();
    if (dist < m_dist)
    {
      m_dist = dist;
      m_t1 = t1;
      m_t2 = ton2;
    }
  }
  // points of curve 2
  for (const auto t2 : std::initializer_list{ tmin2, tmax2 })
  {
    const auto pnt2 = m_curve2.getPoint(t2);

    using namespace geom2d::extrema::point_and_plato;
    distanceCalculator<DataGetterOfX> calculator{ pnt2, tmin1, tmax1, m_curve1 };

    const auto [dist, ton1] = calculator.execute();
    if (dist < m_dist)
    {
      m_dist = dist;
      m_t1 = ton1;
      m_t2 = t2;
    }
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPlatoY1Point2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT2, theT1] = execPoint_and_PlatoY(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPoint1Any2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT1, theT2] = execPoint_and_Any(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPoint1PlatoX2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT1, theT2] = execPoint_and_PlatoX(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPoint1PlatoY2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT1, theT2] = execPoint_and_PlatoY(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPoint1Point2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const std::initializer_list tsof1{ tmin1, tmax1 };
  const std::initializer_list tsof2{ tmin2, tmax2 };

  for (const auto t1 : tsof1)
  {
    const auto P1 = m_curve1.getPoint(t1);
    for (const auto t2 : tsof2)
    {
      const auto P2 = m_curve2.getPoint(t2);
      const auto curDist = point::distance(P1, P2);
      if (curDist < m_dist)
      {
        m_dist = curDist;
        m_t1 = t1;
        m_t2 = t2;
      }
    }
  }
}

//*********************************************************************************************************************

const std::tuple<double, double, double>
  geom2d::curveMutualDistanceCalculator::execPoint_and_PlatoY
  (
    const double tMinOfPoint,
    const double tMaxOfPoint,
    const baseCurve& curvePoint,
    const double tMinOfPlatoY,
    const double tMaxOfPlatoY,
    const baseCurve& curvePlatoY
  )
{
  using namespace geom2d::extrema::point_and_plato;

  double theDist = math::infinite::distance;
  double theTpoint = 0.0;
  double theTplato = 0.0;

  for (const auto tpoint : std::initializer_list{ tMinOfPoint , tMaxOfPoint })
  {
    const auto pnt = curvePoint.getPoint(tpoint);
    distanceCalculator<DataGetterOfX> calculator{ pnt , tMinOfPlatoY , tMaxOfPlatoY , curvePlatoY };
    const auto [dist, tofPlato] = calculator.execute();

    if (dist < theDist)
    {
      theDist = dist;
      theTpoint = tpoint;
      theTplato = tofPlato;
    }
  }
  return std::tuple{theDist, theTpoint, theTplato};
}

//*********************************************************************************************************************

const std::tuple<double, double, double>
  geom2d::curveMutualDistanceCalculator::execPoint_and_PlatoX
  (
    const double tMinOfPoint,
    const double tMaxOfPoint,
    const baseCurve& curvePoint,
    const double tMinOfPlatoX,
    const double tMaxOfPlatoX,
    const baseCurve& curvePlatoX
  )
{
  using namespace geom2d::extrema::point_and_plato;

  double theDist = math::infinite::distance;
  double theTpoint = 0.0;
  double theTplato = 0.0;

  for (const auto tpoint : std::initializer_list{ tMinOfPoint , tMaxOfPoint })
  {
    const auto pnt = curvePoint.getPoint(tpoint);
    distanceCalculator<DataGetterOfY> calculator{ pnt , tMinOfPlatoX , tMaxOfPlatoX , curvePlatoX };
    const auto [dist, tofPlato] = calculator.execute();

    if (dist < theDist)
    {
      theDist = dist;
      theTpoint = tpoint;
      theTplato = tofPlato;
    }
  }
  return std::tuple{ theDist, theTpoint, theTplato };
}

//*********************************************************************************************************************

const std::tuple<double, double, double>
  geom2d::curveMutualDistanceCalculator::execPoint_and_Any
  (
    const double tMinOfPoint,
    const double tMaxOfPoint,
    const baseCurve& curvePoint,
    const double tMinOfAny,
    const double tMaxOfAny,
    const baseCurve& curveAny
  )
{
  double theDist = math::infinite::distance;
  double theTofPoint = 0.0;
  double theTofAny = 0.0;

  for (const auto tp : std::initializer_list{ tMinOfPoint , tMaxOfPoint })
  {
    const auto pnt = curvePoint.getPoint(tp);

    using namespace geom2d::extrema::point_and_curve;
    nearestPoint np{ pnt, tMinOfAny , tMaxOfAny , curveAny };
    const auto [dist, tofAny] = np.execute();

    if (dist < theDist)
    {
      theDist = dist;
      theTofPoint = tp;
      theTofAny = tofAny;
    }
  }
  return std::tuple{ theDist, theTofPoint , theTofAny };
}

//*********************************************************************************************************************

const std::tuple<double, double, double>
  geom2d::curveMutualDistanceCalculator::execPlatoX_and_PlatoY
  (
    const double tminPlatoX,
    const double tmaxPlatoX,
    const baseCurve& curvePlatoX,
    const double tminPlatoY,
    const double tmaxPlatoY,
    const baseCurve& curvePlatoY
  )
{
  double localDist = math::infinite::distance;
  double localTofPlatoX = 0.0;
  double localTofPlatoY = 0.0;

  for (const auto tofPlatoX : std::initializer_list{ tminPlatoX, tmaxPlatoX })
  {
    const auto pnt = curvePlatoX.getPoint(tofPlatoX);

    using namespace geom2d::extrema::point_and_plato;
    distanceCalculator<DataGetterOfX> calculator{ pnt , tminPlatoY , tmaxPlatoY, curvePlatoY };
    const auto [dist, tofPlatoY] = calculator.execute();

    if (dist < localDist)
    {
      localDist = dist;
      localTofPlatoX = tofPlatoX;
      localTofPlatoY = tofPlatoY;
    }
  }
  for (const auto tofPlatoY : std::initializer_list{ tminPlatoY, tmaxPlatoY })
  {
    const auto pnt = curvePlatoY.getPoint(tofPlatoY);

    using namespace geom2d::extrema::point_and_plato;
    distanceCalculator<DataGetterOfY> calculator{ pnt , tminPlatoX , tmaxPlatoX, curvePlatoX };
    const auto [dist, tofPlatoX] = calculator.execute();

    if (dist < localDist)
    {
      localDist = dist;
      localTofPlatoX = tofPlatoX;
      localTofPlatoY = tofPlatoY;
    }
  }
  return std::tuple{ localDist , localTofPlatoX , localTofPlatoY };
}

//*********************************************************************************************************************

const std::tuple<double, double, double>
  geom2d::curveMutualDistanceCalculator::execPlatoX_and_Any
  (
    const double tminPlatoX,
    const double tmaxPlatoX,
    const baseCurve& curvePlatoX,
    const double tminAny,
    const double tmaxAny,
    const baseCurve& curveAny
  )
{
  double theDist = math::infinite::distance;
  double tExtremaOfPlato = 0.0;
  double tExtremaOfAny = 0.0;

  for (const auto tofplato : std::initializer_list{ tminPlatoX , tmaxPlatoX })
  {
    const auto pnt = curvePlatoX.getPoint(tofplato);

    using namespace geom2d::extrema::point_and_curve;
    nearestPoint np{ pnt, tminAny, tmaxAny, curveAny };

    const auto [dist, tonAny] = np.execute();
    if (dist < theDist)
    {
      theDist = dist;
      tExtremaOfPlato = tofplato;
      tExtremaOfAny = tonAny;
    }
  }

  for (const auto tofany : std::initializer_list{ tminAny , tmaxAny })
  {
    const auto pnt = curveAny.getPoint(tofany);

    using namespace geom2d::extrema::point_and_plato;
    distanceCalculator<DataGetterOfY> calculator{ pnt, tminPlatoX, tmaxPlatoX, curvePlatoX };

    const auto [dist, tonPlatoX] = calculator.execute();
    if (dist < theDist)
    {
      theDist = dist;
      tExtremaOfPlato = tonPlatoX;
      tExtremaOfAny = tofany;
    }
  }
  return std::tuple{ theDist , tExtremaOfPlato , tExtremaOfAny };
}

//*********************************************************************************************************************

const std::tuple<double, double, double>
  geom2d::curveMutualDistanceCalculator::execPlatoY_and_Any
  (
    const double tminPlatoY,
    const double tmaxPlatoY,
    const baseCurve& curvePlatoY,
    const double tminAny,
    const double tmaxAny,
    const baseCurve& curveAny
  )
{
  double theDist = math::infinite::distance;
  double tExtremaOfPlato = 0.0;
  double tExtremaOfAny = 0.0;

  for (const auto tofplato : std::initializer_list{ tminPlatoY , tmaxPlatoY })
  {
    const auto pnt = curvePlatoY.getPoint(tofplato);

    using namespace geom2d::extrema::point_and_curve;
    nearestPoint np{ pnt, tminAny, tmaxAny, curveAny };

    const auto [dist, tonAny] = np.execute();
    if (dist < theDist)
    {
      theDist = dist;
      tExtremaOfPlato = tofplato;
      tExtremaOfAny = tonAny;
    }
  }

  for (const auto tofany : std::initializer_list{ tminAny , tmaxAny })
  {
    const auto pnt = curveAny.getPoint(tofany);

    using namespace geom2d::extrema::point_and_plato;
    distanceCalculator<DataGetterOfX> calculator{ pnt, tminPlatoY, tmaxPlatoY, curvePlatoY };

    const auto [dist, tonPlatoX] = calculator.execute();
    if (dist < theDist)
    {
      theDist = dist;
      tExtremaOfPlato = tonPlatoX;
      tExtremaOfAny = tofany;
    }
  }
  return std::tuple{ theDist , tExtremaOfPlato , tExtremaOfAny };
}

//*********************************************************************************************************************

const std::tuple<double, double, double>
  geom2d::curveMutualDistanceCalculator::execAnyOfDifferentOrientation
  (
    const double tmin1,
    const double tmax1,
    const baseCurve& curve1,
    const double tmin2,
    const double tmax2,
    const baseCurve& curve2
  )
{
  double theDist = math::infinite::distance;
  double theT1 = 0.0;
  double theT2 = 0.0;
  
  for (const auto t1 : std::initializer_list{ tmin1, tmax1 })
  {
    const point P1 = curve1.getPoint(t1);

    using namespace geom2d::extrema::point_and_curve;
    nearestPoint np{ P1, tmin2, tmax2, curve2 };

    const auto [dist, ton2] = np.execute();
    if (dist < theDist)
    {
      theDist = dist;
      theT1 = t1;
      theT2 = ton2;
    }
  }

  for (const auto t2 : std::initializer_list{ tmin2, tmax2 })
  {
    const point P2 = curve2.getPoint(t2);

    using namespace geom2d::extrema::point_and_curve;
    nearestPoint np{ P2, tmin1, tmax1, curve1 };

    const auto [dist, ton1] = np.execute();
    if (dist < theDist)
    {
      theDist = dist;
      theT1 = ton1;
      theT2 = t2;
    }
  }
  return std::tuple{ theDist , theT1 , theT2 };
}

//*********************************************************************************************************************

const std::tuple<double, double, double>
  geom2d::curveMutualDistanceCalculator::execAnyOfSameOrientation
  (
    const double tmin1,
    const double tmax1,
    const baseCurve& curve1,
    const double tmin2,
    const double tmax2,
    const baseCurve& curve2
  )
{
  double theDist = math::infinite::distance;
  double theT1 = 0.0;
  double theT2 = 0.0;
  
  {
    using namespace geom2d::extrema::point_and_curve;
    
    for (const auto ton1 : std::initializer_list{tmin1, tmax1})
    {
      const auto pon1 = curve1.getPoint(ton1);
      nearestPoint np{ pon1, tmin2, tmax2, curve2 };
      const auto [dist, ton2] = np.execute();

      if (dist < theDist)
      {
        theDist = dist;
        theT1 = ton1;
        theT2 = ton2;
      }
    }

    for (const auto ton2 : std::initializer_list{ tmin2, tmax2 })
    {
      const auto pon2 = curve2.getPoint(ton2);
      nearestPoint np{ pon2, tmin1, tmax1, curve1 };
      const auto [dist, ton1] = np.execute();

      if (dist < theDist)
      {
        theDist = dist;
        theT1 = ton1;
        theT2 = ton2;
      }
    }
  }

  constexpr size_t n = 1000; // number of split of curves

  aabb aabbOfCurve2{ std::initializer_list{curve2.getPoint(tmin2), curve2.getPoint(tmax2)} };

  const auto xmin2 = aabbOfCurve2.xmin();
  const auto xmax2 = aabbOfCurve2.xmax();
  const auto ymin2 = aabbOfCurve2.ymin();
  const auto ymax2 = aabbOfCurve2.ymax();

  auto stopCriteria = [&curve1](const double ta, const double tb) -> bool
  {
    const auto pa = curve1.getPoint(ta);
    const auto pb = curve1.getPoint(ta);
    return point::isSame(pa, pb);
  };

  // функтор для расчёта скалярного произведения касательного (единичного) вектора нормали к кривой 2
  // и вектора, соединяющего точку на кривой 2 с точкой на кривой 1
  // по сути получаем длину проекции точки на кривой 1 на касательную к кривой 2
  auto calcScalarMult = [&curve1, &curve2](const double t1, const double t2) -> double
  {
    const auto p1 = curve1.getPoint(t1);

    const auto p2 = curve2.getPoint(t2);
    const auto velocity = curve2.getVelocity(t2);
    const auto tau2{ vector{velocity}.normalized() };

    const vector B{ p2, p1 };
    return (tau2, B.normalized());
  };

  auto
    calcPartnerT =
      [&curve1, &curve2, tmin2, tmax2, xmin2, xmax2, ymin2, ymax2]
        (const double t1) -> std::optional<double>
  {
    // строим отрезок нормали к кривой 1 в точке P1
    const auto segment = buildSegmentCurveAsNormalToCurve(t1, curve1, xmin2, xmax2, ymin2, ymax2);
    // получаем гипотетическую точку пересечения этого отрезка с кривой 2
    return getTofIntersectionWithSegment(segment, tmin2, tmax2, curve2);
  };

  auto
    calcScalarMultFunc =
      [&curve1, &curve2, &calcScalarMult, &calcPartnerT/*, tmin2, tmax2, xmin2, xmax2, ymin2, ymax2*/]
        (const double t1) -> double
  {
    const auto solution = calcPartnerT(t1);
    if (not solution)
    {
      throw std::logic_error("calcScalarMultFunc: solution of intersection segment with curve2 does not exists");
    }

    const auto t2 = solution.value();
    return calcScalarMult(t1, t2);
  };
  
  for (size_t j = 0; j < n; ++j)
  {
    const double tP1 = tmin1 + (tmax1 - tmin1) * static_cast<double>(j    ) / static_cast<double>(n);
    const double tQ1 = tmin1 + (tmax1 - tmin1) * static_cast<double>(j + 1) / static_cast<double>(n);

    // строим отрезок нормали к кривой 1 в точке P1
    //const auto segNormForP = buildSegmentCurveAsNormalToCurve(tP1, curve1, xmin2, xmax2, ymin2, ymax2);
    // получаем гипотетическую точку пересечения этого отрезка с кривой 2
    const auto solP1 = calcPartnerT(tP1);
    //const auto solP1 = getTofIntersectionWithSegment(segNormForP, tmin2, tmax2, curve2);

    // строим отрезок нормали к кривой 1 в точке Q1
    //const auto segNormForQ = buildSegmentCurveAsNormalToCurve(tQ1, curve1, xmin2, xmax2, ymin2, ymax2);
    // получаем гипотетическую точку пересечения этого отрезка с кривой 2
    const auto solQ1 = calcPartnerT(tQ1);
    //const auto solQ1 = getTofIntersectionWithSegment(segNormForQ, tmin2, tmax2, curve2);

    if ((not solP1) and (not solQ1))
    {
      continue;
    }

    if (solP1 and solQ1)
    {
      const auto tP2 = solP1.value();
      const auto tQ2 = solQ1.value();

      const double scalarMultP = calcScalarMult(tP1, tP2);
      const double scalarMultQ = calcScalarMult(tQ1, tQ2);

      // проверяем точки P на предмет достижения в них экстремума
      if (std::abs(scalarMultP) <= math::tolerance::tolScalarMult)
      {
        const auto pnt1 = curve1.getPoint(tP1);
        const auto pnt2 = curve2.getPoint(tP2);
        const auto dist = point::distance(pnt1, pnt2);

        if (dist < theDist)
        {
          theDist = dist;
          theT1 = tP1;
          theT2 = tP2;
        }
        // экстрема найдена - заканчиваем шаг цикла
        continue;
      }

      // проверяем точки Q на предмет достижения в них экстремума
      if (std::abs(scalarMultQ) <= math::tolerance::tolScalarMult)
      {
        const auto pnt1 = curve1.getPoint(tQ1);
        const auto pnt2 = curve2.getPoint(tQ2);
        const auto dist = point::distance(pnt1, pnt2);

        if (dist < theDist)
        {
          theDist = dist;
          theT1 = tQ1;
          theT2 = tQ2;
        }
        // экстрема найдена - заканчиваем шаг цикла
        continue;
      }

      if (((scalarMultP > 0.0) and (scalarMultQ > 0.0)) or ((scalarMultP < 0.0) and (scalarMultQ < 0.0)))
      {
        // экстрем нет на отрезке [tP1; tQ1] x [tP2; tQ2]
        continue;
      }
      
      // экстрема где-то между этими точками - решаем нелинейное уравнение относительно парамета t на кривой 1
      // используем для этого функтор одной переменной calcScalarMultFunc
      // вопрос - какую точность здесь мы хотим ?
      // точку пересечения мы определяем с точностью 0.1 от стандартного толеранса точки - это стандартное значение.
      // возможно здесь следует потребовать меньшую точность - 1.0 от стандартного толеранса точки.
      const double tofExtrema1 = math::findUniqueFunctionRoot(tP1, tQ1, calcScalarMultFunc, math::tolerance::tolScalarMult);
        //math::findUniqueFunctionRoot(tP1, tQ1, calcScalarMultFunc);
      
      // строим отрезок нормали к кривой 1
      //const auto segExtrema = buildSegmentCurveAsNormalToCurve(tofExtrema1, curve1, xmin2, xmax2, ymin2, ymax2);
      // получаем гипотетическую точку пересечения этого отрезка с кривой 2
      const auto solExtrema = calcPartnerT(tofExtrema1);
        //getTofIntersectionWithSegment(segExtrema, tmin2, tmax2, curve2);

      if (not solExtrema)
      {
        throw std::logic_error("geom2d::curveMutualDistanceCalculator::execAnyOfSameOrientation: impossible (solExtrema does not exist)");
      }

      const auto tofExtrema2 = solExtrema.value();

      const auto pofExtrema1 = curve1.getPoint(tofExtrema1);
      const auto pofExtrema2 = curve2.getPoint(tofExtrema2);

      const auto dist = point::distance(pofExtrema1, pofExtrema2);

      if (dist < theDist)
      {
        theDist = dist;
        theT1 = tofExtrema1;
        theT2 = tofExtrema2;
      }
      continue;
    }

    if (solP1)
    {
      const auto tP2 = solP1.value();
      const double scalarMultP = calcScalarMult(tP1, tP2);

      double ta = tP1;
      double tb = tQ1;

      auto sola = solP1; // +
      auto solb = solQ1; // -

      auto pa = curve1.getPoint(ta);
      auto pb = curve1.getPoint(tb);

      while (point::distance(pa, pb) <= math::tolerance::tolPoint)
      {
        const double t = 0.5 * (ta + tb);
        const auto sol = calcPartnerT(t);
        if (sol)
        {
          ta = t;
          pa = curve1.getPoint(ta);
          sola = sol;
        }
        else
        {
          tb = t;
          pb = curve1.getPoint(tb);
          solb = sol;
        }
      }
      
      const double tPnew1 = tP1;
      const double tQnew1 = ta;
      const auto pQnew1 = pa;

      const double tPnew2 = tP2;
      const double tQnew2 = sola.value();
      const auto pQnew2 = curve2.getPoint(tQnew2);

      const double scalarMultQnew = calcScalarMult(tQnew1, tQnew2);
      if (std::abs(scalarMultQnew) <= math::tolerance::tolScalarMult)
      {
        const auto dist = point::distance(pQnew1, pQnew2);
        if (dist < theDist)
        {
          theDist = dist;
          theT1 = tQnew1;
          theT2 = tQnew2;
        }

        continue;
      }

      if (((scalarMultP > 0.0) and (scalarMultQnew > 0.0)) or ((scalarMultP < 0.0) and (scalarMultQnew < 0.0)))
      {
        // экстрем на отрезке [tP1; tQnew1] x [tP2; tQnew2]
        continue;
      }

      const auto tofExtremaOn1 = math::findUniqueFunctionRoot(tPnew1, tQnew1, calcScalarMultFunc, math::tolerance::tolScalarMult);
      const auto solExtrema = calcPartnerT(tofExtremaOn1);

      if (not solExtrema)
      {
        throw std::logic_error("geom2d::curveMutualDistanceCalculator::execAnyOfSameOrientation: impossible (solExtrema does not exist)");
      }

      const auto tofExtremaOn2 = solExtrema.value();
      const point pofExtremaOn1 = curve1.getPoint(tofExtremaOn1);
      const point pofExtremaOn2 = curve2.getPoint(tofExtremaOn2);

      const auto dist = point::distance(pofExtremaOn1, pofExtremaOn2);
      if (dist < theDist)
      {
        theDist = dist;
        theT1 = tofExtremaOn1;
        theT2 = tofExtremaOn2;
      }

      continue;
    }

    if (solQ1)
    {
      const auto tQ2 = solQ1.value();
      const double scalarMultQ = calcScalarMult(tQ1, tQ2);

      double ta = tP1;
      double tb = tQ1;

      auto sola = solP1; // -
      auto solb = solQ1; // +

      auto pa = curve1.getPoint(ta);
      auto pb = curve1.getPoint(tb);

      while (point::distance(pa, pb) <= math::tolerance::tolPoint)
      {
        const double t = 0.5 * (ta + tb);
        const auto sol = calcPartnerT(t);
        if (sol)
        {
          tb = t;
          pb = curve1.getPoint(tb);
          solb = sol;
        }
        else
        {
          ta = t;
          pa = curve1.getPoint(ta);
          sola = sol;
        }
      }

      const double tPnew1 = tb;
      const double tQnew1 = tQ1;
      const auto pPnew1 = pb;

      const double tPnew2 = solb.value();
      const double tQnew2 = tQ2;
      const auto pPnew2 = curve2.getPoint(tPnew2);

      const double scalarMultPnew = calcScalarMult(tPnew1, tPnew2);
      if (std::abs(scalarMultPnew) <= math::tolerance::tolScalarMult)
      {
        const auto dist = point::distance(pPnew1, pPnew2);
        if (dist < theDist)
        {
          theDist = dist;
          theT1 = tPnew1;
          theT2 = tPnew2;
        }

        continue;
      }

      if (((scalarMultQ > 0.0) and (scalarMultPnew > 0.0)) or ((scalarMultQ < 0.0) and (scalarMultPnew < 0.0)))
      {
        // экстрем на отрезке [tPnew1; tQ1] x [tPnew2; tQ2]
        continue;
      }

      const auto tofExtremaOn1 = math::findUniqueFunctionRoot(tPnew1, tQnew1, calcScalarMultFunc, math::tolerance::tolScalarMult);
      const auto solExtrema = calcPartnerT(tofExtremaOn1);

      if (not solExtrema)
      {
        throw std::logic_error("geom2d::curveMutualDistanceCalculator::execAnyOfSameOrientation: impossible (solExtrema does not exist)");
      }

      const auto tofExtremaOn2 = solExtrema.value();
      const point pofExtremaOn1 = curve1.getPoint(tofExtremaOn1);
      const point pofExtremaOn2 = curve2.getPoint(tofExtremaOn2);

      const auto dist = point::distance(pofExtremaOn1, pofExtremaOn2);
      if (dist < theDist)
      {
        theDist = dist;
        theT1 = tofExtremaOn1;
        theT2 = tofExtremaOn2;
      }

      continue;
    }
  }
  return std::tuple{ theDist , theT1 , theT2 };
}

//*********************************************************************************************************************

std::optional<double>
  geom2d::curveMutualDistanceCalculator::getTofIntersectionWithSegment
  (
    const segmentCurve& segment,
    const double tmin,
    const double tmax,
    const baseCurve& curve
  )
{
  const double tminSeg = segment.parameterMin();
  const double tmaxSeg = segment.parameterMax();
  const auto theIntersection = curveIntersector::exec_Any_and_Any_Unique_Intersection(tminSeg, tmaxSeg, segment, tmin, tmax, curve);
  if (not theIntersection)
  {
    return std::nullopt;
  }

  const auto [intPoint, tiOnSegment, tiOnCurve] = theIntersection.value();
  return tiOnCurve;
}

//*********************************************************************************************************************

geom2d::segmentCurve
  geom2d::curveMutualDistanceCalculator::buildSegmentCurveAsNormalToCurve
  (
    const double tOnCurve,
    const baseCurve& curve,
    const double xmin,
    const double xmax,
    const double ymin,
    const double ymax
  )
{
  const auto pc = curve.getPoint(tOnCurve);
  const auto tau = vector{ curve.getVelocity(tOnCurve) }.normalized(); // единичный вектор касательной к кривой в заданной точке
  const vector norm{ tau.y, -tau.x };

  // уравнение отрезка будет таким:
  //   vec{r(t)} = vec{pc} + t x vec{norm}

  const double tofmin = (std::abs(norm.x) > std::abs(norm.y)) ? ((xmin - pc.x) / norm.x) : ((ymin - pc.y) / norm.y);
  const double tofmax = (std::abs(norm.x) > std::abs(norm.y)) ? ((xmax - pc.x) / norm.x) : ((ymax - pc.y) / norm.y);

  const vector vP = vector{ pc } + (tofmin * norm);
  const point P{ vP.x, vP.y };

  const vector vQ = vector{ pc } + (tofmax * norm);
  const point Q{ vQ.x, vQ.y };

  return segmentCurve{ P, Q };
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performScreen1Screen2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  {
    const auto [theDist, theT1, theT2] = execAnyOfSameOrientation(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
    if (theDist < m_dist)
    {
      m_dist = theDist;
      m_t1 = theT1;
      m_t2 = theT2;
    }
  }
  {
    const auto [theDist, theT2, theT1] = execAnyOfSameOrientation(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
    if (theDist < m_dist)
    {
      m_dist = theDist;
      m_t1 = theT1;
      m_t2 = theT2;
    }
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performScreen1Normal2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT1, theT2] = execAnyOfDifferentOrientation(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performScreen1PlatoX2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT2, theT1] = execPlatoX_and_Any(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performScreen1PlatoY2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT2, theT1] = execPlatoY_and_Any(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performScreen1Point2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT2, theT1] = execPoint_and_Any(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performNormal1Screen2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT1, theT2] = execAnyOfDifferentOrientation(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performNormal1Normal2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  {
    const auto [theDist, theT1, theT2] = execAnyOfSameOrientation(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
    if (theDist < m_dist)
    {
      m_dist = theDist;
      m_t1 = theT1;
      m_t2 = theT2;
    }
  }
  {
    const auto [theDist, theT2, theT1] = execAnyOfSameOrientation(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
    if (theDist < m_dist)
    {
      m_dist = theDist;
      m_t1 = theT1;
      m_t2 = theT2;
    }
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performNormal1PlatoX2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT2, theT1] = execPlatoX_and_Any(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performNormal1PlatoY2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT2, theT1] = execPlatoY_and_Any(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performNormal1Point2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT2, theT1] = execPoint_and_Any(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPlatoX1Any2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT1, theT2] = execPlatoX_and_Any(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPlatoX1PlatoX2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  // points of curve 1
  for (const auto t1 : std::initializer_list{ tmin1, tmax1 })
  {
    const auto pnt1 = m_curve1.getPoint(t1);

    using namespace geom2d::extrema::point_and_plato;
    distanceCalculator<DataGetterOfY> calculator{ pnt1, tmin2, tmax2, m_curve2 };

    const auto [dist, ton2] = calculator.execute();
    if (dist < m_dist)
    {
      m_dist = dist;
      m_t1 = t1;
      m_t2 = ton2;
    }
  }
  // points of curve 2
  for (const auto t2 : std::initializer_list{ tmin2, tmax2 })
  {
    const auto pnt2 = m_curve2.getPoint(t2);

    using namespace geom2d::extrema::point_and_plato;
    distanceCalculator<DataGetterOfY> calculator{ pnt2, tmin1, tmax1, m_curve1 };

    const auto [dist, ton1] = calculator.execute();
    if (dist < m_dist)
    {
      m_dist = dist;
      m_t1 = ton1;
      m_t2 = t2;
    }
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPlatoX1PlatoY2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT1, theT2] = execPlatoX_and_PlatoY(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);

  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPlatoX1Point2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT2, theT1] = execPoint_and_PlatoX(tmin2, tmax2, m_curve2, tmin1, tmax1, m_curve1);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************

void
  geom2d::curveMutualDistanceCalculator::performPlatoY1Any2
  (
    const double tmin1,
    const double tmax1,
    const double tmin2,
    const double tmax2
  )
{
  const auto [theDist, theT1, theT2] = execPlatoY_and_Any(tmin1, tmax1, m_curve1, tmin2, tmax2, m_curve2);
  if (theDist < m_dist)
  {
    m_dist = theDist;
    m_t1 = theT1;
    m_t2 = theT2;
  }
}

//*********************************************************************************************************************
