#include "pch.h"
#include "curveMutualDistanceCalculator.h"
#include "pointCurve.h"
#include "segmentCurve.h"
#include "arbitraryCurve2d.h"
#include "bezierCurve2d.h"

#include <random>
#include <cmath>
#include <algorithm>
#include <numbers>

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, Point_and_Point)
{
  using namespace geom2d;
  pointCurve curPnt1{ point{ 3.0, 1.0 } };
  pointCurve curPnt2{ point{ 6.0, 5.0 } };
  curveMutualDistanceCalculator extremaCalculator{ curPnt1, curPnt2 };
  extremaCalculator.fulfill();
  const auto [dist, t1, t2] = extremaCalculator.getExtrema();

  ASSERT_NEAR(dist, 5.0, 1.0e-9);
};

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, Point_and_PlatoY)
{
  using namespace geom2d;
  segmentCurve seg{ point{1.0, 2.0}, point{5.0, 2.0} };

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-6.0, 6.0);

  constexpr unsigned short n = 798;
  for (unsigned short j = 0; j < n; ++j)
  {
    pointCurve pntCurve{ point{dis(gen), dis(gen)} };

    const auto distNativeSolution = seg.distanceTo(pntCurve.thePoint());
    const auto t_extremaNativeSolution = seg.nearestTo(pntCurve.thePoint());

    if (distNativeSolution < 1.0e-3)
    {
      continue;
    }

    {
      curveMutualDistanceCalculator calculator{ pntCurve, seg };
      calculator.fulfill();

      const auto [distNumericalSolution, tpnt, t_extremaNumericalSolution] = calculator.getExtrema();

      ASSERT_NEAR(distNativeSolution, distNumericalSolution, 1.0e-9);
      ASSERT_NEAR(t_extremaNativeSolution, t_extremaNumericalSolution, 1.0e-7);
    }
    {
      curveMutualDistanceCalculator calculator{ seg, pntCurve };
      calculator.fulfill();

      const auto [distNumericalSolution, t_extremaNumericalSolution, tpnt ] = calculator.getExtrema();

      ASSERT_NEAR(distNativeSolution, distNumericalSolution, 1.0e-9);
      ASSERT_NEAR(t_extremaNativeSolution, t_extremaNumericalSolution, 1.0e-7);
    }

  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, Point_and_PlatoX)
{
  using namespace geom2d;
  segmentCurve seg{ point{2.0, 1.0}, point{2.0, 5.0} };

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-6.0, 6.0);

  constexpr unsigned short n = 987;
  for (unsigned short j = 0; j < n; ++j)
  {
    pointCurve pntCurve{ point{dis(gen), dis(gen)} };

    const auto distNativeSolution = seg.distanceTo(pntCurve.thePoint());
    const auto t_extremaNativeSolution = seg.nearestTo(pntCurve.thePoint());

    if (distNativeSolution < 1.0e-3)
    {
      continue;
    }

    {
      curveMutualDistanceCalculator calculator{ pntCurve, seg };
      calculator.fulfill();

      const auto [distNumericalSolution, tpnt, t_extremaNumericalSolution] = calculator.getExtrema();

      ASSERT_NEAR(distNativeSolution, distNumericalSolution, 1.0e-9);
      ASSERT_NEAR(t_extremaNativeSolution, t_extremaNumericalSolution, 1.0e-7);
    }
    {
      curveMutualDistanceCalculator calculator{ seg, pntCurve };
      calculator.fulfill();

      const auto [distNumericalSolution, t_extremaNumericalSolution, tpnt ] = calculator.getExtrema();

      ASSERT_NEAR(distNativeSolution, distNumericalSolution, 1.0e-9);
      ASSERT_NEAR(t_extremaNativeSolution, t_extremaNumericalSolution, 1.0e-7);
    }
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, Point_and_Arbitrary_Segment)
{
  using namespace geom2d;

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis4point(-6.0, 6.0);
  std::uniform_real_distribution<> dis4segment(-5.0, 5.0);

  constexpr size_t n = 19855;
  size_t passed = 0;
  while (passed < n)
  {
    segmentCurve seg{ point{dis4segment(gen), dis4segment(gen)}, point{dis4segment(gen), dis4segment(gen)}};
    {
      const auto length = seg.length();
      if (length < 0.01)
      {
        continue;
      }

      const auto tau = seg.getTau();
      if ((std::abs(tau.x) < 1.0e-2) or (std::abs(tau.y) < 1.0e-2))
      {
        continue;
      }
    }
    
    pointCurve pntCurve{ point{dis4point(gen), dis4point(gen)} };

    const auto distNativeSolution = seg.distanceTo(pntCurve.thePoint());
    const auto t_extremaNativeSolution = seg.nearestTo(pntCurve.thePoint());

    if (distNativeSolution < 1.0e-3)
    {
      continue;
    }

    {
      curveMutualDistanceCalculator calculator{ pntCurve, seg };
      calculator.fulfill();

      const auto [distNumericalSolution, tpnt, t_extremaNumericalSolution] = calculator.getExtrema();

      ASSERT_NEAR(distNativeSolution, distNumericalSolution, 1.0e-9);
      ASSERT_NEAR(t_extremaNativeSolution, t_extremaNumericalSolution, 1.0e-7);
    }
    {
      curveMutualDistanceCalculator calculator{ seg, pntCurve };
      calculator.fulfill();

      const auto [distNumericalSolution, t_extremaNumericalSolution, tpnt] = calculator.getExtrema();

      ASSERT_NEAR(distNativeSolution, distNumericalSolution, 1.0e-9);
      ASSERT_NEAR(t_extremaNativeSolution, t_extremaNumericalSolution, 1.0e-7);
    }
    ++passed;
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, PlatoX_and_PlatoX_1)
{
  using namespace geom2d;

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-6.0, 6.0);

  constexpr size_t n = 8544;
  size_t passed = 0;

  while (passed < n)
  {
    const double xof1 = dis(gen);
    const double yP1 = dis(gen);
    const double yQ1 = dis(gen);
    const double ymin1 = std::min(yP1, yQ1);
    const double ymax1 = std::max(yP1, yQ1);

    const double xof2 = dis(gen);
    const double yP2 = dis(gen);
    const double yQ2 = dis(gen);
    const double ymin2 = std::min(yP2, yQ2);
    const double ymax2 = std::max(yP2, yQ2);

    if (not (ymax1 < ymin2))
    {
      continue;
    }

    if ((ymax1 - ymin1 <= 0.001) or (ymax2 - ymin2 <= 0.001))
    {
      continue;
    }
    
    segmentCurve seg1{ point{xof1, ymin1}, point{xof1, ymax1} };
    segmentCurve seg2{ point{xof2, ymin2}, point{xof2, ymax2} };

    const auto [dist, tt1, tt2] = segmentCurve::distanceBetween(seg1, seg2);

    if (dist <= 0.001)
    {
      continue;
    }

    curveMutualDistanceCalculator calculator{ seg1, seg2 };
    calculator.fulfill();

    const auto [dnum, t1num, t2num] = calculator.getExtrema();

    ASSERT_NEAR(dnum, dist, 1.0e-9);
    ASSERT_NEAR(t1num, 1.0, 1.0e-7);
    ASSERT_NEAR(t2num, 0.0, 1.0e-7);

    ++passed;
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, PlatoX_and_PlatoX_2)
{
  // здесь мы проверяем только dist так как в данной задаче есть неоднозначность в значении параметров решения
  using namespace geom2d;

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-6.0, 6.0);

  constexpr size_t n = 29854;
  size_t passed = 0;

  while (passed < n)
  {
    const double xof1 = dis(gen);
    const double yP1 = dis(gen);
    const double yQ1 = dis(gen);
    const double ymin1 = std::min(yP1, yQ1);
    const double ymax1 = std::max(yP1, yQ1);

    const double xof2 = dis(gen);
    const double yP2 = dis(gen);
    const double yQ2 = dis(gen);
    const double ymin2 = std::min(yP2, yQ2);
    const double ymax2 = std::max(yP2, yQ2);

    if ((ymax1 - ymin1 <= 2.0) or (ymax2 - ymin2 <= 2.0))
    {
      continue;
    }

    segmentCurve seg1{ point{xof1, ymin1}, point{xof1, ymax1} };
    segmentCurve seg2{ point{xof2, ymin2}, point{xof2, ymax2} };

    const auto [dist, tt1, tt2] = segmentCurve::distanceBetween(seg1, seg2);

    if (dist <= 0.001)
    {
      continue;
    }

    curveMutualDistanceCalculator calculator{ seg1, seg2 };
    calculator.fulfill();

    const auto [dnum, t1num, t2num] = calculator.getExtrema();

    ASSERT_NEAR(dnum, dist, 1.0e-9);

    ++passed;
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, PlatoY_and_PlatoY_1)
{
  using namespace geom2d;

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-6.0, 6.0);

  constexpr size_t n = 9865;
  size_t passed = 0;

  while (passed < n)
  {
    const double yof1 = dis(gen);
    const double xP1 = dis(gen);
    const double xQ1 = dis(gen);
    const double xmin1 = std::min(xP1, xQ1);
    const double xmax1 = std::max(xP1, xQ1);

    const double yof2 = dis(gen);
    const double xP2 = dis(gen);
    const double xQ2 = dis(gen);
    const double xmin2 = std::min(xP2, xQ2);
    const double xmax2 = std::max(xP2, xQ2);

    if (not (xmax1 < xmin2))
    {
      continue;
    }

    if ((xmax1 - xmin1 <= 0.001) or (xmax2 - xmin2 <= 0.001))
    {
      continue;
    }

    segmentCurve seg1{ point{xmin1, yof1}, point{xmax1, yof1} };
    segmentCurve seg2{ point{xmin2, yof2}, point{xmax2, yof2} };

    const auto [dist, tt1, tt2] = segmentCurve::distanceBetween(seg1, seg2);

    if (dist <= 0.001)
    {
      continue;
    }

    curveMutualDistanceCalculator calculator{ seg1, seg2 };
    calculator.fulfill();

    const auto [dnum, t1num, t2num] = calculator.getExtrema();

    ASSERT_NEAR(dnum, dist, 1.0e-9);
    ASSERT_NEAR(t1num, 1.0, 1.0e-7);
    ASSERT_NEAR(t2num, 0.0, 1.0e-7);

    ++passed;
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, PlatoY_and_PlatoY_2)
{
  // здесь мы проверяем только dist так как в данной задаче есть неоднозначность в значении параметров решения
  using namespace geom2d;

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-6.0, 6.0);

  constexpr size_t n = 17896;
  size_t passed = 0;

  while (passed < n)
  {
    const double yof1 = dis(gen);
    const double xP1 = dis(gen);
    const double xQ1 = dis(gen);
    const double xmin1 = std::min(xP1, xQ1);
    const double xmax1 = std::max(xP1, xQ1);

    const double yof2 = dis(gen);
    const double xP2 = dis(gen);
    const double xQ2 = dis(gen);
    const double xmin2 = std::min(xP2, xQ2);
    const double xmax2 = std::max(xP2, xQ2);

    if ((xmax1 - xmin1 <= 2.0) or (xmax2 - xmin2 <= 2.0))
    {
      continue;
    }

    segmentCurve seg1{ point{xmin1, yof1}, point{xmax1, yof1} };
    segmentCurve seg2{ point{xmin2, yof2}, point{xmax2, yof2} };

    const auto [dist, tt1, tt2] = segmentCurve::distanceBetween(seg1, seg2);

    if (dist <= 0.001)
    {
      continue;
    }

    curveMutualDistanceCalculator calculator{ seg1, seg2 };
    calculator.fulfill();

    const auto [dnum, t1num, t2num] = calculator.getExtrema();

    ASSERT_NEAR(dnum, dist, 1.0e-9);

    ++passed;
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, PlatoX_and_PlatoY)
{
  using namespace geom2d;

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-6.0, 6.0);

  constexpr size_t n = 19875;
  size_t passed = 0;

  while (passed < n)
  {
    // seg1 - PlatoX
    const double xof1 = dis(gen);
    segmentCurve seg1{ point{xof1, dis(gen)}, point{xof1, dis(gen)} };
    const double yof2 = dis(gen);
    segmentCurve seg2{ point{dis(gen), yof2}, point{dis(gen), yof2}};

    if ((seg1.length() < 1.0e-3) or (seg2.length() < 1.0e-3))
    {
      continue;
    }

    const auto [dist, tt1, tt2] = segmentCurve::distanceBetween(seg1, seg2);

    if (dist <= 1.0e-2)
    {
      continue;
    }

    curveMutualDistanceCalculator calculator{ seg1, seg2 };
    calculator.fulfill();

    const auto [dnum, t1num, t2num] = calculator.getExtrema();

    ASSERT_NEAR(dnum, dist, 1.0e-9);
    ASSERT_NEAR(t1num, tt1, 1.0e-7);
    ASSERT_NEAR(t2num, tt2, 1.0e-7);

    ++passed;
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, PlatoX_Any_as_segment)
{
  // Здесь мы моделируем произвольную кривую Any отрезком, поэтому в некотором смысле это не полноценное тестирование
  using namespace geom2d;

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-6.0, 6.0);

  constexpr size_t n = 95784;
  size_t passed = 0;

  while (passed < n)
  {
    const auto xof1 = dis(gen);
    segmentCurve seg_PlatoX{ point{xof1, dis(gen)}, point {xof1, dis(gen)} };
    if (seg_PlatoX.length() < 1.0e-2)
    {
      continue;
    }

    // теперь получаем другой сегмент моделирующий кривую класса не PlatoX
    segmentCurve seg_Any{ point{dis(gen), dis(gen)}, point{dis(gen), dis(gen)} };
    if (seg_Any.length() < 1.0e-2)
    {
      continue;
    }

    const auto tauOfAny = seg_Any.getTau();

    if (std::abs(tauOfAny.x) < 1.0e-2)
    {
      continue;
    }

    const auto [dist, tseg, tany] = segmentCurve::distanceBetween(seg_PlatoX, seg_Any);

    if (dist <= 1.0e-2)
    {
      continue;
    }

    curveMutualDistanceCalculator calculator{ seg_PlatoX, seg_Any };
    calculator.fulfill();

    const auto [dnum, t_PlatoX, t_Any] = calculator.getExtrema();

    ASSERT_NEAR(dist, dnum    , 1.0e-9);
    ASSERT_NEAR(tseg, t_PlatoX, 1.0e-7);
    ASSERT_NEAR(tany, t_Any   , 1.0e-7);

    ++passed;
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, PlatoY_Any_as_segment)
{
  // Здесь мы моделируем произвольную кривую Any отрезком, поэтому в некотором смысле это не полноценное тестирование
  using namespace geom2d;

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-6.0, 6.0);

  constexpr size_t n = 78965;
  size_t passed = 0;

  while (passed < n)
  {
    const auto yof1 = dis(gen);
    segmentCurve seg_PlatoY{ point{dis(gen), yof1}, point {dis(gen), yof1} };
    if (seg_PlatoY.length() < 1.0e-2)
    {
      continue;
    }

    // теперь получаем другой сегмент моделирующий кривую класса не PlatoX
    segmentCurve seg_Any{ point{dis(gen), dis(gen)}, point{dis(gen), dis(gen)} };
    if (seg_Any.length() < 1.0e-2)
    {
      continue;
    }

    const auto tauOfAny = seg_Any.getTau();

    if (std::abs(tauOfAny.y) < 1.0e-2)
    {
      continue;
    }

    const auto [dist, tseg, tany] = segmentCurve::distanceBetween(seg_PlatoY, seg_Any);

    if (dist <= 1.0e-2)
    {
      continue;
    }

    curveMutualDistanceCalculator calculator{ seg_PlatoY, seg_Any };
    calculator.fulfill();

    const auto [dnum, t_PlatoY, t_Any] = calculator.getExtrema();

    ASSERT_NEAR(dist, dnum, 1.0e-9);
    ASSERT_NEAR(tseg, t_PlatoY, 1.0e-7);
    ASSERT_NEAR(tany, t_Any, 1.0e-7);

    ++passed;
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, Any_Segments)
{
  // Здесь мы моделируем произвольную кривую Any отрезком, поэтому в некотором смысле это не полноценное тестирование
  using namespace geom2d;

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-6.0, 6.0);

  constexpr size_t n = 15897;
  size_t passed = 0;

  while (passed < n)
  {
    segmentCurve seg1{ point{dis(gen), dis(gen)}, point {dis(gen), dis(gen)} };
    if (seg1.length() < 1.0)
    {
      continue;
    }

    segmentCurve seg2{ point{dis(gen), dis(gen)}, point{dis(gen), dis(gen)} };
    if (seg2.length() < 1.0)
    {
      continue;
    }

    const auto tau1 = seg1.getTau();
    const auto tau2 = seg2.getTau();

    if ((tau1, tau2) > 0.99)
    {
      continue;
    }

    const auto [dist, t1, t2] = segmentCurve::distanceBetween(seg1, seg2);

    if (dist <= 1.0e-2)
    {
      continue;
    }

    curveMutualDistanceCalculator calculator{ seg1, seg2 };
    calculator.fulfill();

    const auto [dnum, tnum1, tnum2] = calculator.getExtrema();

    ASSERT_NEAR(dist, dnum, 1.0e-9);
    ASSERT_NEAR(t1, tnum1, 1.0e-7);
    ASSERT_NEAR(t2, tnum2, 1.0e-7);

    ++passed;
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, Point_and_Circle)
{
  auto px = [](const double fi) -> double
  {
    return std::cos(fi);
  };

  auto py = [](const double fi) -> double
  {
    return std::sin(fi);
  };

  auto ux = [](const double fi) -> double
  {
    return -std::sin(fi);
  };

  auto uy = [](const double fi) -> double
  {
    return std::cos(fi);
  };

  using namespace geom2d;
  using namespace std::numbers;
  arbitraryCurve circ(px, py, ux, uy, -pi, pi);

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-3.0, 3.0);

  constexpr size_t n = 75896;
  size_t passed = 0;

  while (passed < n)
  {
    pointCurve pntCurve{ point{dis(gen), dis(gen)}};
    const double R = point::distance(point{ 0.0, 0.0 }, pntCurve.thePoint());
    if (R < 1.01)
    {
      continue;
    }

    curveMutualDistanceCalculator calc(pntCurve, circ);
    calc.fulfill();
    const auto [dist, tofPnt, tofCirc] = calc.getExtrema();

    ASSERT_NEAR(dist, R - 1.0, 1.0e-9);

    double trueTofCirc = std::atan2(pntCurve.thePoint().y, pntCurve.thePoint().x);
    ASSERT_NEAR(tofCirc, trueTofCirc, 1.0e-5);
    ++passed;    
  }
}

// oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

TEST(EXTREMA, Two_Circles)
{
  using namespace geom2d;
  using namespace std::numbers;
  point pc1{0.0, 0.0};
  double R1 = 1.0;

  point pc2{3.0, 4.0};
  double R2 = 2.0;

  auto px1 = [&pc1, &R1](const double fi) -> double { return pc1.x + R1 * std::cos(fi); };
  auto py1 = [&pc1, &R1](const double fi) -> double { return pc1.y + R1 * std::sin(fi); };
  auto ux1 = [&R1](const double fi) -> double { return -R1 * std::sin(fi); };
  auto uy1 = [&R1](const double fi) -> double { return R1 * std::cos(fi); };
  arbitraryCurve circ1(px1, py1, ux1, uy1, -pi, pi);

  auto px2 = [&pc2, &R2](const double fi) -> double { return pc2.x + R2 * std::cos(fi); };
  auto py2 = [&pc2, &R2](const double fi) -> double { return pc2.y + R2 * std::sin(fi); };
  auto ux2 = [&R2](const double fi) -> double { return -R2 * std::sin(fi); };
  auto uy2 = [&R2](const double fi) -> double { return R2 * std::cos(fi); };
  arbitraryCurve circ2(px2, py2, ux2, uy2, -pi, pi);

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> disPos(-5.0, 5.0);
  std::uniform_real_distribution<> disRad(0.5, 4.0);

  constexpr size_t n = 1024;
  size_t passed = 0;
  while (passed < n)
  {
    pc1.x = /*3.11506;*/ disPos(gen);
    pc1.y = /*-4.97108;*/ disPos(gen);
    R1 = /*3.0162;*/ disRad(gen);

    pc2.x = /*2.60125;*/ disPos(gen);
    pc2.y = /*4.55262;*/ disPos(gen);
    R2 = /*3.33606;*/ disRad(gen);

    const double NativeDistance = point::distance(pc1, pc2) - R1 - R2;
    

    if (NativeDistance < 0.1)
    {
      continue;
    }

    curveMutualDistanceCalculator calc(circ1, circ2);
    calc.fulfill();

    const auto [dist, t1, t2] = calc.getExtrema();

    ASSERT_NEAR(dist, NativeDistance, 1.0e-4);
    // std::cout << "passed: " << passed << std::endl;
    ++passed;
  }
}

TEST(EXTREMA, BEZIER_AND_BEZIER)
{
  using namespace geom2d;
  
  bezierCurve bc1{ std::initializer_list{point{-1.0, 0.0}, point{0.0, 0.0}, point{0.0, 1.0}} };
  bezierCurve bc2{ std::initializer_list{point{0.0, -1.0}, point{0.0, 0.0}, point{1.0, 0.0}} };

  curveMutualDistanceCalculator calc{ bc1, bc2 };
  calc.fulfill();

  const auto [dist, t1, t2] = calc.getExtrema();

  ASSERT_NEAR(t1, 0.5, 1.0e-10);
  ASSERT_NEAR(t2, 0.5, 1.0e-10);
}