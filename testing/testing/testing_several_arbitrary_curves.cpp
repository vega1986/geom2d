#include "pch.h"
#include "arbitraryCurve2d.h"
#include "segmentCurve.h"
#include "curveIntersector.h"

#include <iostream>
#include <iomanip>

// ********************************************************************************************************************

TEST(ArbitraryCurves, parabola_and_segment_01)
{
  auto px = [](const double t) -> double
  {
    return t;
  };
  auto py = [](const double t) -> double
  {
    return t * t;
  };
  auto ux = [](const double t) -> double
  {
    return 1.0;
  };
  auto uy = [](const double t) -> double
  {
    return 2.0 * t;
  };
  using namespace geom2d;
  arbitraryCurve<decltype(px), decltype(py), decltype(ux), decltype(uy)> parabola(px, py, ux, uy, -1.0, 3.0);
  segmentCurve segment{ point{-1.0, -2.0}, point {3.0, 6.0} };

  curveIntersector inter{ parabola, segment };
  inter.fulfill();

  const auto points = inter.getSolutionPoints();
  const auto t1s = inter.getSolutionT1();
  const auto t2s = inter.getSolutionT2();

  ASSERT_EQ(points.size(), 2);

  ASSERT_NEAR(t1s[0], 0.00, 1.0e-6);
  ASSERT_NEAR(t2s[0], 0.25, 1.0e-6);
  ASSERT_NEAR(t1s[1], 2.00, 1.0e-6);
  ASSERT_NEAR(t2s[1], 0.75, 1.0e-6);

  for (size_t j = 0; j < 2; ++j)
  {
    const auto t1 = t1s[j];
    const auto t2 = t2s[j];

    const auto p1 = parabola.getPoint(t1);
    const auto p2 = segment.getPoint(t2);

    ASSERT_TRUE(point::isSame(p1, p2));
  }
}

// ********************************************************************************************************************

TEST(ArbitraryCurves, parabola_and_segment_02)
{
  auto px = [](const double t) -> double
  {
    return t;
  };
  auto py = [](const double t) -> double
  {
    return t * t;
  };
  auto ux = [](const double t) -> double
  {
    return 1.0;
  };
  auto uy = [](const double t) -> double
  {
    return 2.0 * t;
  };
  using namespace geom2d;
  arbitraryCurve<decltype(px), decltype(py), decltype(ux), decltype(uy)> parabola(px, py, ux, uy, -1.0, 3.0);
  segmentCurve segment{ point{-1.0, -2.8}, point {3.0, 5.2} };

  curveIntersector inter{ parabola, segment };
  inter.fulfill();

  const auto points = inter.getSolutionPoints();
  const auto t1s = inter.getSolutionT1();
  const auto t2s = inter.getSolutionT2();

  ASSERT_EQ(points.size(), 2);

  ASSERT_NEAR(t1s[0], 0.552786, 1.0e-6);
  ASSERT_NEAR(t2s[0], 0.388197, 1.0e-6);
  ASSERT_NEAR(t1s[1], 1.447213, 1.0e-6);
  ASSERT_NEAR(t2s[1], 0.611803, 1.0e-6);

  for (size_t j = 0; j < 2; ++j)
  {
    const auto t1 = t1s[j];
    const auto t2 = t2s[j];

    const auto p1 = parabola.getPoint(t1);
    const auto p2 = segment.getPoint(t2);

    ASSERT_TRUE(point::isSame(p1, p2));
  }
}

// ********************************************************************************************************************

TEST(ArbitraryCurves, parabola_and_segment_03)
{
  auto px = [](const double t) -> double
  {
    return t;
  };
  auto py = [](const double t) -> double
  {
    return t * t;
  };
  auto ux = [](const double t) -> double
  {
    return 1.0;
  };
  auto uy = [](const double t) -> double
  {
    return - 2.0 * t;
  };
  using namespace geom2d;
  arbitraryCurve<decltype(px), decltype(py), decltype(ux), decltype(uy)> parabola(px, py, ux, uy, -3.0, 1.0);
  segmentCurve segment{ point{1.0, -2.0}, point {-3.0, 6.0} };

  curveIntersector inter{ parabola, segment };
  inter.fulfill();

  const auto points = inter.getSolutionPoints();
  const auto t1s = inter.getSolutionT1();
  const auto t2s = inter.getSolutionT2();

  ASSERT_EQ(points.size(), 2);

  ASSERT_NEAR(t1s[0], -2.00, 1.0e-6);
  ASSERT_NEAR(t2s[0],  0.75, 1.0e-6);
  ASSERT_NEAR(t1s[1],  0.00, 1.0e-6);
  ASSERT_NEAR(t2s[1],  0.25, 1.0e-6);

  for (size_t j = 0; j < 2; ++j)
  {
    const auto t1 = t1s[j];
    const auto t2 = t2s[j];

    const auto p1 = parabola.getPoint(t1);
    const auto p2 = segment.getPoint(t2);

    ASSERT_TRUE(point::isSame(p1, p2));
  }
}

// ********************************************************************************************************************

TEST(ArbitraryCurves, parabola_and_segment_04)
{
  auto px = [](const double t) -> double
  {
    return t;
  };
  auto py = [](const double t) -> double
  {
    return t * t;
  };
  auto ux = [](const double t) -> double
  {
    return 1.0;
  };
  auto uy = [](const double t) -> double
  {
    return 2.0 * t;
  };
  using namespace geom2d;
  arbitraryCurve<decltype(px), decltype(py), decltype(ux), decltype(uy)> parabola(px, py, ux, uy, -3.0, 1.0);
  segmentCurve segment{ point{1.0, -2.8}, point {-3.0, 5.2} };

  curveIntersector inter{ parabola, segment };
  inter.fulfill();

  const auto points = inter.getSolutionPoints();
  const auto t1s = inter.getSolutionT1();
  const auto t2s = inter.getSolutionT2();

  ASSERT_EQ(points.size(), 2);
  ASSERT_NEAR(t1s[0], -0.552786, 1.0e-6);
  ASSERT_NEAR(t2s[0],  0.388197, 1.0e-6);
  ASSERT_NEAR(t1s[1], -1.447214, 1.0e-6);
  ASSERT_NEAR(t2s[1],  0.611803, 1.0e-6);

  for (size_t j = 0; j < points.size(); ++j)
  {
    const auto t1 = t1s[j];
    const auto t2 = t2s[j];

    const auto p1 = parabola.getPoint(t1);
    const auto p2 = segment.getPoint(t2);

    ASSERT_TRUE(point::isSame(p1, p2));
  }
}

TEST(ArbitraryCurves, parabola_and_segment_as_arbitraries_01)
{
  // parabola
  auto px1 = [](const double t) -> double { return t;       };
  auto py1 = [](const double t) -> double { return t * t;   };
  auto ux1 = [](const double t) -> double { return 1.0;     };
  auto uy1 = [](const double t) -> double { return 2.0 * t; };

  // segment
  auto px2 = [](const double t) -> double { return       t        ; };
  auto py2 = [](const double t) -> double { return 0.1 * t - 0.002; };
  auto ux2 = [](const double t) -> double { return 1.0            ; };
  auto uy2 = [](const double t) -> double { return 0.1            ; };

  using namespace geom2d;
  arbitraryCurve<decltype(px1), decltype(py1), decltype(ux1), decltype(uy1)> parabola(px1, py1, ux1, uy1, -0.005, 0.08);
  arbitraryCurve<decltype(px2), decltype(py2), decltype(ux2), decltype(uy2)>  segment(px2, py2, ux2, uy2, -0.005, 0.08);

  curveIntersector inter{ parabola, segment };
  inter.fulfill();

  const auto points = inter.getSolutionPoints();
  const auto t1s = inter.getSolutionT1();
  const auto t2s = inter.getSolutionT2();

  ASSERT_EQ(points.size(), 2);
  ASSERT_NEAR(t1s[0], 0.027639, 1.0e-6);
  ASSERT_NEAR(t2s[0], 0.027639, 1.0e-6);
  ASSERT_NEAR(t1s[1], 0.072361, 1.0e-6);
  ASSERT_NEAR(t2s[1], 0.072361, 1.0e-6);

  for (size_t j = 0; j < points.size(); ++j)
  {
    const auto t1 = t1s[j];
    const auto t2 = t2s[j];

    const auto p1 = parabola.getPoint(t1);
    const auto p2 = segment.getPoint(t2);

    ASSERT_TRUE(point::isSame(p1, p2));
  }
}

TEST(ArbitraryCurves, parabola_and_segment_as_arbitraries_02)
{
  // parabola
  auto px1 = [](const double t) -> double { return t;       };
  auto py1 = [](const double t) -> double { return t * t;   };
  auto ux1 = [](const double t) -> double { return 1.0;     };
  auto uy1 = [](const double t) -> double { return 2.0 * t; };

  // segment
  auto px2 = [](const double t) -> double { return         t; };
  auto py2 = [](const double t) -> double { return - 0.1 * t - 0.002; };
  auto ux2 = [](const double t) -> double { return   1.0; };
  auto uy2 = [](const double t) -> double { return - 0.1; };

  using namespace geom2d;
  arbitraryCurve<decltype(px1), decltype(py1), decltype(ux1), decltype(uy1)> parabola(px1, py1, ux1, uy1, -0.08, 0.005);
  arbitraryCurve<decltype(px2), decltype(py2), decltype(ux2), decltype(uy2)>  segment(px2, py2, ux2, uy2, -0.08, 0.005);

  curveIntersector inter{ parabola, segment };
  inter.fulfill();

  const auto points = inter.getSolutionPoints();
  const auto t1s = inter.getSolutionT1();
  const auto t2s = inter.getSolutionT2();

  ASSERT_EQ(points.size(), 2);
  ASSERT_NEAR(t1s[0], -0.072361, 1.0e-6);
  ASSERT_NEAR(t2s[0], -0.072361, 1.0e-6);
  ASSERT_NEAR(t1s[1], -0.027639, 1.0e-6);
  ASSERT_NEAR(t2s[1], -0.027639, 1.0e-6);

  for (size_t j = 0; j < points.size(); ++j)
  {
    const auto t1 = t1s[j];
    const auto t2 = t2s[j];

    const auto p1 = parabola.getPoint(t1);
    const auto p2 = segment.getPoint(t2);

    ASSERT_TRUE(point::isSame(p1, p2));
  }
}
