#include "pch.h"
#include "segmentCurve.h"

TEST(Segment_Segment_Native_Solver, general_case)
{
  using namespace geom2d;
  segmentCurve seg1{ point{-3.0, -2.3}, point {4.0, 4.0} };
  segmentCurve seg2{ point{-4.0,  2.0}, point {4.1423435, -1.764675} };

  const auto result = segmentCurve::intersectionOf(seg1, seg2);
  ASSERT_TRUE(result);

  const auto [t1, t2] = result.value();
  const auto p1 = seg1.getPoint(t1);
  const auto p2 = seg2.getPoint(t2);

  ASSERT_NEAR(t1, 0.402416, 1.0e-6);
  ASSERT_NEAR(t2, 0.468773, 1.0e-6);

  ASSERT_NEAR(p1.x, -0.183088, 1.0e-6);
  ASSERT_NEAR(p1.y,  0.235221, 1.0e-6);

  ASSERT_NEAR(p2.x, -0.183088, 1.0e-6);
  ASSERT_NEAR(p2.y,  0.235221, 1.0e-6);

  ASSERT_TRUE(point::isSame(p1, p2));
}

TEST(Segment_Segment_Native_Solver, special_cases)
{
  using namespace geom2d;
  // # subcase no. 1
  {
    segmentCurve seg1{ point{1.0, 1.0}, point {2.0, 2.0} };
    segmentCurve seg2{ point{3.0, 3.0}, point {4.0, 4.0} };

    const auto resultOne = segmentCurve::intersectionOf(seg1, seg2);
    ASSERT_TRUE(not resultOne);

    const auto resultTwo = segmentCurve::intersectionOf(seg2, seg1);
    ASSERT_TRUE(not resultOne);
  }
  // # subcase no. 2
  {
    segmentCurve seg1{ point{1.0, 1.0}, point {2.0, 2.0} };
    segmentCurve seg2{ point{3.0, 2.0}, point {4.0, 3.0} };

    const auto resultOne = segmentCurve::intersectionOf(seg1, seg2);
    ASSERT_TRUE(not resultOne);

    const auto resultTwo = segmentCurve::intersectionOf(seg2, seg1);
    ASSERT_TRUE(not resultOne);
  }
  // # subcase no. 3
  {
    segmentCurve seg1{ point{1.0, 1.0}, point {2.0, 2.0} };
    segmentCurve seg2{ point{2.0, 2.0}, point {3.0, 3.0} };
    {
      const auto result = segmentCurve::intersectionOf(seg1, seg2);
      ASSERT_TRUE(result);
      const auto [t1, t2] = result.value();
      ASSERT_NEAR(t1, 1.0, 1.0e-6);
      ASSERT_NEAR(t2, 0.0, 1.0e-6);
    }
    {
      const auto result = segmentCurve::intersectionOf(seg2, seg1);
      ASSERT_TRUE(result);
      const auto [t2, t1] = result.value();
      ASSERT_NEAR(t1, 1.0, 1.0e-6);
      ASSERT_NEAR(t2, 0.0, 1.0e-6);
    }
  }
  // # subcase no. 4
  {
    segmentCurve seg1{ point{1.0, 1.0}, point {2.0, 2.0} };
    segmentCurve seg2{ point{1.5, 1.5}, point {2.5, 2.5} };
    {
      const auto result = segmentCurve::intersectionOf(seg1, seg2);
      ASSERT_TRUE(result);
      const auto [t1, t2] = result.value();
      ASSERT_NEAR(t1, 1.0, 1.0e-6);
      ASSERT_NEAR(t2, 0.5, 1.0e-6);
    }
    {
      const auto result = segmentCurve::intersectionOf(seg2, seg1);
      ASSERT_TRUE(result);
      const auto [t2, t1] = result.value();
      ASSERT_NEAR(t1, 0.5, 1.0e-6);
      ASSERT_NEAR(t2, 0.0, 1.0e-6);
    }
  }
}