#include "pch.h"
#include "segmentCurve.h"
#include "curveIntersector.h"

#include <iostream>
#include <array>

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_hor_01)
{
  using namespace geom2d;
  
  {
    const double y = 0.0;
    // coordinates of first segment
    const double x11 = 0.0;
    const double x12 = 1.0;
    // coordinates of second segment
    const double x21 = 2.0;
    const double x22 = 3.0;

    segmentCurve seg1{ point{x11, y}, point{x12, y} };
    segmentCurve seg2{ point{x21, y}, point{x22, y} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.perform();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 0);
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_hor_02)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;

  {
    const double y = 0.0;
    // coordinates of first segment
    const double x11 = 0.0;
    const double x12 = 1.0;
    // coordinates of second segment
    const double x21 = 1.0;
    const double x22 = 2.0;

    segmentCurve seg1{ point{x11, y}, point{x12, y} };
    segmentCurve seg2{ point{x21, y}, point{x22, y} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.perform();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 1);
    ASSERT_TRUE(point::isSame(point{ x12, y }, points[0]));
    ASSERT_NEAR(t1s[0], 1.0, tabserror);
    ASSERT_NEAR(t2s[0], 0.0, tabserror);
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_hor_03)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;

  {
    const double y = 0.0;
    // coordinates of first segment
    const double x11 = 0.0;
    const double x12 = 1.0;
    // coordinates of second segment
    const double x21 = 0.5;
    const double x22 = 1.5;

    segmentCurve seg1{ point{x11, y}, point{x12, y} };
    segmentCurve seg2{ point{x21, y}, point{x22, y} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.perform();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 1);

    const auto pnt = points[0];
    const auto t1 = t1s[0];
    const auto t2 = t2s[0];

    ASSERT_TRUE(point::isSame(pnt, point{ 1.0, y }));
    ASSERT_NEAR(t1, 1.0, tabserror);
    ASSERT_NEAR(t2, 0.5, tabserror);
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_hor_04)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;

  {
    const double y = 0.0;
    // coordinates of first segment
    const double x11 = 0.0;
    const double x12 = 1.0;
    // coordinates of second segment
    const double x21 = -0.5;
    const double x22 = 0.5;

    segmentCurve seg1{ point{x11, y}, point{x12, y} };
    segmentCurve seg2{ point{x21, y}, point{x22, y} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.perform();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 1);

    const auto pnt = points[0];
    const auto t1 = t1s[0];
    const auto t2 = t2s[0];

    ASSERT_TRUE(point::isSame(pnt, point{ 0.0, y }));
    ASSERT_NEAR(t1, 0.0, tabserror);
    ASSERT_NEAR(t2, 0.5, tabserror);
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_hor_05)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;

  {
    const double y = 0.0;
    // coordinates of first segment
    const double x11 = 0.0;
    const double x12 = 1.0;
    // coordinates of second segment
    const double x21 = -1.0;
    const double x22 = 0.0;

    segmentCurve seg1{ point{x11, y}, point{x12, y} };
    segmentCurve seg2{ point{x21, y}, point{x22, y} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.perform();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 1);

    const auto pnt = points[0];
    const auto t1 = t1s[0];
    const auto t2 = t2s[0];

    ASSERT_TRUE(point::isSame(pnt, point{ 0.0, y }));
    ASSERT_NEAR(t1, 0.0, tabserror);
    ASSERT_NEAR(t2, 1.0, tabserror);
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_hor_06)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;

  {
    const double y = 0.0;
    // coordinates of first segment
    const double x11 = 0.0;
    const double x12 = 1.0;
    // coordinates of second segment
    const double x21 = -2.0;
    const double x22 = -1.0;

    segmentCurve seg1{ point{x11, y}, point{x12, y} };
    segmentCurve seg2{ point{x21, y}, point{x22, y} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.perform();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 0);
  }
}