#include "pch.h"
#include "segmentCurve.h"
#include "curveIntersector.h"

#include <iostream>
#include <array>

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_01)
{
  using namespace geom2d;

  {
    const double x = 0.0;
    // coordinates of first segment
    const double y11 = 0.0;
    const double y12 = 1.0;
    // coordinates of second segment
    const double y21 = 2.0;
    const double y22 = 3.0;

    segmentCurve seg1{ point{x, y11}, point{x, y12} };
    segmentCurve seg2{ point{x, y21}, point{x, y22} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.fulfill();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 0);
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_02)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;

  {
    const double x = 0.0;
    // coordinates of first segment
    const double y11 = 0.0;
    const double y12 = 1.0;
    // coordinates of second segment
    const double y21 = 1.0;
    const double y22 = 2.0;

    segmentCurve seg1{ point{x, y11}, point{x, y12} };
    segmentCurve seg2{ point{x, y21}, point{x, y22} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.fulfill();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 1);
    ASSERT_TRUE(point::isSame(point{ x, y12 }, points[0]));
    ASSERT_NEAR(t1s[0], 1.0, tabserror);
    ASSERT_NEAR(t2s[0], 0.0, tabserror);
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_03)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;

  {
    const double x = 0.0;
    // coordinates of first segment
    const double y11 = 0.0;
    const double y12 = 1.0;
    // coordinates of second segment
    const double y21 = 0.5;
    const double y22 = 1.5;

    segmentCurve seg1{ point{x, y11}, point{x, y12} };
    segmentCurve seg2{ point{x, y21}, point{x, y22} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.fulfill();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 1);

    const auto pnt = points[0];
    const auto t1 = t1s[0];
    const auto t2 = t2s[0];

    ASSERT_TRUE(point::isSame(pnt, point{ x, 1.0 }));
    ASSERT_NEAR(t1, 1.0, tabserror);
    ASSERT_NEAR(t2, 0.5, tabserror);
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_04)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;

  {
    const double x = 0.0;
    // coordinates of first segment
    const double y11 = 0.0;
    const double y12 = 1.0;
    // coordinates of second segment
    const double y21 = -0.5;
    const double y22 = 0.5;

    segmentCurve seg1{ point{x, y11}, point{x, y12} };
    segmentCurve seg2{ point{x, y21}, point{x, y22} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.fulfill();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 1);

    const auto pnt = points[0];
    const auto t1 = t1s[0];
    const auto t2 = t2s[0];

    ASSERT_TRUE(point::isSame(pnt, point{ x, 0.0 }));
    ASSERT_NEAR(t1, 0.0, tabserror);
    ASSERT_NEAR(t2, 0.5, tabserror);
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_05)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;

  {
    const double x = 0.0;
    // coordinates of first segment
    const double y11 = 0.0;
    const double y12 = 1.0;
    // coordinates of second segment
    const double y21 = -1.0;
    const double y22 = 0.0;

    segmentCurve seg1{ point{x, y11}, point{x, y12} };
    segmentCurve seg2{ point{x, y21}, point{x, y22} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.fulfill();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 1);

    const auto pnt = points[0];
    const auto t1 = t1s[0];
    const auto t2 = t2s[0];

    ASSERT_TRUE(point::isSame(pnt, point{ x, 0.0 }));
    ASSERT_NEAR(t1, 0.0, tabserror);
    ASSERT_NEAR(t2, 1.0, tabserror);
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_06)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;

  {
    const double x = 0.0;
    // coordinates of first segment
    const double y11 = 0.0;
    const double y12 = 1.0;
    // coordinates of second segment
    const double y21 = -2.0;
    const double y22 = -1.0;

    segmentCurve seg1{ point{x, y11}, point{x, y12} };
    segmentCurve seg2{ point{x, y21}, point{x, y22} };

    curveIntersector theIntersector{ seg1, seg2 };
    theIntersector.fulfill();

    const auto points = theIntersector.getSolutionPoints();
    const auto t1s = theIntersector.getSolutionT1();
    const auto t2s = theIntersector.getSolutionT2();

    ASSERT_EQ(points.size(), 0);
  }
}