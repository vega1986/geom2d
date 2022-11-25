#include "pch.h"
#include "pointCurve.h"
#include "segmentCurve.h"
#include "curveIntersector.h"

#include <iostream>
#include <array>

//*********************************************************************************************************************

TEST(Intersector, Point_Segment)
{
  using namespace geom2d;
  constexpr double tabserror = 1.0e-9;
  
  auto segmentTesting = [](const segmentCurve& segment)
  {
    const auto t1 = segment.parameterMin();
    const auto t2 = segment.parameterMax();

    const auto p1 = segment.getPoint(t1);
    const auto p2 = segment.getPoint(t2);

    auto pointOnSegment = [p1, p2](const double t)->point
    {
      return (p1 * (1.0 - t)) + (p2 * t);
    };
    

    std::array<point, 5> allPoints
    {
      pointOnSegment(-0.5), // not intersected
      pointOnSegment( 0.0), // intersected | t of segment = 0.0
      pointOnSegment( 0.5), // intersected | t of segment = 0.5
      pointOnSegment( 1.0), // intersected | t of segment = 1.0
      pointOnSegment( 1.5)  // not intersected
    };

    // THIS IS the referent parameters !!!!!!
    constexpr std::array<const size_t, 5> countofIntersections{ 0, 1, 1, 1, 0 };
    constexpr std::array<const double, 5> tofSegment{-0.5, 0.0, 0.5, 1.0, 1.5};

    for (size_t j = 0; j < 5; ++j)
    {
      pointCurve pntCurve{ allPoints[j]};
      curveIntersector theIntersector{ segment , pntCurve };
      theIntersector.perform();

      const auto allIntersections = theIntersector.getSolutionPoints();
      const auto allTofsegment = theIntersector.getSolutionT1();

      ASSERT_EQ(countofIntersections[j], allIntersections.size());
      if (allIntersections.size() > 0)
      {
        ASSERT_NEAR(tofSegment[j], allTofsegment[0], 1.0e-9);
      }
    }
  };
  
  std::array<segmentCurve, 12> segCurves
  {
    segmentCurve{point{ 1.0,  0.0}, point{ 2.0,  0.0}},
    segmentCurve{point{ 2.0,  0.0}, point{ 1.0,  0.0}}, // inversed

    segmentCurve{point{-1.0,  0.0}, point{-2.0,  0.0}},
    segmentCurve{point{-2.0,  0.0}, point{-1.0,  0.0}}, // inversed

    segmentCurve{point{ 0.0,  1.0}, point{ 0.0,  2.0}},
    segmentCurve{point{ 0.0,  2.0}, point{ 0.0,  1.0}}, // inversed

    segmentCurve{point{ 0.0, -1.0}, point{ 0.0, -2.0}},
    segmentCurve{point{ 0.0, -2.0}, point{ 0.0, -1.0}}, // inversed

    segmentCurve{point{-1.0, -1.0}, point{ 1.0,  1.0}},
    segmentCurve{point{ 1.0,  1.0}, point{-1.0, -1.0}}, // inversed

    segmentCurve{point{ 1.0, -1.0}, point{-1.0,  1.0}},
    segmentCurve{point{-1.0,  1.0}, point{ 1.0, -1.0}}  // inversed
  };
  for (const auto& segCurve : segCurves)
  {
    segmentTesting(segCurve);
  }
}