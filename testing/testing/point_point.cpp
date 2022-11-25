#include "pch.h"
#include "pointCurve.h"
#include "curveIntersector.h"

//*********************************************************************************************************************

TEST(Intersector, Point_Point_Not_Intersected)
{
  using namespace geom2d;

  pointCurve point1{ point{0.0, 0.0} };
  pointCurve point2{ point{1.0, 0.0} };

  curveIntersector theIntersector{ point1, point2 };
  theIntersector.perform();

  // theIntersector.dumpIntersections(std::cout);
  const auto points = theIntersector.getSolutionPoints();
  ASSERT_EQ(points.size(), 0) << "two different points intersected!";
}

//*********************************************************************************************************************

TEST(Intersector, Point_Point_Intersected)
{
  using namespace geom2d;

  pointCurve curve1{ point{1.0, 0.0} };
  pointCurve curve2{ point{1.0, 0.0} };

  curveIntersector theIntersector{ curve1, curve2 };
  theIntersector.perform();

  // theIntersector.dumpIntersections(std::cout);
  const auto points = theIntersector.getSolutionPoints();
  ASSERT_EQ(points.size(), 1) << "two equal points do not intersected!";

  const auto solutionIsCorrect = point::isSame(points[0], curve1.thePoint()) and point::isSame(points[0], curve2.thePoint());
  ASSERT_TRUE(solutionIsCorrect);
}

//*********************************************************************************************************************
