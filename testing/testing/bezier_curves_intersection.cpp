#include "pch.h"
#include "bezierCurve2d.h"
#include "curveIntersector.h"

TEST(BezierCurvesIntersection, test_01)
{
  using namespace geom2d;

  bezierCurve bez1(std::initializer_list<point>{ {0.0, 0.0}, { 1.0, 0.0 }, { 1.0, 1.0 }});
  bezierCurve bez2(std::initializer_list<point>{ {0.0, 1.0}, { 0.0, 0.0 }, { 1.0, 0.0 }});

  curveIntersector inter{ bez1, bez2 };
  inter.fulfill();

  const auto pnts = inter.getSolutionPoints();
  const auto t1s = inter.getSolutionT1();
  const auto t2s = inter.getSolutionT2();

  ASSERT_EQ(pnts.size(), 1);
  
  for (size_t j = 0; j < pnts.size(); ++j)
  {
    const auto t1 = t1s[j];
    const auto t2 = t2s[j];

    const auto p1 = bez1.getPoint(t1);
    const auto p2 = bez2.getPoint(t2);

    ASSERT_TRUE(point::isSame(p1, p2));
  }
  
  ASSERT_NEAR(pnts.front().x, 0.5000000, 1.0e-7);
  ASSERT_NEAR(pnts.front().y, 0.0857864, 1.0e-7);
}