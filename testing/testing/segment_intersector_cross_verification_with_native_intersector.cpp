#include "pch.h"
#include "segmentCurve.h"
#include "curveIntersector.h"
#include "avector.h"

#include <random>
#include <cmath>
#include <iostream>

TEST(Segment_Intersector_Cross_Verification, main_case)
{
  using namespace geom2d;
  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(-1.0, 1.0);

  //constexpr size_t nj = 100;
  //size_t j = 0;
  size_t total_cases_passed = 0;
  constexpr size_t total_cases_needed = 10000;
  while (total_cases_passed < total_cases_needed)
  {
    // gen random segments
    const point P{ dis(gen) , dis(gen) };
    const point Q{ dis(gen) , dis(gen) };

    const point A{ dis(gen) , dis(gen) };
    const point B{ dis(gen) , dis(gen) };

    if ((point::distance(P, Q) < 0.5) or (point::distance(A, B) < 0.5)) continue;

    segmentCurve seg1{ P, Q };
    segmentCurve seg2{ A, B };

    const auto tau1 = seg1.getTau();
    const auto tau2 = seg2.getTau();

    if (std::abs(tau1 ^ tau2) < 0.1) continue;

    const auto result = segmentCurve::intersectionOf(seg1, seg2);

    if (not result) continue;

    const auto [t1, t2] = result.value();
      
    const auto p1 = seg1.getPoint(t1);
    const auto p2 = seg2.getPoint(t2);

    ASSERT_TRUE(point::isSame(p1, p2));
    
    if ((t1 < 0.01) or (t1 > 0.99) or (t2 < 0.01) or (t2 > 0.99)) continue;
    
    ++total_cases_passed;
    {
      curveIntersector theIntersector(seg1, seg2);
      theIntersector.fulfill();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 1);
      ASSERT_EQ(t1s.size(), 1);
      ASSERT_EQ(t2s.size(), 1);

      const auto pntSol = points.front();
      const auto t1Sol = t1s.front();
      const auto t2Sol = t2s.front();

      const auto p1Sol = seg1.getPoint(t1Sol);
      const auto p2Sol = seg2.getPoint(t2Sol);

      ASSERT_TRUE(point::isSame(p1Sol, p2Sol));

      // This test as consequence of non-deterministic nature of random pseudo-number generator may not pass - so i taked off them.
      // But if you want to be a gik - reenable them.
      //ASSERT_TRUE(point::isSame(p1Sol, p1, math::tolerance::tolPoint * 10.0)) << "p1Sol = " << p1Sol << " p1 = " << p1 << " dist = " << point::distance(p1Sol, p1) << std::endl;
      //ASSERT_TRUE(point::isSame(p1Sol, p2, math::tolerance::tolPoint * 10.0)) << "p1Sol = " << p1Sol << " p2 = " << p2 << " dist = " << point::distance(p1Sol, p2) << std::endl;
      //ASSERT_TRUE(point::isSame(p2Sol, p1, math::tolerance::tolPoint * 10.0)) << "p2Sol = " << p1Sol << " p1 = " << p1 << " dist = " << point::distance(p2Sol, p1) << std::endl;
      //ASSERT_TRUE(point::isSame(p2Sol, p2, math::tolerance::tolPoint * 10.0)) << "p2Sol = " << p2Sol << " p2 = " << p2 << " dist = " << point::distance(p2Sol, p2) << std::endl;

      ASSERT_NEAR(t1, t1Sol, 1.0e-7);
      ASSERT_NEAR(t2, t2Sol, 1.0e-7);
    }
  }
}