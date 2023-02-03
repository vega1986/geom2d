#include "pch.h"
#include "segmentCurve.h"
#include "arbitraryCurve2d.h"
#include "findFunctionRoots.h"
#include "curveIntersector.h"

#include <numbers>
#include <random>
#include <cmath>

TEST(INTERSECTOR, circle_and_segment)
{
  constexpr double R = 1.0;
  constexpr double deltaR = 1.0e-2;
  constexpr size_t n = 1024 * 128;

  using namespace geom2d;
  using namespace std::numbers;

  constexpr double fiMin = 0.0;
  constexpr double fiMax = pi / 2.0;
  constexpr double deltaFi = pi * 1.0e-4;

  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution disR1(deltaR, R - deltaR);
  std::uniform_real_distribution disR2(R + deltaR, R * 5.0);
  std::uniform_real_distribution disFi(fiMin + deltaFi, fiMax - deltaFi);

  auto px = [R](const double fi) -> const double { return R * std::cos(fi); };
  auto py = [R](const double fi) -> const double { return R * std::sin(fi); };

  auto ux = [R](const double fi) -> const double { return - R * std::sin(fi); };
  auto uy = [R](const double fi) -> const double { return   R * std::cos(fi); };

  arbitraryCurve arc{ px, py, ux, uy, fiMin , fiMax };

  for (size_t j = 0; j < n; ++j)
  {
    const double r1 = disR1(gen);
    const double fi1 = disFi(gen);

    const point P{ r1 * std::cos(fi1), r1 * std::sin(fi1) };

    const double r2 = disR2(gen);
    const double fi2 = disFi(gen);

    const point Q{ r2 * std::cos(fi2), r2 * std::sin(fi2) };

    const segmentCurve seg{ P, Q };

    auto func = [&seg, R](const double t) -> double
    {
      return seg.getPoint(t).length() - R;
    };

    const auto tofR = math::findUniqueFunctionRoot(seg.parameterMin(), seg.parameterMax(), func);
    const auto pointOfR = seg.getPoint(tofR);

    ASSERT_NEAR(pointOfR.length(), R, math::tolerance::tolPoint);

    curveIntersector intersector{ seg, arc };
    intersector.fulfill();

    const auto points = intersector.getSolutionPoints();
    ASSERT_EQ(points.size(), 1);

    const auto intPoint = points.front();
    ASSERT_NEAR(intPoint.length(), R, math::tolerance::tolPoint);
    ASSERT_TRUE(point::distance(pointOfR, intPoint) <= (math::tolerance::tolPoint * 1.0e+2)) << std::endl << " $ dist = " << point::distance(pointOfR, intPoint) << std::endl;
  }
}
