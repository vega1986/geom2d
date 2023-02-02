#include "pch.h"
#include "segmentCurve.h"
#include "curveIntersector.h"

#include <iostream>
#include <array>
#include <cmath>
#include <numbers>

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_01)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    const point po = Q;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ po, point{ 1.0, po.y} });
    curves.push_back(segmentCurve{ po, point{-1.0, po.y} });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 197; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po };
        const point pend{ po.x + theXo, po.y + theYo };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 1);

      const auto intPoint = points.front();
      const auto t1 = t1s.front();
      const auto t2 = t2s.front();

      ASSERT_TRUE(point::isSame(po, intPoint));
      ASSERT_NEAR(t1, 1.0, tabserror);
      ASSERT_NEAR(t2, 0.0, tabserror);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_02)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    const point po = Q;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ point{-0.5, po.y}, point{ 0.5, po.y} });
    curves.push_back(segmentCurve{ point{ 0.5, po.y}, point{-0.5, po.y} });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 163; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po.x - theXo * 0.5, po.y - theYo * 0.5 };
        const point pend{ po.x + theXo * 0.5, po.y + theYo * 0.5 };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 1);

      const auto intPoint = points.front();
      const auto t1 = t1s.front();
      const auto t2 = t2s.front();

      ASSERT_TRUE(point::isSame(po, intPoint));
      ASSERT_NEAR(t1, 1.0, tabserror);
      ASSERT_NEAR(t2, 0.5, tabserror);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_03)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    const point po = Q;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ point{-1.0, po.y}, po });
    curves.push_back(segmentCurve{ point{ 1.0, po.y}, po });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 149; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po.x - theXo, po.y - theYo };
        const point pend{ po };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 1);

      const auto intPoint = points.front();
      const auto t1 = t1s.front();
      const auto t2 = t2s.front();

      ASSERT_TRUE(point::isSame(po, intPoint));
      ASSERT_NEAR(t1, 1.0, tabserror);
      ASSERT_NEAR(t2, 1.0, tabserror);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_04)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5};
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    const point po = (P + Q) / 2.0;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ po, point{ 1.0, po.y} });
    curves.push_back(segmentCurve{ po, point{-1.0, po.y} });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 239; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po };
        const point pend{ po.x + theXo, po.y + theYo };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 1);

      const auto intPoint = points.front();
      const auto t1 = t1s.front();
      const auto t2 = t2s.front();

      ASSERT_TRUE(point::isSame(po, intPoint));
      ASSERT_NEAR(t1, 0.5, tabserror);
      ASSERT_NEAR(t2, 0.0, tabserror);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_05)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    const point po = (P + Q) / 2.0;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ point{-0.5, po.y}, point{ 0.5, po.y} });
    curves.push_back(segmentCurve{ point{ 0.5, po.y}, point{-0.5, po.y} });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 113; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po.x - theXo * 0.5, po.y - theYo * 0.5 };
        const point pend{ po.x + theXo * 0.5, po.y + theYo * 0.5 };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 1);

      const auto intPoint = points.front();
      const auto t1 = t1s.front();
      const auto t2 = t2s.front();

      ASSERT_TRUE(point::isSame(po, intPoint));
      ASSERT_NEAR(t1, 0.5, tabserror);
      ASSERT_NEAR(t2, 0.5, tabserror);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_06)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    const point po = (P + Q) / 2.0;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ point{-1.0, po.y }, po });
    curves.push_back(segmentCurve{ point{ 1.0, po.y }, po });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 239; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po.x - theXo, po.y - theYo };
        const point pend{ po };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 1);

      const auto intPoint = points.front();
      const auto t1 = t1s.front();
      const auto t2 = t2s.front();

      ASSERT_TRUE(point::isSame(po, intPoint));
      ASSERT_NEAR(t1, 0.5, tabserror);
      ASSERT_NEAR(t2, 1.0, tabserror);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_07)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    const point po = P;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ po, point{ 1.0, po.y} });
    curves.push_back(segmentCurve{ po, point{-1.0, po.y} });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 107; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po };
        const point pend{ po.x + theXo, po.y + theYo };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 1);

      const auto intPoint = points.front();
      const auto t1 = t1s.front();
      const auto t2 = t2s.front();

      ASSERT_TRUE(point::isSame(po, intPoint));
      ASSERT_NEAR(t1, 0.0, tabserror);
      ASSERT_NEAR(t2, 0.0, tabserror);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_08)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    const point po = P;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ point{-0.5, po.y}, point{ 0.5, po.y} });
    curves.push_back(segmentCurve{ point{ 0.5, po.y}, point{-0.5, po.y} });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 113; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po.x - theXo * 0.5, po.y - theYo * 0.5 };
        const point pend{ po.x + theXo * 0.5, po.y + theYo * 0.5 };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 1);

      const auto intPoint = points.front();
      const auto t1 = t1s.front();
      const auto t2 = t2s.front();

      ASSERT_TRUE(point::isSame(po, intPoint));
      ASSERT_NEAR(t1, 0.0, tabserror);
      ASSERT_NEAR(t2, 0.5, tabserror);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_09)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    const point po = P;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ point{-1.0, po.y}, po });
    curves.push_back(segmentCurve{ point{ 1.0, po.y}, po });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 257; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po.x - theXo, po.y - theYo };
        const point pend{ po };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 1);

      const auto intPoint = points.front();
      const auto t1 = t1s.front();
      const auto t2 = t2s.front();

      ASSERT_TRUE(point::isSame(po, intPoint));
      ASSERT_NEAR(t1, 0.0, tabserror);
      ASSERT_NEAR(t2, 1.0, tabserror);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_10)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    constexpr double alpha = -0.1;
    constexpr double beta = 1.0 - alpha;
    const point po = P * beta + Q * alpha;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ po, point{ 1.0, po.y} });
    curves.push_back(segmentCurve{ po, point{-1.0, po.y} });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 107; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po };
        const point pend{ po.x + theXo, po.y + theYo };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 0);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_11)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    constexpr double alpha = -0.1;
    constexpr double beta = 1.0 - alpha;
    const point po = P * beta + Q * alpha;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ point{-0.5, po.y}, point{ 0.5, po.y} });
    curves.push_back(segmentCurve{ point{ 0.5, po.y}, point{-0.5, po.y} });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 113; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po.x - theXo * 0.5, po.y - theYo * 0.5 };
        const point pend{ po.x + theXo * 0.5, po.y + theYo * 0.5 };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 0);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_12)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    constexpr double alpha = -0.1;
    constexpr double beta = 1.0 - alpha;
    const point po = P * beta + Q * alpha;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ point{-1.0, po.y}, po });
    curves.push_back(segmentCurve{ point{ 1.0, po.y}, po });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 257; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po.x - theXo, po.y - theYo };
        const point pend{ po };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 0);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_13)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    constexpr double alpha = 1.1;
    constexpr double beta = 1.0 - alpha;
    const point po = P * beta + Q * alpha;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ po, point{ 1.0, po.y} });
    curves.push_back(segmentCurve{ po, point{-1.0, po.y} });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 107; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po };
        const point pend{ po.x + theXo, po.y + theYo };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 0);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_14)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    constexpr double alpha = 1.1;
    constexpr double beta = 1.0 - alpha;
    const point po = P * beta + Q * alpha;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ point{-0.5, po.y}, point{ 0.5, po.y} });
    curves.push_back(segmentCurve{ point{ 0.5, po.y}, point{-0.5, po.y} });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 113; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po.x - theXo * 0.5, po.y - theYo * 0.5 };
        const point pend{ po.x + theXo * 0.5, po.y + theYo * 0.5 };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 0);
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment_ver_any_15)
{
  using namespace geom2d;
  using namespace std::numbers;
  constexpr double tabserror = 1.0e-9;
  const point P{ 0.0, -0.5 };
  const point Q{ 0.0,  0.5 };
  segmentCurve horSegment{ P, Q };
  {
    constexpr double alpha = 1.1;
    constexpr double beta = 1.0 - alpha;
    const point po = P * beta + Q * alpha;
    std::vector<segmentCurve> curves;
    curves.push_back(segmentCurve{ point{-1.0, po.y}, po });
    curves.push_back(segmentCurve{ point{ 1.0, po.y}, po });

    {
      constexpr double alpha0 = pi / sqrt3;
      constexpr double dalpha = pi / ln10 / phi / sqrt2 / log2e / 7.0;
      double alpha = alpha0;
      for (size_t j = 0; j < 257; ++j)
      {
        const double theXo = sin(alpha);
        const double theYo = cos(alpha);

        const point pbegin{ po.x - theXo, po.y - theYo };
        const point pend{ po };

        segmentCurve otherSegment{ pbegin, pend };
        curves.push_back(otherSegment);
      }
    }

    for (const auto& otherSegment : curves)
    {
      curveIntersector theIntersector{ horSegment, otherSegment };
      theIntersector.perform();

      const auto points = theIntersector.getSolutionPoints();
      const auto t1s = theIntersector.getSolutionT1();
      const auto t2s = theIntersector.getSolutionT2();

      ASSERT_EQ(points.size(), 0);
    }
  }
}

//*********************************************************************************************************************
