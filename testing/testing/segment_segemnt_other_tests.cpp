#include "pch.h"
#include "segmentCurve.h"
#include "curveIntersector.h"

#include <iostream>
#include <array>
#include <cmath>
#include <numbers>
#include <vector>
#include <cmath>

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment__main_other_01)
{
  using namespace geom2d;
  using namespace std::numbers;

  constexpr double angleTolerance = pi * 1.0e-4;
  constexpr double tabserror = 1.0e-7;
  
  std::vector<double> mainSegmentAngles;
  constexpr double angleInitMain = pi / ln10 / phi / sqrt2 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    mainSegmentAngles.push_back(angleInitMain + 2.0 * pi * j / 100.0);


  std::vector<double> otherSegmentAngles;
  constexpr double angleInitOther = pi / ln10 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    otherSegmentAngles.push_back(angleInitOther + 2.0 * pi * j / 100.0);

  constexpr point omain{ 0.0, 0.0 };
  for (const auto mainAng : mainSegmentAngles)
  {
    const point beginMain{ omain.x - 0.5 * cos(mainAng), omain.y - 0.5 * sin(mainAng) };
    const point endMain{ omain.x + 0.5 * cos(mainAng), omain.y + 0.5 * sin(mainAng) };
    for (const auto otherAng : otherSegmentAngles)
    {
      if (std::fabs(mainAng - otherAng) <= angleTolerance) continue;
      
      const point beginOther{ endMain };
      const point endOther{ beginOther.x + cos(otherAng), beginOther.y + sin(otherAng) };
      
      segmentCurve mainCurve{ beginMain, endMain };
      segmentCurve otherCurve{ beginOther, endOther };
      
      curveIntersector intersector{ mainCurve, otherCurve };
      intersector.fulfill();

      const auto intPoints = intersector.getSolutionPoints();
      const auto T1s = intersector.getSolutionT1();
      const auto T2s = intersector.getSolutionT2();
      ASSERT_EQ(intPoints.size(), 1);

      ASSERT_NEAR(T1s[0], 1.0, tabserror);
      ASSERT_NEAR(T2s[0], 0.0, tabserror);

      const auto pon1 = mainCurve.getPoint(T1s[0]);
      const auto pon2 = otherCurve.getPoint(T2s[0]);

      ASSERT_TRUE(point::isSame(pon1, pon2));
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment__main_other_02)
{
  using namespace geom2d;
  using namespace std::numbers;

  constexpr double angleTolerance = pi * 1.0e-4;
  constexpr double tabserror = 1.0e-7;

  std::vector<double> mainSegmentAngles;
  constexpr double angleInitMain = pi / ln10 / phi / sqrt2 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    mainSegmentAngles.push_back(angleInitMain + 2.0 * pi * j / 100.0);


  std::vector<double> otherSegmentAngles;
  constexpr double angleInitOther = pi / ln10 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    otherSegmentAngles.push_back(angleInitOther + 2.0 * pi * j / 100.0);

  constexpr point omain{ 0.0, 0.0 };
  for (const auto mainAng : mainSegmentAngles)
  {
    const point beginMain{ omain.x - 0.5 * cos(mainAng), omain.y - 0.5 * sin(mainAng) };
    const point endMain{ omain.x + 0.5 * cos(mainAng), omain.y + 0.5 * sin(mainAng) };
    for (const auto otherAng : otherSegmentAngles)
    {
      if (std::fabs(mainAng - otherAng) <= angleTolerance) continue;

      const point beginOther{ endMain.x - 0.5 * cos(otherAng), endMain.y - 0.5 * sin(otherAng) };
      const point endOther{ endMain.x + 0.5 * cos(otherAng), endMain.y + 0.5 * sin(otherAng) };

      segmentCurve mainCurve{ beginMain, endMain };
      segmentCurve otherCurve{ beginOther, endOther };

      curveIntersector intersector{ mainCurve, otherCurve };
      intersector.fulfill();

      const auto intPoints = intersector.getSolutionPoints();
      const auto T1s = intersector.getSolutionT1();
      const auto T2s = intersector.getSolutionT2();
      ASSERT_EQ(intPoints.size(), 1);

      ASSERT_NEAR(T1s[0], 1.0, tabserror);
      ASSERT_NEAR(T2s[0], 0.5, tabserror);

      const auto pon1 = mainCurve.getPoint(T1s[0]);
      const auto pon2 = otherCurve.getPoint(T2s[0]);

      ASSERT_TRUE(point::isSame(pon1, pon2));
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment__main_other_03)
{
  using namespace geom2d;
  using namespace std::numbers;

  constexpr double angleTolerance = pi * 1.0e-4;
  constexpr double tabserror = 1.0e-7;

  std::vector<double> mainSegmentAngles;
  constexpr double angleInitMain = pi / ln10 / phi / sqrt2 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    mainSegmentAngles.push_back(angleInitMain + 2.0 * pi * j / 100.0);


  std::vector<double> otherSegmentAngles;
  constexpr double angleInitOther = pi / ln10 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    otherSegmentAngles.push_back(angleInitOther + 2.0 * pi * j / 100.0);

  constexpr point omain{ 0.0, 0.0 };
  for (const auto mainAng : mainSegmentAngles)
  {
    const point beginMain{ omain.x - 0.5 * cos(mainAng), omain.y - 0.5 * sin(mainAng) };
    const point endMain{ omain.x + 0.5 * cos(mainAng), omain.y + 0.5 * sin(mainAng) };
    for (const auto otherAng : otherSegmentAngles)
    {
      if (std::fabs(mainAng - otherAng) <= angleTolerance) continue;

      const point beginOther{ endMain.x - cos(otherAng), endMain.y - sin(otherAng) };
      const point endOther{ endMain };

      segmentCurve mainCurve{ beginMain, endMain };
      segmentCurve otherCurve{ beginOther, endOther };

      curveIntersector intersector{ mainCurve, otherCurve };
      intersector.fulfill();

      const auto intPoints = intersector.getSolutionPoints();
      const auto T1s = intersector.getSolutionT1();
      const auto T2s = intersector.getSolutionT2();
      ASSERT_EQ(intPoints.size(), 1);

      ASSERT_NEAR(T1s[0], 1.0, tabserror);
      ASSERT_NEAR(T2s[0], 1.0, tabserror);

      const auto pon1 = mainCurve.getPoint(T1s[0]);
      const auto pon2 = otherCurve.getPoint(T2s[0]);

      ASSERT_TRUE(point::isSame(pon1, pon2));
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment__main_other_04)
{
  using namespace geom2d;
  using namespace std::numbers;

  constexpr double angleTolerance = pi * 1.0e-4;
  constexpr double tabserror = 1.0e-7;

  std::vector<double> mainSegmentAngles;
  constexpr double angleInitMain = pi / ln10 / phi / sqrt2 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    mainSegmentAngles.push_back(angleInitMain + 2.0 * pi * j / 100.0);


  std::vector<double> otherSegmentAngles;
  constexpr double angleInitOther = pi / ln10 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    otherSegmentAngles.push_back(angleInitOther + 2.0 * pi * j / 100.0);

  constexpr point omain{ 0.0, 0.0 };
  for (const auto mainAng : mainSegmentAngles)
  {
    const point beginMain{ omain.x - 0.5 * cos(mainAng), omain.y - 0.5 * sin(mainAng) };
    const point endMain{ omain.x + 0.5 * cos(mainAng), omain.y + 0.5 * sin(mainAng) };
    for (const auto otherAng : otherSegmentAngles)
    {
      if (std::fabs(mainAng - otherAng) <= angleTolerance) continue;

      const point poOnMain{ (beginMain + endMain) / 2.0 };
      const point beginOther{ poOnMain };
      const point endOther{ poOnMain.x + cos(otherAng), poOnMain.y + sin(otherAng) };

      segmentCurve mainCurve{ beginMain, endMain };
      segmentCurve otherCurve{ beginOther, endOther };

      curveIntersector intersector{ mainCurve, otherCurve };
      intersector.fulfill();

      const auto intPoints = intersector.getSolutionPoints();
      const auto T1s = intersector.getSolutionT1();
      const auto T2s = intersector.getSolutionT2();
      ASSERT_EQ(intPoints.size(), 1);

      ASSERT_NEAR(T1s[0], 0.5, tabserror);
      ASSERT_NEAR(T2s[0], 0.0, tabserror);

      const auto pon1 = mainCurve.getPoint(T1s[0]);
      const auto pon2 = otherCurve.getPoint(T2s[0]);

      ASSERT_TRUE(point::isSame(pon1, pon2));
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment__main_other_05)
{
  using namespace geom2d;
  using namespace std::numbers;

  constexpr double angleTolerance = pi * 1.0e-4;
  constexpr double tabserror = 1.0e-7;

  std::vector<double> mainSegmentAngles;
  constexpr double angleInitMain = pi / ln10 / phi / sqrt2 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    mainSegmentAngles.push_back(angleInitMain + 2.0 * pi * j / 100.0);


  std::vector<double> otherSegmentAngles;
  constexpr double angleInitOther = pi / ln10 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    otherSegmentAngles.push_back(angleInitOther + 2.0 * pi * j / 100.0);

  constexpr point omain{ 0.0, 0.0 };
  size_t index = 0;
  for (const auto mainAng : mainSegmentAngles)
  {
    const point beginMain{ omain.x - 0.5 * cos(mainAng), omain.y - 0.5 * sin(mainAng) };
    const point endMain{ omain.x + 0.5 * cos(mainAng), omain.y + 0.5 * sin(mainAng) };
    for (const auto otherAng : otherSegmentAngles)
    {
      if (std::fabs(mainAng - otherAng) <= angleTolerance) continue;

      const point poOnMain{ (beginMain + endMain) / 2.0 };
      const point beginOther{ poOnMain.x - 0.5 * cos(otherAng), poOnMain.y - 0.5 * sin(otherAng) };
      const point endOther{ poOnMain.x + 0.5 * cos(otherAng), poOnMain.y + 0.5 * sin(otherAng) };

      segmentCurve mainCurve{ beginMain, endMain };
      segmentCurve otherCurve{ beginOther, endOther };

      curveIntersector intersector{ mainCurve, otherCurve };
      ++index;
      intersector.fulfill();

      const auto intPoints = intersector.getSolutionPoints();
      const auto T1s = intersector.getSolutionT1();
      const auto T2s = intersector.getSolutionT2();
      
      ASSERT_EQ(intPoints.size(), 1);
      ASSERT_NEAR(T1s[0], 0.5, tabserror);
      ASSERT_NEAR(T2s[0], 0.5, tabserror);

      const auto pon1 = mainCurve.getPoint(T1s[0]);
      const auto pon2 = otherCurve.getPoint(T2s[0]);

      ASSERT_TRUE(point::isSame(pon1, pon2));
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment__main_other_06)
{
  using namespace geom2d;
  using namespace std::numbers;

  constexpr double angleTolerance = pi * 1.0e-4;
  constexpr double tabserror = 1.0e-7;

  std::vector<double> mainSegmentAngles;
  constexpr double angleInitMain = pi / ln10 / phi / sqrt2 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    mainSegmentAngles.push_back(angleInitMain + 2.0 * pi * j / 100.0);


  std::vector<double> otherSegmentAngles;
  constexpr double angleInitOther = pi / ln10 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    otherSegmentAngles.push_back(angleInitOther + 2.0 * pi * j / 100.0);

  constexpr point omain{ 0.0, 0.0 };
  for (const auto mainAng : mainSegmentAngles)
  {
    const point beginMain{ omain.x - 0.5 * cos(mainAng), omain.y - 0.5 * sin(mainAng) };
    const point endMain{ omain.x + 0.5 * cos(mainAng), omain.y + 0.5 * sin(mainAng) };
    for (const auto otherAng : otherSegmentAngles)
    {
      if (std::fabs(mainAng - otherAng) <= angleTolerance) continue;

      const point poOnMain{ (beginMain + endMain) / 2.0 };
      const point beginOther{ poOnMain.x - cos(otherAng), poOnMain.y - sin(otherAng) };
      const point endOther{ poOnMain };

      segmentCurve mainCurve{ beginMain, endMain };
      segmentCurve otherCurve{ beginOther, endOther };

      curveIntersector intersector{ mainCurve, otherCurve };
      intersector.fulfill();

      const auto intPoints = intersector.getSolutionPoints();
      const auto T1s = intersector.getSolutionT1();
      const auto T2s = intersector.getSolutionT2();
      ASSERT_EQ(intPoints.size(), 1);

      ASSERT_NEAR(T1s[0], 0.5, tabserror);
      ASSERT_NEAR(T2s[0], 1.0, tabserror);

      const auto pon1 = mainCurve.getPoint(T1s[0]);
      const auto pon2 = otherCurve.getPoint(T2s[0]);

      ASSERT_TRUE(point::isSame(pon1, pon2));
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment__main_other_07)
{
  using namespace geom2d;
  using namespace std::numbers;

  constexpr double angleTolerance = pi * 1.0e-4;
  constexpr double tabserror = 1.0e-7;

  std::vector<double> mainSegmentAngles;
  constexpr double angleInitMain = pi / ln10 / phi / sqrt2 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    mainSegmentAngles.push_back(angleInitMain + 2.0 * pi * j / 100.0);


  std::vector<double> otherSegmentAngles;
  constexpr double angleInitOther = pi / ln10 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    otherSegmentAngles.push_back(angleInitOther + 2.0 * pi * j / 100.0);

  constexpr point omain{ 0.0, 0.0 };
  for (const auto mainAng : mainSegmentAngles)
  {
    const point beginMain{ omain.x - 0.5 * cos(mainAng), omain.y - 0.5 * sin(mainAng) };
    const point endMain{ omain.x + 0.5 * cos(mainAng), omain.y + 0.5 * sin(mainAng) };
    for (const auto otherAng : otherSegmentAngles)
    {
      if (std::fabs(mainAng - otherAng) <= angleTolerance) continue;

      const point poOnMain{ beginMain };
      const point beginOther{ poOnMain };
      const point endOther{ poOnMain.x + cos(otherAng), poOnMain.y + sin(otherAng) };

      segmentCurve mainCurve{ beginMain, endMain };
      segmentCurve otherCurve{ beginOther, endOther };

      curveIntersector intersector{ mainCurve, otherCurve };
      intersector.fulfill();

      const auto intPoints = intersector.getSolutionPoints();
      const auto T1s = intersector.getSolutionT1();
      const auto T2s = intersector.getSolutionT2();
      ASSERT_EQ(intPoints.size(), 1);

      ASSERT_NEAR(T1s[0], 0.0, tabserror);
      ASSERT_NEAR(T2s[0], 0.0, tabserror);

      const auto pon1 = mainCurve.getPoint(T1s[0]);
      const auto pon2 = otherCurve.getPoint(T2s[0]);

      ASSERT_TRUE(point::isSame(pon1, pon2));
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment__main_other_08)
{
  using namespace geom2d;
  using namespace std::numbers;

  constexpr double angleTolerance = pi * 1.0e-4;
  constexpr double tabserror = 1.0e-7;

  std::vector<double> mainSegmentAngles;
  constexpr double angleInitMain = pi / ln10 / phi / sqrt2 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    mainSegmentAngles.push_back(angleInitMain + 2.0 * pi * j / 100.0);


  std::vector<double> otherSegmentAngles;
  constexpr double angleInitOther = pi / ln10 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    otherSegmentAngles.push_back(angleInitOther + 2.0 * pi * j / 100.0);

  constexpr point omain{ 0.0, 0.0 };
  size_t index = 0;
  for (const auto mainAng : mainSegmentAngles)
  {
    const point beginMain{ omain.x - 0.5 * cos(mainAng), omain.y - 0.5 * sin(mainAng) };
    const point endMain{ omain.x + 0.5 * cos(mainAng), omain.y + 0.5 * sin(mainAng) };
    for (const auto otherAng : otherSegmentAngles)
    {
      if (std::fabs(mainAng - otherAng) <= angleTolerance) continue;

      const point poOnMain{ beginMain };
      const point beginOther{ poOnMain.x - 0.5 * cos(otherAng), poOnMain.y - 0.5 * sin(otherAng) };
      const point endOther{ poOnMain.x + 0.5 * cos(otherAng), poOnMain.y + 0.5 * sin(otherAng) };

      segmentCurve mainCurve{ beginMain, endMain };
      segmentCurve otherCurve{ beginOther, endOther };

      curveIntersector intersector{ mainCurve, otherCurve };
      ++index;
      intersector.fulfill();

      const auto intPoints = intersector.getSolutionPoints();
      const auto T1s = intersector.getSolutionT1();
      const auto T2s = intersector.getSolutionT2();

      ASSERT_EQ(intPoints.size(), 1);
      ASSERT_NEAR(T1s[0], 0.0, tabserror);
      ASSERT_NEAR(T2s[0], 0.5, tabserror);

      const auto pon1 = mainCurve.getPoint(T1s[0]);
      const auto pon2 = otherCurve.getPoint(T2s[0]);

      ASSERT_TRUE(point::isSame(pon1, pon2));
    }
  }
}

//*********************************************************************************************************************

TEST(Intersector, Segment_Segment__main_other_09)
{
  using namespace geom2d;
  using namespace std::numbers;

  constexpr double angleTolerance = pi * 1.0e-4;
  constexpr double tabserror = 1.0e-7;

  std::vector<double> mainSegmentAngles;
  constexpr double angleInitMain = pi / ln10 / phi / sqrt2 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    mainSegmentAngles.push_back(angleInitMain + 2.0 * pi * j / 100.0);


  std::vector<double> otherSegmentAngles;
  constexpr double angleInitOther = pi / ln10 / log2e / 7.0;
  for (size_t j = 0; j < 100; ++j)
    otherSegmentAngles.push_back(angleInitOther + 2.0 * pi * j / 100.0);

  constexpr point omain{ 0.0, 0.0 };
  for (const auto mainAng : mainSegmentAngles)
  {
    const point beginMain{ omain.x - 0.5 * cos(mainAng), omain.y - 0.5 * sin(mainAng) };
    const point endMain{ omain.x + 0.5 * cos(mainAng), omain.y + 0.5 * sin(mainAng) };
    for (const auto otherAng : otherSegmentAngles)
    {
      if (std::fabs(mainAng - otherAng) <= angleTolerance) continue;

      const point poOnMain{ beginMain };
      const point beginOther{ poOnMain.x - cos(otherAng), poOnMain.y - sin(otherAng) };
      const point endOther{ poOnMain };

      segmentCurve mainCurve{ beginMain, endMain };
      segmentCurve otherCurve{ beginOther, endOther };

      curveIntersector intersector{ mainCurve, otherCurve };
      intersector.fulfill();

      const auto intPoints = intersector.getSolutionPoints();
      const auto T1s = intersector.getSolutionT1();
      const auto T2s = intersector.getSolutionT2();
      ASSERT_EQ(intPoints.size(), 1);

      ASSERT_NEAR(T1s[0], 0.0, tabserror);
      ASSERT_NEAR(T2s[0], 1.0, tabserror);

      const auto pon1 = mainCurve.getPoint(T1s[0]);
      const auto pon2 = otherCurve.getPoint(T2s[0]);

      ASSERT_TRUE(point::isSame(pon1, pon2));
    }
  }
}

//*********************************************************************************************************************
