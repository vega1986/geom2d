// geom2d.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <chrono>

//#include "amath.h"
#include "bezierCurve2d.h"
#include "curveIntersector.h"
#include "curveMutualDistanceCalculator.h"
//#include "findFunctionRoots.h"
//#include "curveIntersector.h"
//#include "pointCurve.h"
//#include "segmentCurve.h"

int main()
{
  using namespace geom2d;

  while (true)
  {
    std::vector<point> points1;
    std::vector<point> points2;

    std::cout << "input Bezier curve 1: ";
    
    char ch = 0;
    std::cin >> ch;
    if (ch == 'q') break;

    std::cin.unget();
    while (std::cin)
    {
      point p;
      std::cin >> p;
      points1.push_back(p);

      char symPost = 0;
      std::cin >> symPost;
      if (symPost == '.')
      {
        break;
      }
      else
      {
        std::cin.unget();
      }
    }

    std::cout << "input Bezier curve 2: ";
    while (std::cin)
    {
      point p;
      std::cin >> p;
      points2.push_back(p);

      char symPost = 0;
      std::cin >> symPost;
      if (symPost == '.')
      {
        break;
      }
      else
      {
        std::cin.unget();
      }
    }

    bezierCurve bc1{ points1 };
    bezierCurve bc2{ points2 };

    curveIntersector inter{ bc1, bc2 };
    inter.fulfill();

    const auto ips = inter.getSolutionPoints();
    const auto its1 = inter.getSolutionT1();
    const auto its2 = inter.getSolutionT2();

    if (not ips.empty())
    {
      const auto n = ips.size();
      for (size_t j = 0; j < n; ++j)
      {
        std::cout << std::endl << "  intersection point: (" << j + 1 << "): " << ips[j] << " | t1 = " << its1[j] << " | t2 = " << its2[j] << std::endl;
      }
      std::cout << std::endl;
    }
    else
    {
      curveMutualDistanceCalculator dc{ bc1, bc2 };
      dc.fulfill();
      const auto [dist, t1, t2] = dc.getExtrema();
      std::cout << std::endl << "  mutual nearest point: " << dist << " | t1 = " << t1 << " | t2 = " << t2 << std::endl;
      std::cout << std::endl;
    }
  }

  return 0;
}
