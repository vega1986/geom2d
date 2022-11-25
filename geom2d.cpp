// geom2d.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <chrono>

#include "math.h"
#include "bezierCurve2d.h"
#include "findFunctionRoots.h"

int main()
{
  //using namespace math;
  //auto fac5 = factorial(5);
  //std::cout << "factorial(5) = " << fac5 << std::endl;

  //auto combx = comb(35, 35);
  //std::cout << "combx = " << combx << std::endl;

  //auto pp = std::pow(0.0, 1.0);
  //std::cout << "pow = " << pp << std::endl;

  using namespace geom2d;
  bezierCurve myCurve(std::initializer_list<point>{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}});

  constexpr double tcheck = 0.5;
  auto pp = myCurve.getPoint(tcheck);
  auto ppDerviv = myCurve.getVelocity(tcheck);
  std::cout << "pp.x = " << pp.x << " pp.y = " << pp.y << std::endl;
  std::cout << "DER.x = " << ppDerviv.x << " DER.y = " << ppDerviv.y << std::endl;

  // тестирование функции math::power
  // std::cout << "3.0^11 = " << math::power(3.0, 11);
  bool qq = false;
  
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (size_t i = 0; i < 10000000; ++i)
  {
    const double w = math::power(3.0, i % 13);
    //const double w = std::powf(3.0, i % 13);
    qq = w > 658.0;
    //std::cout << w << std::endl;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << qq << std::endl;

  std::cout << "Time difference = " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " [ms]" << std::endl;

  // протестируем поиск корней функции на отрезке
  // корни 1 и 2
  auto func1 = [](double t) -> double
  {
    return t * t - 3.0 * t + 2;
  };
  // корни 1, 2, 5, 7
  auto func2 = [](double t) -> double
  {
    return t * t * t * t - 15 * t * t * t + 73.0 * t * t - 129.0 * t + 70;
  };
  // корни 1 2 5 7 9  -3 -5
  auto func3 = [](double t) -> double
  {
    return                  std::pow(t, 7.0)
                -    16.0 * std::pow(t, 6.0)
                +    31.0 * std::pow(t, 5.0)
                +   518.0 * std::pow(t, 4.0)
                -  1937.0 * std::pow(t, 3.0)
                -  2572.0 * std::pow(t, 2.0)
                + 13425.0 * t
                -  9450.0;
  };
  auto roots = math::findFunctionRoots(-5.0 - 0.1, 9.0 + 0.1, func3);
  std::cout << std::endl << "ROOTS:" << std::endl << std::endl;
  for (const auto r : roots)
  {
    std::cout << r << std::endl;
  }

  char ch = 0;
  std::cin >> ch;
  return 1;
}
