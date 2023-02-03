#include "pch.h"
#include "findFunctionRoots.h"

#include <iomanip>

TEST(OneArgumentFunctionRoots, Polynomial)
{
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
      - 16.0 * std::pow(t, 6.0)
      + 31.0 * std::pow(t, 5.0)
      + 518.0 * std::pow(t, 4.0)
      - 1937.0 * std::pow(t, 3.0)
      - 2572.0 * std::pow(t, 2.0)
      + 13425.0 * t
      - 9450.0;
  };
  auto roots = math::findFunctionRoots(-5.0 - 0.1, 9.0 + 0.1, func3);

  ASSERT_EQ(roots.size(), 7);

  auto it = roots.cbegin();
  ASSERT_NEAR(*it, -5.0, 1.0e-6); ++it;
  ASSERT_NEAR(*it, -3.0, 1.0e-6); ++it;
  ASSERT_NEAR(*it,  1.0, 1.0e-6); ++it;
  ASSERT_NEAR(*it,  2.0, 1.0e-6); ++it;
  ASSERT_NEAR(*it,  5.0, 1.0e-6); ++it;
  ASSERT_NEAR(*it,  7.0, 1.0e-6); ++it;
  ASSERT_NEAR(*it,  9.0, 1.0e-6); ++it;
}