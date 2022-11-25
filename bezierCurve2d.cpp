#include <cmath>
#include "bezierCurve2d.h"
#include "math.h"

//-----------------------------------------------------------------------------

void geom2d::bezierCurve::initialize_collateral()
{
  // здесь мы инициализируем множители для формулы расчёта координат
  // кривой Безье и их производных (для дальнейшего расчёта их по параметру t)
  const auto n = points.size();
  beta.reserve(n);
  gamma.reserve(n);
  ksi.reserve(n);
  for (size_t j = 0; j < n; ++j)
  {
    const auto currentPoint = points[j];
    
    // число сочетаний из n-1 по j
    const auto combj = math::comb(n - 1, j);

    const double betax = currentPoint.x * combj;
    const double betay = currentPoint.y * combj;
    beta.push_back({betax, betay});

    const double gammax = j * betax;
    const double gammay = j * betay;
    gamma.push_back({gammax, gammay});

    const double ksix = (n - 1 - j) * betax;
    const double ksiy = (n - 1 - j) * betay;
    ksi.push_back({ksix, ksiy});
  }
}

//-----------------------------------------------------------------------------

geom2d::point geom2d::bezierCurve::getPoint(double t) const
{
  point res{ 0.0, 0.0 };

  const auto n = points.size();
  for (size_t j = 0; j < n; ++j)
  {
    const auto coeff = std::pow(t, j) * std::pow(1.0 - t, n - 1 - j);
    res.x += beta[j].x * coeff;
    res.y += beta[j].y * coeff;
  }
  return res;
}

//-----------------------------------------------------------------------------

geom2d::point geom2d::bezierCurve::getVelocity(double t) const
{
  point res{ 0.0, 0.0 };

  const auto n = points.size();
  for (size_t j = 1; j < n; ++j)
  {
    const auto coeff = std::pow(t, j - 1) * std::pow(1.0 - t, n - 1 - j);
    res.x += gamma[j].x * coeff;
    res.y += gamma[j].y * coeff;
  }
  for (size_t j = 0; j < n - 1; ++j)
  {
    const auto coeff = std::pow(t, j) * std::pow(1.0 - t, n - 2 - j);
    res.x -= ksi[j].x * coeff;
    res.y -= ksi[j].y * coeff;
  }
  return res;
}

//-----------------------------------------------------------------------------
