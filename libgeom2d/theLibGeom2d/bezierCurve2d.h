#pragma once

#include "baseCurve.h"
#include "point.h"
#include "theLibGeom2d.h"

#include <vector>

namespace geom2d
{

  // Класс криввых Безье на плоскости
  class bezierCurve final : public baseCurve
  {
  private:

    // вспомогательная функция,
    // вызываемая из конструктора
    // Функция вычисляет вспомогательные коэффициенты
    // используемые для вычисления положение точки на кривой
    // для заданного параметра t
    THELIBGEOM2D_API
    void initialize_collateral ();

  public:

    // Принимает в качестве входного параметра произвольный итерабельный контейнер (такой как std::list, std::vector, std::array, std::initializer_list)
    template<class pointsContainer>
    bezierCurve (const pointsContainer& container)
      : baseCurve()
    {
      for (const auto p : container)
      {
        points.push_back(p);
      }
      initialize_collateral();
    }

    THELIBGEOM2D_API
    virtual ~bezierCurve() = default;

    // получить точку на кривой для заданного параметра t
    // концептуальная формула: sum(j from 0 to n - 1) { beta(j) * t^j * (1 - t)^(n - 1 - j) }
    THELIBGEOM2D_API
    virtual point getPoint (double t) const;

    // получить производные по t от точки на кривой для заданного параметра t
    // концептуальная формаула: sum(j from 0 to n - 1)
    // { j * beta(j) * t^(j - 1) * (1 - t)^(n - 1 - j) - 
    //  - (n - 1 - j) * beta(j) * t^j * (1 - t)^(n - 2 - j) }
    THELIBGEOM2D_API
    virtual point getVelocity (double t) const;

    // минимальное значение параметра t
    THELIBGEOM2D_API
    virtual double parameterMin() const
    {
      return 0.0;
    }

    // максимальное значение параметра t
    THELIBGEOM2D_API
    virtual double parameterMax() const
    {
      return 1.0;
    }
    
  private:
    std::vector<point> points;

    // n         : количество точек сплайна, j = [0; n - 1]
    // x[j]      : x-координата j-ой точки сплайна
    // y[j]      : у-координата j-ой точки сплайна
    // C(n, k)   : число сочетаний из n по k
    // beta[j].x : x[j] * C(n - 1, j)
    // beta[j].y : y[j] * C(n - 1, j)
    std::vector<point> beta;

    // gamma[j]  : beta[j] * j
    std::vector<point> gamma;

    // ksi[j]    : beta[j] * (n - 1 - j)
    std::vector<point> ksi;
  };

}
