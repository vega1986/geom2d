#pragma once

#include "baseCurve.h"
#include "point.h"

#include <set>
#include <vector>
#include <optional>


namespace geom2d
{
  // классифицируем поведение участка кривой на участке монотонности
  enum class curveClass : unsigned short int
  {
    Screen = 1, // (dx / dt > 0 && dy / dt < 0) || (dx / dt < 0 && dy / dt > 0)
    Normal = 2, // (dx / dt > 0 && dy / dt > 0) || (dx / dt < 0 && dy / dt < 0)
    PlatoX = 3, // (dx / dt == 0 && (dy / dt>0 || dy / dt < 0)
    PlatoY = 4, // (dy / dt == 0 && (dx / dt>0 || dx / dt < 0)
    Point  = 5  // (dx / dt == 0 && dy / dt == 0)
  };

  // класс предназначен для поиска точек пересечения двух кривых типа baseCurve,
  class curveIntersector
  {
  private:
    enum class parseAxis : unsigned short int
    {
      X = 1, // along x
      Y = 2  // along y
    };
  public:
    curveIntersector (const baseCurve & curve1, const baseCurve & curve2)
      : m_curve1(curve1), m_curve2(curve2)
    {}

    // ищем точки пересечения
    void perform();
  private:

    // ищем точки пересечения на заданных участках монотонности обеих кривых
    void perform(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    // ищем значения параметра при которых направление вдижения точки
    // меняется на противоположное вдоль одной из координатных осей
    static std::set<double> rootsOfCurveVelocity(const baseCurve& curve);

    ///////////////////
    //               //
    // Точка и точка //
    //               //
    ///////////////////

    // анализируем пересечение кривых, когда они обе сводятся к точке
    static bool performPointXPoint(const point p, const point q);

    // Исследуем решение, если первая и вторая кривая - точки
    static
      std::optional<std::tuple<geom2d::point, double, double>>
        execPointAndPoint
        (
          const double tminOfPoint1,
          const double tmaxOfPoint1,
          const baseCurve& curvePoint1,
          const double tminOfPoint2,
          const double tmaxOfPoint2,
          const baseCurve& curvePoint2
        );

    //////////////////////
    //                  //
    // Точка и НЕ точка //
    //                  //
    //////////////////////

    // анализируем пересечение кривых, когда первая - точка, а вторая - плато по x (x ~ const)
    // но функция может применяться даже если x - не константа, главное, чтоб по y кривая была
    // не константа
    static std::optional<double>
      performPointVSAnyAlongY(
        const point pnt,
        const double tmin,
        const double tmax,
        const baseCurve& curve);

    // анализируем пересечение кривых, когда первая - точка, а вторая - плато по y (y ~ const)
    // но функция может применяться даже если y - не константа, главное, чтоб по x кривая была
    // не константа
    static std::optional<double>
      performPointVSAnyAlongX(
        const point pnt,
        const double tmin,
        const double tmax,
        const baseCurve& curve);

    // Исследуем решение, если первая кривая - точка,
    // а вторая 'X = const'
    static
      std::optional<std::tuple<geom2d::point, double, double>>
        execPointAndPlatoX
        (
          const double tminOfPoint,
          const double tmaxOfPoint,
          const baseCurve& curvePoint,
          const double tminOfPlatoX,
          const double tmaxOfPlatoX,
          const baseCurve& curvePlatoX
        );

    // Исследуем решение, если первая кривая - точка,
    // а вторая 'Y = const'
    static
      std::optional<std::tuple<geom2d::point, double, double>>
        execPointAndPlatoY
        (
          const double tminOfPoint,
          const double tmaxOfPoint,
          const baseCurve& curvePoint,
          const double tminOfPlatoY,
          const double tmaxOfPlatoY,
          const baseCurve& curvePlatoY
        );

    // Исследуем решение, если первая кривая - точка,
    // а вторая не плато, т.е. изменчива как по X, так и по Y.
    static
      std::optional<std::tuple<geom2d::point, double, double>>
        execPointAndAny
        (
          const double tminOfPoint,
          const double tmaxOfPoint,
          const baseCurve& curvePoint,
          const double tminOfAny,
          const double tmaxOfAny,
          const baseCurve& curveAny
        );

    // Исследуем решение, если первая и вторая кривая - PlatoY
    static
      std::optional<std::tuple<geom2d::point, double, double>>
        execPlatoYAndPlatoY
        (
          const double tmin1,
          const double tmax1,
          const baseCurve& curvePlatoY1,
          const double tmin2,
          const double tmax2,
          const baseCurve& curvePlatoY2
        );

    // Исследуем решение, если первая и вторая кривая - PlatoX
    static
      std::optional<std::tuple<geom2d::point, double, double>>
      execPlatoXAndPlatoX
      (
        const double tmin1,
        const double tmax1,
        const baseCurve& curvePlatoX1,
        const double tmin2,
        const double tmax2,
        const baseCurve& curvePlatoX2
      );

  private:

    const baseCurve & m_curve1;
    const baseCurve & m_curve2;

    std::vector<point> solutionPoints;
    std::vector<double> solutionParameterOnCurve1;
    std::vector<double> solutionParameterOnCurve2;
  };
}
