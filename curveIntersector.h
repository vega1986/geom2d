#pragma once

#include "baseCurve.h"
#include "point.h"
#include "IntersecctionSolutionType.h"

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

    // Исследуем решение, если первая и вторая кривая - точки
    static
      std::optional<geom2d::IntersecctionSolutionType>
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

    // Исследуем решение, если первая кривая - точка,
    // а вторая 'X = const'
    static
      std::optional<geom2d::IntersecctionSolutionType>
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
      std::optional<geom2d::IntersecctionSolutionType>
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
      std::optional<geom2d::IntersecctionSolutionType>
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
      std::optional<geom2d::IntersecctionSolutionType>
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
      std::optional<geom2d::IntersecctionSolutionType>
      execPlatoXAndPlatoX
      (
        const double tmin1,
        const double tmax1,
        const baseCurve& curvePlatoX1,
        const double tmin2,
        const double tmax2,
        const baseCurve& curvePlatoX2
      );

    // Исследуем решение, если первая кривая - PlatoX, вторая - PlatoY
    static
      std::optional<geom2d::IntersecctionSolutionType>
      execPlatoXAndPlatoY
      (
        const double tminOfPlatox,
        const double tmaxOfPlatox,
        const baseCurve& curvePlatoX,
        const double tminOfPlatoy,
        const double tmaxOfPlatoy,
        const baseCurve& curvePlatoY
      );

    // Находим точку двух кривых
    // в предположении что они точно пересекаются 1 раз
    // (если одна кривая убывает, а другая возрастает вдоль оси X или Y то такие кривые точно пересекаются 1 раз)
    // так же эта функция будет справедлива когда кривые квази параллельны осям координат
    static
      std::optional<geom2d::IntersecctionSolutionType>
      findUniqueIntersection
      (
        const double tmin1,
        const double tmax1,
        const baseCurve& curve1,
        const double tmin2,
        const double tmax2,
        const baseCurve& curve2
      );

    // Находим точку двух кривых
    // в предположении что они точно пересекаются 1 раз
    // (если одна кривая убывает, а другая возрастает вдоль оси X или Y то такие кривые точно пересекаются 1 раз)
    // так же эта функция будет справедлива когда кривые квази параллельны осям координат
    // Данная функция отличается от предыдущей тем, что известна референстная кривая
    // Метод деления отрезка пополам при поиске точки пересечения применяется ко второй кривой (не референсной)
    // Наибольший размер референсной кривой - вдоль оси X
    static
      std::optional<geom2d::IntersecctionSolutionType>
      findUniqueIntersectionRefAlongX
      (
        const double trefmin,
        const double trefmax,
        const baseCurve& referenceCurve,
        const double tothmin,
        const double tothmax,
        const baseCurve& otherCurve
      );

    // -//-
    // Наибольший размер референсной кривой - вдоль оси Y
    static
      std::optional<geom2d::IntersecctionSolutionType>
      findUniqueIntersectionRefAlongY
      (
        const double trefmin,
        const double trefmax,
        const baseCurve& referenceCurve,
        const double tothmin,
        const double tothmax,
        const baseCurve& otherCurve
      );

    // Исследуем решение, если первая кривая - PlatoX, вторая - Any
    static
      std::optional<geom2d::IntersecctionSolutionType>
      execPlatoXAndAny
      (
        const double tminOfPlatox,
        const double tmaxOfPlatox,
        const baseCurve& curvePlatoX,
        const double tminOfAny,
        const double tmaxOfAny,
        const baseCurve& curveAny
      );

    // Исследуем решение, если первая кривая - PlatoY, вторая - Any
    static
      std::optional<geom2d::IntersecctionSolutionType>
      execPlatoYAndAny
      (
        const double tminOfPlatoy,
        const double tmaxOfPlatoy,
        const baseCurve& curvePlatoY,
        const double tminOfAny,
        const double tmaxOfAny,
        const baseCurve& curveAny
      );

    // Исследуем решение, если первая кривая - Normal, вторая - Screen
    static
      std::optional<geom2d::IntersecctionSolutionType>
      execNormalAndScreen
      (
        const double tminOfNormal,
        const double tmaxOfNormal,
        const baseCurve& curveNormal,
        const double tminOfScreen,
        const double tmaxOfScreen,
        const baseCurve& curveScreen
      );

    // Исследуем решение, если первая и вторая кривые - Normal
    static
      std::vector<geom2d::IntersecctionSolutionType>
      execNormalAndNormal
      (
        const double tmin1,
        const double tmax1,
        const baseCurve& curve1,
        const double tmin2,
        const double tmax2,
        const baseCurve& curve2
      );

    // Исследуем решение, если первая и вторая кривые - Screen
    static
      std::vector<geom2d::IntersecctionSolutionType>
      execScreenAndScreen
      (
        const double tmin1,
        const double tmax1,
        const baseCurve& curve1,
        const double tmin2,
        const double tmax2,
        const baseCurve& curve2
      );

  private:

    const baseCurve & m_curve1;
    const baseCurve & m_curve2;

    std::vector<point> solutionPoints;
    std::vector<double> solutionParameterOnCurve1;
    std::vector<double> solutionParameterOnCurve2;
  };
}
