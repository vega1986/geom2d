#pragma once

#include "baseCurve.h"
#include "point.h"
#include "IntersecctionSolutionType.h"
#include "theLibGeom2d.h"
#include "curveAnalizerBase.h"

#include <set>
#include <vector>
#include <optional>
#include <ostream>


namespace geom2d
{
  // класс предназначен для поиска точек пересечения двух кривых типа baseCurve,
  class curveIntersector final : private curveAnalizerBase
  {
  private:
    enum class parseAxis : unsigned short int
    {
      X = 1, // along x
      Y = 2  // along y
    };
  public:
    THELIBGEOM2D_API curveIntersector (const baseCurve & curve1, const baseCurve & curve2)
      : curveAnalizerBase(curve1, curve2){}

    // выполнить анализ на пересечение - заполнить массив точек пересечения и значений параметров
    THELIBGEOM2D_API void fulfill();

  private:
    
    // постпроцессинг - вызывает excludeDuplicatesFromSolution
    virtual void postProcessing();

    // исключение из массива решения повторяющихся точек
    void excludeDuplicatesFromSolution();

    // см. базовый класс
    // семейство методов ниже определяет точки пересечения сегментов кривых на участках монотонности
    // точки пересечения и соответствующие параметры добавляются в массив решения
    // вызов этих методов зафиксирован в базовом классе, здесь они просто перегружены
    virtual void performScreen1Screen2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performScreen1Normal2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performScreen1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performScreen1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performScreen1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performNormal1Screen2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performNormal1Normal2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performNormal1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performNormal1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performNormal1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPlatoX1Any2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPlatoX1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPlatoX1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPlatoX1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPlatoY1Any2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPlatoY1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPlatoY1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPlatoY1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPoint1Any2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPoint1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPoint1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    virtual void performPoint1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

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
    // а вторая - либо Normal либо Screen
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

    // Находим точку пересечения двух кривых на участках монотонности
    // если одна кривая убывает, а другая возрастает вдоль оси X или Y.
    // Так же эта функция будет справедлива когда кривые квази параллельны разным осям координат.
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

    // Находим единственную точку пересечения двух кривых на участках монотонности, выделяя пересечение их ОДЗ вдоль оси X и исследуя взаимное положение концов кривых на концах ОДЗ.
    // Для нахождения точки пересечения, ОДЗ по t otherCurve делится поплам на каждом шаге цикла до уменьшения разности y-координаты точек на кривых до толеранса.
    // Метод деления отрезка пополам при поиске точки пересечения применяется ко второй кривой (не референсной).
    // Наибольший размер AABB референсной кривой - вдоль оси X.
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

    // Находим единственную точку пересечения двух кривых на участках монотонности, выделяя пересечение их ОДЗ вдоль оси Y и исследуя взаимное положение концов кривых на концах ОДЗ.
    // Для нахождения точки пересечения, ОДЗ по t otherCurve делится поплам на каждом шаге цикла до уменьшения разности x-координаты точек на кривых до толеранса.
    // Метод деления отрезка пополам при поиске точки пересечения применяется ко второй кривой (не референсной).
    // Наибольший размер AABB референсной кривой - вдоль оси Y.
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

    // Исследуем решение, если первая кривая - PlatoX, вторая - Normal или Screen.
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

    // Исследуем решение, если первая кривая - PlatoY, вторая - Normal или Screen.
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

  public:

    // печать решения в выходной поток
    THELIBGEOM2D_API void dumpIntersections(std::ostream& ost) const;

    // возвращаем массив точек пересечения
    THELIBGEOM2D_API std::vector<point> getSolutionPoints() const
    {
      return solutionPoints;
    }

    // возвращаем массив параметров на кривой 1 точек пересечения
    THELIBGEOM2D_API std::vector<double> getSolutionT1() const
    {
      return solutionParameterOnCurve1;
    }

    // возвращаем массив параметров на кривой 2 точек пересечения
    THELIBGEOM2D_API std::vector<double> getSolutionT2() const
    {
      return solutionParameterOnCurve2;
    }

  public:
    
    // Исследуем решение для двух кривых на участках монотонности если заведомо известно, что пересечение единственно.
    // Функция применима для кривых:
    //  - Screen & Normal и наоборот.
    //  - PlatoX & (Screen | Normal) и наоборот.
    //  - PlatoY & (Screen | Normal) и наоборот.
    //  - PlatoX & PlatoY и наоборот.
    // Функция вызывает findUniqueIntersection.
    // Данная функция необходима для алгоритма поиска ближайших точек двух кривых:
    // С её помощью ищется точка пересечения нормали к одной кривой со второй кривой.
    THELIBGEOM2D_API
    static
      std::optional<geom2d::IntersecctionSolutionType>
      exec_Any_and_Any_Unique_Intersection
      (
        const double tmin1,
        const double tmax1,
        const baseCurve& curve1,
        const double tmin2,
        const double tmax2,
        const baseCurve& curve2
      );


  private:

    std::vector<point> solutionPoints;
    std::vector<double> solutionParameterOnCurve1;
    std::vector<double> solutionParameterOnCurve2;

  };
}
