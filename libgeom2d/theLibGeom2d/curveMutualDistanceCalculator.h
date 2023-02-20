#pragma once
#include "theLibGeom2d.h"
#include "curveAnalizerBase.h"
#include "amath.h"
#include "segmentCurve.h"

#include <tuple>
#include <optional>

namespace geom2d
{
  // Класс для расчёта ближайших точек между двумя кривыми - использует концепцию базового класса curveAnalizerBase для поиска решения
  // (путём деления ОДЗ каждой кривой на участки монотонности и исследвоания решения для двух кривых для каждой пары ОДЗ монотонности
  class curveMutualDistanceCalculator final : private curveAnalizerBase
  {
  public:
    THELIBGEOM2D_API
      curveMutualDistanceCalculator(const baseCurve& curve1, const baseCurve& curve2)
      :curveAnalizerBase(curve1, curve2), m_dist{ math::infinite::distance }, m_t1{ 0.0 }, m_t2{0.0} {}

    // поиск ближайших точек двух кривых
    THELIBGEOM2D_API void fulfill();

    // возвращает найденное решение
    THELIBGEOM2D_API
      const std::tuple<double, double, double> getExtrema() const;

  private:
    // Функции ниже ищут решение для исходной пары кривых для различных сочетаний классов кривых на участках монотонности.
    // Вызов этих методов зафиксирован в базовом классе. Здесь они просто переопределяются
    virtual void performScreen1Screen2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//
    virtual void performScreen1Normal2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//
    virtual void performScreen1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performScreen1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performScreen1Point2 (const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performNormal1Screen2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//
    virtual void performNormal1Normal2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//
    virtual void performNormal1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performNormal1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performNormal1Point2 (const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPlatoX1Any2   (const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPlatoX1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPlatoX1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPlatoX1Point2 (const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPlatoY1Any2   (const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPlatoY1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPlatoY1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPlatoY1Point2 (const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPoint1Any2    (const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPoint1PlatoX2 (const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPoint1PlatoY2 (const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+
    virtual void performPoint1Point2  (const double tmin1, const double tmax1, const double tmin2, const double tmax2);//+

  private:

    // Расстояние между кривыми класса Point и PlatoY.
    static
      const std::tuple<double, double, double>
      execPoint_and_PlatoY
      (
        const double tMinOfPoint,
        const double tMaxOfPoint,
        const baseCurve& curvePoint,
        const double tMinOfPlatoY,
        const double tMaxOfPlatoY,
        const baseCurve& curvePlatoY
      );

    // Расстояние между кривыми класса Point и PlatoX.
    static
      const std::tuple<double, double, double>
      execPoint_and_PlatoX
      (
        const double tMinOfPoint,
        const double tMaxOfPoint,
        const baseCurve& curvePoint,
        const double tMinOfPlatoX,
        const double tMaxOfPlatoX,
        const baseCurve& curvePlatoX
      );

    // Расстояние между кривыми класса Point и (Screen | Normal).
    static
      const std::tuple<double, double, double>
      execPoint_and_Any
      (
        const double tMinOfPoint,
        const double tMaxOfPoint,
        const baseCurve& curvePoint,
        const double tMinOfAny,
        const double tMaxOfAny,
        const baseCurve& curveAny
      );

    // Расстояние между кривыми класса PlatoX и PlatoY.
    static
      const std::tuple<double, double, double>
      execPlatoX_and_PlatoY
      (
        const double tminPlatoX,
        const double tmaxPlatoX,
        const baseCurve& curvePlatoX,

        const double tminPlatoY,
        const double tmaxPlatoY,
        const baseCurve& curvePlatoY
      );

    // Расстояние между кривыми класса PlatoX и (Screen | Normal).
    static
      const std::tuple<double, double, double>
      execPlatoX_and_Any
      (
        const double tminPlatoX,
        const double tmaxPlatoX,
        const baseCurve& curvePlatoX,

        const double tminAny,
        const double tmaxAny,
        const baseCurve& curveAny
      );

    // Расстояние между кривыми класса PlatoY и (Screen | Normal).
    static
      const std::tuple<double, double, double>
      execPlatoY_and_Any
      (
        const double tminPlatoY,
        const double tmaxPlatoY,
        const baseCurve& curvePlatoY,

        const double tminAny,
        const double tmaxAny,
        const baseCurve& curveAny
      );

    // Расстояние между кривыми класса:
    //  - Screen и Normal и наоборот.
    static
      const std::tuple<double, double, double>
      execAnyOfDifferentOrientation
      (
        const double tmin1,
        const double tmax1,
        const baseCurve& curve1,
        const double tmin2,
        const double tmax2,
        const baseCurve& curve2
      );

    // Расстояние между кривыми класса:
    //  - Screen и Screen.
    //  - Normal и Normal.
    static
      const std::tuple<double, double, double>
      execAnyOfSameOrientation
      (
        const double tmin1,
        const double tmax1,
        const baseCurve& curve1,
        const double tmin2,
        const double tmax2,
        const baseCurve& curve2
      );

    // Возвращает t curve точки пересечения кривой (на участке монотонности) и отрезка.
    // Данная функция используется для определения точки пересечения нормали к первой кривой с другой кривой.
    // Условие экстремума выполняется когда эта нормаль будет перпендикулярна другой кривой.
    static
      std::optional<double>
      getTofIntersectionWithSegment
      (
        const segmentCurve& segment,
        const double tmin,
        const double tmax,
        const baseCurve& curve
      );

    // Строим отрезок, перпендикулярный касательной к кривой в данной точке,
    // такой, чтобы его концы по оси X были бы [xmin; xmax] или по оси Y [ymin; ymax]
    // в зависимости от того, каков наклон перпендикуляра.
    static
      geom2d::segmentCurve
      buildSegmentCurveAsNormalToCurve
      (
        const double tOnCurve,
        const baseCurve& curve,
        const double xmin,
        const double xmax,
        const double ymin,
        const double ymax
      );

  private:
    double m_dist;
    double m_t1;
    double m_t2;
  };

}

