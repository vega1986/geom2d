#pragma once
#include "theLibGeom2d.h"
#include "curveAnalizerBase.h"
#include "amath.h"
#include "segmentCurve.h"

#include <tuple>
#include <optional>

namespace geom2d
{
  class curveMutualDistanceCalculator final : private curveAnalizerBase
  {
  public:
    THELIBGEOM2D_API
      curveMutualDistanceCalculator(const baseCurve& curve1, const baseCurve& curve2)
      :curveAnalizerBase(curve1, curve2), m_dist{ math::infinite::distance }, m_t1{ 0.0 }, m_t2{0.0} {}

    THELIBGEOM2D_API void fulfill();

    THELIBGEOM2D_API
      const std::tuple<double, double, double> getExtrema() const;

  private:
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
    // в зависимости от того, каков наклон перпендикуляра
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

