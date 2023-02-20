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
  // ����� ������������ ��� ������ ����� ����������� ���� ������ ���� baseCurve,
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

    // ��������� ������ �� ����������� - ��������� ������ ����� ����������� � �������� ����������
    THELIBGEOM2D_API void fulfill();

  private:
    
    // �������������� - �������� excludeDuplicatesFromSolution
    virtual void postProcessing();

    // ���������� �� ������� ������� ������������� �����
    void excludeDuplicatesFromSolution();

    // ��. ������� �����
    // ��������� ������� ���� ���������� ����� ����������� ��������� ������ �� �������� ������������
    // ����� ����������� � ��������������� ��������� ����������� � ������ �������
    // ����� ���� ������� ������������ � ������� ������, ����� ��� ������ �����������
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
    // ����� � ����� //
    //               //
    ///////////////////

    // ��������� �������, ���� ������ � ������ ������ - �����
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
    // ����� � �� ����� //
    //                  //
    //////////////////////

    // ��������� �������, ���� ������ ������ - �����,
    // � ������ 'X = const'
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

    // ��������� �������, ���� ������ ������ - �����,
    // � ������ 'Y = const'
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

    // ��������� �������, ���� ������ ������ - �����,
    // � ������ - ���� Normal ���� Screen
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

    // ��������� �������, ���� ������ � ������ ������ - PlatoY
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

    // ��������� �������, ���� ������ � ������ ������ - PlatoX
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

    // ��������� �������, ���� ������ ������ - PlatoX, ������ - PlatoY
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

    // ������� ����� ����������� ���� ������ �� �������� ������������
    // ���� ���� ������ �������, � ������ ���������� ����� ��� X ��� Y.
    // ��� �� ��� ������� ����� ����������� ����� ������ ����� ����������� ������ ���� ���������.
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

    // ������� ������������ ����� ����������� ���� ������ �� �������� ������������, ������� ����������� �� ��� ����� ��� X � �������� �������� ��������� ������ ������ �� ������ ���.
    // ��� ���������� ����� �����������, ��� �� t otherCurve ������� ������ �� ������ ���� ����� �� ���������� �������� y-���������� ����� �� ������ �� ���������.
    // ����� ������� ������� ������� ��� ������ ����� ����������� ����������� �� ������ ������ (�� �����������).
    // ���������� ������ AABB ����������� ������ - ����� ��� X.
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

    // ������� ������������ ����� ����������� ���� ������ �� �������� ������������, ������� ����������� �� ��� ����� ��� Y � �������� �������� ��������� ������ ������ �� ������ ���.
    // ��� ���������� ����� �����������, ��� �� t otherCurve ������� ������ �� ������ ���� ����� �� ���������� �������� x-���������� ����� �� ������ �� ���������.
    // ����� ������� ������� ������� ��� ������ ����� ����������� ����������� �� ������ ������ (�� �����������).
    // ���������� ������ AABB ����������� ������ - ����� ��� Y.
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

    // ��������� �������, ���� ������ ������ - PlatoX, ������ - Normal ��� Screen.
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

    // ��������� �������, ���� ������ ������ - PlatoY, ������ - Normal ��� Screen.
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

    // ��������� �������, ���� ������ ������ - Normal, ������ - Screen
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

    // ��������� �������, ���� ������ � ������ ������ - Normal
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

    // ��������� �������, ���� ������ � ������ ������ - Screen
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

    // ������ ������� � �������� �����
    THELIBGEOM2D_API void dumpIntersections(std::ostream& ost) const;

    // ���������� ������ ����� �����������
    THELIBGEOM2D_API std::vector<point> getSolutionPoints() const
    {
      return solutionPoints;
    }

    // ���������� ������ ���������� �� ������ 1 ����� �����������
    THELIBGEOM2D_API std::vector<double> getSolutionT1() const
    {
      return solutionParameterOnCurve1;
    }

    // ���������� ������ ���������� �� ������ 2 ����� �����������
    THELIBGEOM2D_API std::vector<double> getSolutionT2() const
    {
      return solutionParameterOnCurve2;
    }

  public:
    
    // ��������� ������� ��� ���� ������ �� �������� ������������ ���� �������� ��������, ��� ����������� �����������.
    // ������� ��������� ��� ������:
    //  - Screen & Normal � ��������.
    //  - PlatoX & (Screen | Normal) � ��������.
    //  - PlatoY & (Screen | Normal) � ��������.
    //  - PlatoX & PlatoY � ��������.
    // ������� �������� findUniqueIntersection.
    // ������ ������� ���������� ��� ��������� ������ ��������� ����� ���� ������:
    // � � ������� ������ ����� ����������� ������� � ����� ������ �� ������ ������.
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
