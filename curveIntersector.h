#pragma once

#include "baseCurve.h"
#include "point.h"
#include "IntersecctionSolutionType.h"

#include <set>
#include <vector>
#include <optional>


namespace geom2d
{
  // �������������� ��������� ������� ������ �� ������� ������������
  enum class curveClass : unsigned short int
  {
    Screen = 1, // (dx / dt > 0 && dy / dt < 0) || (dx / dt < 0 && dy / dt > 0)
    Normal = 2, // (dx / dt > 0 && dy / dt > 0) || (dx / dt < 0 && dy / dt < 0)
    PlatoX = 3, // (dx / dt == 0 && (dy / dt>0 || dy / dt < 0)
    PlatoY = 4, // (dy / dt == 0 && (dx / dt>0 || dx / dt < 0)
    Point  = 5  // (dx / dt == 0 && dy / dt == 0)
  };

  // ����� ������������ ��� ������ ����� ����������� ���� ������ ���� baseCurve,
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

    // ���� ����� �����������
    void perform();
  private:

    // ���� ����� ����������� �� �������� �������� ������������ ����� ������
    void perform(const double tmin1, const double tmax1, const double tmin2, const double tmax2);

    // ���� �������� ��������� ��� ������� ����������� �������� �����
    // �������� �� ��������������� ����� ����� �� ������������ ����
    static std::set<double> rootsOfCurveVelocity(const baseCurve& curve);

    ///////////////////
    //               //
    // ����� � ����� //
    //               //
    ///////////////////

    // ����������� ����������� ������, ����� ��� ��� �������� � �����
    static bool performPointXPoint(const point p, const point q);

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

    // ����������� ����������� ������, ����� ������ - �����, � ������ - ����� �� x (x ~ const)
    // �� ������� ����� ����������� ���� ���� x - �� ���������, �������, ���� �� y ������ ����
    // �� ���������
    static std::optional<double>
      performPointVSAnyAlongY(
        const point pnt,
        const double tmin,
        const double tmax,
        const baseCurve& curve);

    // ����������� ����������� ������, ����� ������ - �����, � ������ - ����� �� y (y ~ const)
    // �� ������� ����� ����������� ���� ���� y - �� ���������, �������, ���� �� x ������ ����
    // �� ���������
    static std::optional<double>
      performPointVSAnyAlongX(
        const point pnt,
        const double tmin,
        const double tmax,
        const baseCurve& curve);

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
    // � ������ �� �����, �.�. ��������� ��� �� X, ��� � �� Y.
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

    // ������� ����� ���� ������
    // � ������������� ��� ��� ����� ������������ 1 ���
    // (���� ���� ������ �������, � ������ ���������� ����� ��� X ��� Y �� ����� ������ ����� ������������ 1 ���)
    // ��� �� ��� ������� ����� ����������� ����� ������ ����� ����������� ���� ���������
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

    // ������� ����� ���� ������
    // � ������������� ��� ��� ����� ������������ 1 ���
    // (���� ���� ������ �������, � ������ ���������� ����� ��� X ��� Y �� ����� ������ ����� ������������ 1 ���)
    // ��� �� ��� ������� ����� ����������� ����� ������ ����� ����������� ���� ���������
    // ������ ������� ���������� �� ���������� ���, ��� �������� ������������ ������
    // ����� ������� ������� ������� ��� ������ ����� ����������� ����������� �� ������ ������ (�� �����������)
    // ���������� ������ ����������� ������ - ����� ��� X
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
    // ���������� ������ ����������� ������ - ����� ��� Y
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

    // ������� ���������� �������� ���������� ��� ������ (����������� �� ����) ������ ���� ��� �� X
    // I - commonXmin - max of minX of two curves
    // II - tof commonXmin of first curve
    // III - tof commonXmin of second curve
    static
      std::tuple<double, double, double>
      findLargestCommonXmin
      (
        const double tmin1,
        const double tmax1,
        const baseCurve & curve1,
        const double tmin2,
        const double tmax2,
        const baseCurve & curve2
      );


  private:

    const baseCurve & m_curve1;
    const baseCurve & m_curve2;

    std::vector<point> solutionPoints;
    std::vector<double> solutionParameterOnCurve1;
    std::vector<double> solutionParameterOnCurve2;
  };
}
