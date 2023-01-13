#pragma once

#include "baseCurve.h"
#include "point.h"

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

    // ��������� �������, ���� ������ ������ - �����,
    // � ������ 'Y = const'
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

    // ��������� �������, ���� ������ ������ - �����,
    // � ������ �� �����, �.�. ��������� ��� �� X, ��� � �� Y.
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

    // ��������� �������, ���� ������ � ������ ������ - PlatoY
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

    // ��������� �������, ���� ������ � ������ ������ - PlatoX
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
