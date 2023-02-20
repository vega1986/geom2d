#pragma once
#include "point.h"
#include "theLibGeom2d.h"

namespace geom2d
{
  // �������� ����� �������� ��� ������� ������ ������ �� ������� ������������.
  class THELIBGEOM2D_API StatOfCurvePiece
  {
  public:
    StatOfCurvePiece(const double tp, const point p, const double tq, const point q);
    
    // ����������� �������� ��������� �� ������� ������������
    const double tmin() const;
    
    // ����� �� ������, ������� ������������� ������������ �������� ��������� �� ������� ������������
    const point pointOftmin() const;
    
    // ������������ �������� ��������� �� ������� ������������
    const double tmax() const;
    
    // ����� �� ������, ������� ������������� ������������� �������� ��������� �� ������� ������������
    const point pointOftmax() const;

    // �������� ���������, �������� ������������� ����� �� ������, x-���������� ������� ���������� �� ������� ������������
    const double tOfxmin() const;

    // ����� �� ������, �������� x-���������� ������� ���������� �� �������� ������� ������������
    const point pointOfxmin() const;

    // �������� ���������, �������� ������������� ����� �� ������, x-���������� ������� ����������� �� ������� ������������
    const double tOfxmax() const;

    // ����� �� ������, �������� x-���������� ������� ����������� �� �������� ������� ������������
    const point pointOfxmax() const;

    // �������� ���������, �������� ������������� ����� �� ������, y-���������� ������� ���������� �� ������� ������������
    const double tOfymin() const;

    // ����� �� ������, �������� y-���������� ������� ���������� �� �������� ������� ������������
    const point pointOfymin() const;

    // �������� ���������, �������� ������������� ����� �� ������, y-���������� ������� ����������� �� ������� ������������
    const double tOfymax() const;

    // ����� �� ������, �������� y-���������� ������� ����������� �� �������� ������� ������������
    const point pointOfymax() const;

  private:
    double m_tmin = 0.0;
    point m_pointOfTmin;

    double m_tmax = 0.0;
    point m_pointOfTmax;

    double m_tOfPointOfMinX = 0.0;
    point m_pointOfMinX;

    double m_tOfPointOfMaxX = 0.0;
    point m_pointOfMaxX;

    double m_tOfPointOfMinY = 0.0;
    point m_pointOfMinY;

    double m_tOfPointOfMaxY = 0.0;   
    point m_pointOfMaxY;
  };

}
