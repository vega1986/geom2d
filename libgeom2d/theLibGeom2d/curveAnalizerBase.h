#pragma once
#include "theLibGeom2d.h"
#include "baseCurve.h"

#include <set>

namespace geom2d
{

  class curveAnalizerBase
  {
  private:

    // ������ �� 2 ������ �� �������� [tmin1, tmax1] � [tmin2, tmax2] ������� ������ ���� ���������
    // ����������� ���� ������ �� �������� �������� ������������
    // ���������� ����� ������ ������ � � ����������� �� ��������� ���� ������� ��������
    // ��������������� �����.
    // ����� ���� 5 ������� ������ �� �������� ������������.
    // ������� ����� 25 ��������� ���������.
    // ��������� ��� ��������� ���� ������� ������ ����� ��������������� � ����� �������.
    // ����� ���� ���� ������ ������ (���� ������) �������� ������� - ���� ��������� �� ��������,
    // ������� � ������ ������ ���������� ������������ ��� ����������� �������, � ������� ���������� ���� ������ - � ������ ������� (1 � 2 ��� 2 � 1).
    void perform(const double tmin1, const double tmax1, const double tmin2, const double tmax2);
    
    // ����� ���������� ��������� ��������� �� ���������� �������������� - ��� ������� Derived ������ ����
    // �� ��������� �� �� ������ ������!
    virtual void postProcessing();

    // ���� ��� ��� �������. ������ ������ ���� ������������ � ���������� ������� ������������.
    // ������� ���������� ���������: perform + 'curveClass of 1 curve' + 1 + 'curve class of 2 curve' + 2
    // ���� � �������� 'curve class ...' ����� Any, �� ��� ������ ��� ������� ���������� � ��� Screen � ��� Normal
    // (�������� ��� ����� ������ ������, ����������� � ������ �� ���� �������).

    virtual void performScreen1Screen2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performScreen1Normal2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performScreen1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performScreen1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performScreen1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performNormal1Screen2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performNormal1Normal2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performNormal1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performNormal1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performNormal1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoX1Any2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoX1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoX1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoX1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoY1Any2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoY1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoY1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoY1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPoint1Any2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPoint1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPoint1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPoint1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    // ���� �������� ��������� ��� ������� ����������� �������� �����
    // �������� �� ��������������� ����� ����� �� ������������ ����
    static
      std::set<double>
        rootsOfCurveVelocity(const baseCurve& curve);

  protected:
    
    const baseCurve& m_curve1;
    const baseCurve& m_curve2;

    curveAnalizerBase(const baseCurve& curve1, const baseCurve& curve2)
      : m_curve1(curve1), m_curve2(curve2) {}

    ~curveAnalizerBase() = default; // ���������� �� �����������, ��� ��� ������������ �� ����� ������ ����� private!

  protected:
    // ����� ������ �� �������� ������������ � ��� ������ ���� ��� ������������ ��������� ��������� ����������� �������
    // ������� ����� ���� ��� ������� Derived ������
    void perform();

  };

}