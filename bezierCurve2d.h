#pragma once

#include "baseCurve.h"
#include "point.h"

#include <vector>

namespace geom2d
{

  class bezierCurve : public baseCurve
  {
  private:

    // ��������������� �������,
    // ���������� �� ������������
    void initialize_collateral ();

  public:

    template<class pointsContainer>
    bezierCurve (const pointsContainer& container)
      : baseCurve()
    {
      for (const auto p : container)
      {
        points.push_back(p);
      }
      initialize_collateral();
    }

    // �������� ����� �� ������ ��� ��������� ��������� t
    // �������������� �������: sum(j from 0 to n - 1) { beta(j) * t^j * (1 - t)^(n - 1 - j) }
    virtual point getPoint (double t) const;

    // �������� ����������� �� t �� ����� �� ������ ��� ��������� ��������� t
    // �������������� ��������: sum(j from 0 to n - 1)
    // { j * beta(j) * t^(j - 1) * (1 - t)^(n - 1 - j) - 
    //  - (n - 1 - j) * beta(j) * t^j * (1 - t)^(n - 2 - j) }
    virtual point getVelocity (double t) const;

    // ����������� �������� ��������� t
    virtual double parameterMin() const
    {
      return 0.0;
    }

    // ������������ �������� ��������� t
    virtual double parameterMax() const
    {
      return 1.0;
    }
    
  private:
    std::vector<point> points;

    // n         : ���������� ����� �������, j = [0; n - 1]
    // x[j]      : x-���������� j-�� ����� �������
    // y[j]      : �-���������� j-�� ����� �������
    // C(n, k)   : ����� ��������� �� n �� k
    // beta[j].x : x[j] * C(n - 1, j)
    // beta[j].y : y[j] * C(n - 1, j)
    std::vector<point> beta;

    // gamma[j]  : beta[j] * j
    std::vector<point> gamma;

    // ksi[j]    : beta[j] * (n - 1 - j)
    std::vector<point> ksi;
  };

}
