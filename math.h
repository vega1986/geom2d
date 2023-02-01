#pragma once
#include "libgeom2d/theLibGeom2d/theLibGeom2d.h"

namespace math
{
  THELIBGEOM2D_API size_t factorial(size_t n);

  // ����� ���������� �� n �� k
  THELIBGEOM2D_API size_t place(size_t n, size_t k);

  // ����� ��������� �� n �� k
  THELIBGEOM2D_API size_t comb(size_t n, size_t k);

  // ����������� ���������� ������� - �������� � 2 ���� ����� �����������, ��� std::pow
  // ��� ������ ���������
  THELIBGEOM2D_API double power(double value, size_t p);

  struct THELIBGEOM2D_API infinite
  {
    // ����������� ���������� ����� �������
    static constexpr double distance = 1.0e+30;
  };

  struct THELIBGEOM2D_API tolerance
  {
    // �������� �����
    static constexpr double tolPoint = 1.0e-9;

    // ����� ��������
    static constexpr double tolNumeric = 1.0e-9;

  };
}
