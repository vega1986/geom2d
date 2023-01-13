#pragma once

namespace math
{
  size_t factorial(size_t n);

  // ����� ���������� �� n �� k
  size_t place(size_t n, size_t k);

  // ����� ��������� �� n �� k
  size_t comb(size_t n, size_t k);

  // ����������� ���������� ������� - �������� � 2 ���� ����� �����������, ��� std::pow
  // ��� ������ ���������
  double power(double value, size_t p);

  struct infinite
  {
    // ����������� ���������� ����� �������
    static constexpr double distance = 1.0e+30;
  };

  struct tolerance
  {
    // �������� �����
    static constexpr double tolPoint = 1.0e-9;

    // ����� ��������
    static constexpr double tolNumeric = 1.0e-9;

  };
}
