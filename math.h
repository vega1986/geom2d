#pragma once

namespace math
{
  size_t factorial(size_t n);

  // число размещений из n по k
  size_t place(size_t n, size_t k);

  // число сочетаний из n по k
  size_t comb(size_t n, size_t k);

  // оптимальное вычисление степени - примерно в 2 раза более эффективнее, чем std::pow
  // нет смысла применять
  double power(double value, size_t p);

  struct infinite
  {
    // бесконечное расстояние между точками
    static constexpr double distance = 1.0e+30;
  };

  struct tolerance
  {
    // толеранс точки
    static constexpr double tolPoint = 1.0e-9;

    // общий толеранс
    static constexpr double tolNumeric = 1.0e-9;

  };
}
