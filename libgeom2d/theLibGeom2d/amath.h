#pragma once
#include "theLibGeom2d.h"

namespace math
{
  THELIBGEOM2D_API size_t factorial(size_t n);

  // число размещений из n по k
  THELIBGEOM2D_API size_t place(size_t n, size_t k);

  // число сочетаний из n по k
  THELIBGEOM2D_API size_t comb(size_t n, size_t k);

  // оптимальное вычисление степени - примерно в 2 раза более эффективнее, чем std::pow
  // нет смысла применять
  THELIBGEOM2D_API double power(double value, size_t p);

  struct THELIBGEOM2D_API infinite
  {
    // бесконечное расстояние между точками
    static constexpr double distance = 1.0e+30;
  };

  struct THELIBGEOM2D_API tolerance
  {
    // толеранс точки
    static constexpr double tolPoint = 1.0e-9;

    // общий толеранс
    static constexpr double tolNumeric = 1.0e-9;

    // толеранс скалярного произведения единичных векторов
    static constexpr double tolScalarMult = 1.0e-7;

    // толеранс параметра
    static constexpr double tolParameter = 1.0e-10;

  };
}
