#include "pch.h"
#include "Eigen/Dense"

#include <iostream>

TEST(Eigen, Solve_BLAS_basic_example)
{
  // решаем СЛАУ 2x2
  // by fruitfull manual: https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
  using namespace Eigen;
  Matrix2d A;
  Vector2d right_part;
  Vector2d solution;
  A << 1, 2, 3, 4;
  right_part << 17, 41;

  solution = A.fullPivLu().solve(right_part);

  ASSERT_NEAR(solution(0), 7.0, 1.0e-9);
  ASSERT_NEAR(solution(1), 5.0, 1.0e-9);
}