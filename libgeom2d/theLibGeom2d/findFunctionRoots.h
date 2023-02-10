#pragma once

#include <set>

#include "amath.h"

namespace math
{
  // ���� ������������ ������ ������� func �� �������� ������� [a; b]
  // ������� ��������� ���������� �� ���� �������
  // ����������� ���� �������: func(a) * func(b) < 0
  // ���� � ��������� �� tol

  struct needStopDummy
  {
    needStopDummy() = default;

    bool operator()(const double t1, const double t2) const
    {
      return false;
    }
  };

  template<class theFunc, class NeedStop>
  double findUniqueFunctionRootUsingCriteria(double a, double b, theFunc& func, const NeedStop& ns)
  {
    static double func_a = func(a);
    static double func_b = func(b);
    constexpr double funcTolerance = math::tolerance::tolPoint * 1.0e-2;
    if (std::abs(func_a) <= funcTolerance)
    {
      return a;
    }
    if (std::abs(func_b) <= funcTolerance)
    {
      return b;
    }
    while (not ns(a, b))
    {
      const double funcLeft = func(a);
      const double funcRight = func(b);
      const double c = 0.5 * (a + b);
      const double funcMiddle = func(c);

      // ���������
      if (std::abs(funcMiddle) <= funcTolerance)
      {
        return c;
      }

      if (funcMiddle > 0.0)
      {
        if (funcLeft < 0.0)
        {
          b = c;
          continue;
        }
        else
        {
          a = c;
          continue;
        }
      }
      else if (funcMiddle < 0.0)
      {
        if (funcLeft > 0.0)
        {
          b = c;
          continue;
        }
        else
        {
          a = c;
          continue;
        }
      }
      else
      {
        // � ������ ������ ������ ����� ���������� �������
        break;
      }
    }
    return 0.5 * (a + b);
  }

  template<class theFunc>
  double findUniqueFunctionRoot(double a, double b, theFunc& func, const double funcTolerance = tolerance::tolPoint * 0.1)
  {
    const double func_a = func(a);
    const double func_b = func(b);
    if (std::abs(func_a) <= funcTolerance)
    {
      return a;
    }
    if (std::abs(func_b) <= funcTolerance)
    {
      return b;
    }
    while (true)
    {
      const double funcLeft = func(a);
      const double funcRight = func(b);
      const double c = 0.5 * (a + b);
      const double funcMiddle = func(c);

      // ���������
      if (std::abs(funcMiddle) <= funcTolerance)
      {
        return c;
      }
     
      if (funcMiddle > 0.0)
      {
        if (funcLeft < 0.0)
        {
          b = c;
          continue;
        }
        else
        {
          a = c;
          continue;
        }
      }
      else if (funcMiddle < 0.0)
      {
        if (funcLeft > 0.0)
        {
          b = c;
          continue;
        }
        else
        {
          a = c;
          continue;
        }
      }
      else
      {
        // � ������ ������ ������ ����� ���������� �������
        break;
      }
    }
    return 0.5 * (a + b);   
  }

  // ���� ��� ����� ������� �� ��������
  // ������� ��� (a; b), �������� �����
  // � �������� �� � ��������� result
  template<class theFunc>
  std::set<double> findFunctionRoots(const double a, const double b, theFunc & func)
  {
    std::set<double> result;
    // ����� ������� ������������ �� 100 ���������� ������
    constexpr size_t segmentNumber = 100;

    // ������� �������� ��� x (a; b) ������� �� segmentNumber ���������� ��������
    
    // ���������� �������� ������ ����� ������� j
    // � ����������� �� ��� ������ [0; segmentNumber)
    auto argLeftInSegment = [a, b](size_t j) -> double
    {
      const double beta = static_cast<double>(j) / static_cast<double>(segmentNumber);
      const double alpha = 1.0 - beta;
      return (a * alpha) + (b * beta);
    };

    // ���������� �������� ������� ����� ������� j
    // � ����������� �� ��� ������ [0; segmentNumber)
    auto argRightInSegment = [a, b](size_t j) -> double
    {
      const double beta = static_cast<double>(j + 1) / static_cast<double>(segmentNumber);
      const double alpha = 1.0 - beta;
      return (a * alpha) + (b * beta);
    };

    // ������� ��� ���������� � findRoot
    //auto curveXderiv = [&curve](double t) -> double
    //{
    //  return curve.getPointDerivative(t).x;
    //};

    // ���� �� ������� ������� �������
    for (size_t j = 0; j < segmentNumber; ++j)
    {
      const bool lastSegment = (j == (segmentNumber - 1));

      // ����� ��������� [aj; bj]:
      const double aj = argLeftInSegment(j);
      const double bj = argRightInSegment(j);

      const double funcaj = func(aj);
      const double funcbj = func(bj);

      if (funcaj > 0.0)
      {
        if (funcbj > 0.0)
        {
          // ����������� ��� �� ������ �������
          continue;
        }
        else if (funcbj < 0.0)
        {
          // ���� ����������� �� ������ �������
          const double foundRoot = findUniqueFunctionRoot(aj, bj, func);
          result.insert(foundRoot);
          continue;
        }
        else
        {
          if (not lastSegment)
          {
            // bj - ������� ������
            result.insert(bj);
          }
          continue;
        }
      }
      else if (funcaj < 0.0)
      {
        if (funcbj < 0.0)
        {
          // ����������� ��� �� ������ �������
          continue;
        }
        else if (funcbj > 0.0)
        {
          // ���� ����������� �� ������ �������
          const double foundRoot = findUniqueFunctionRoot(aj, bj, func);
          result.insert(foundRoot);
          continue;
        }
        else
        {
          if (not lastSegment)
          {
            // bj - ������� ������
            result.insert(bj);
          }
          continue;
        }
      }
      else
      {
        // ������� ��������� �������� 0 ��� �������� ��������� aj
        // ������ �������� ������ ���� ���� ������������ �� ���������� ����
        // ���� ���������������, ���� ��� ������ ���, ��� ��� ����� ��� �� ���������������
        continue;
      }
    }
    return result;
  }
}
