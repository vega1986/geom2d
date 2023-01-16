#pragma once

#include <set>

#include "math.h"

namespace math
{
  // Ищем единственный корень функции func на заданном отрезке [a; b]
  // Функция считается монотонной на этом отрезке
  // предусловие этой функции: func(a) * func(b) < 0
  // ищем с точностью до tol
  template<class theFunc>
  double findUniqueFunctionRoot(double a, double b, theFunc& func, const double funcTolerance = tolerance::tolPoint * 0.1)
  {
    while (true)
    {
      const double funcLeft = func(a);
      const double funcRight = func(b);
      const double c = 0.5 * (a + b);
      const double funcMiddle = func(c);
      // если функция достаточно близка к 0, то считаем, что нашли решение
      // левая граница интервала
      if (std::abs(funcLeft) <= funcTolerance)
      {
        return a;
      }
      // правая граница интервала
      if (std::abs(funcRight) <= funcTolerance)
      {
        return b;
      }
      // серединка
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
        // в данном случае корень точно посередине отрезка
        break;
      }
    }
    return 0.5 * (a + b);   
  }

  // ищем все корни функции на заданном
  // отрезке ОДЗ (a; b), исключая концы
  // и помещаем их в результат result
  template<class theFunc>
  std::set<double> findFunctionRoots(const double a, const double b, theFunc & func)
  {
    std::set<double> result;
    // делим область отпределения на 100 одинаковых частей
    constexpr size_t segmentNumber = 100;

    // отрезок числовой оси x (a; b) поделен на segmentNumber одинаковых кусочков
    
    // возвращает значение левого конца отрезка j
    // в зависимости от его номера [0; segmentNumber)
    auto argLeftInSegment = [a, b](size_t j) -> double
    {
      const double beta = static_cast<double>(j) / static_cast<double>(segmentNumber);
      const double alpha = 1.0 - beta;
      return (a * alpha) + (b * beta);
    };

    // возвращает значение правого конца отрезка j
    // в зависимости от его номера [0; segmentNumber)
    auto argRightInSegment = [a, b](size_t j) -> double
    {
      const double beta = static_cast<double>(j + 1) / static_cast<double>(segmentNumber);
      const double alpha = 1.0 - beta;
      return (a * alpha) + (b * beta);
    };

    // адаптор для применения в findRoot
    //auto curveXderiv = [&curve](double t) -> double
    //{
    //  return curve.getPointDerivative(t).x;
    //};

    // цикл по каждому отрезку деления
    for (size_t j = 0; j < segmentNumber; ++j)
    {
      const bool lastSegment = (j == (segmentNumber - 1));

      // имеет отрезочек [aj; bj]:
      const double aj = argLeftInSegment(j);
      const double bj = argRightInSegment(j);

      const double funcaj = func(aj);
      const double funcbj = func(bj);

      if (funcaj > 0.0)
      {
        if (funcbj > 0.0)
        {
          // пересечений нет на данном отрезке
          continue;
        }
        else if (funcbj < 0.0)
        {
          // ищем пересечение на данном отрезке
          const double foundRoot = findUniqueFunctionRoot(aj, bj, func);
          result.insert(foundRoot);
          continue;
        }
        else
        {
          if (not lastSegment)
          {
            // bj - искомый корень
            result.insert(bj);
          }
          continue;
        }
      }
      else if (funcaj < 0.0)
      {
        if (funcbj < 0.0)
        {
          // пересечений нет на данном отрезке
          continue;
        }
        else if (funcbj > 0.0)
        {
          // ищем пересечение на данном отрезке
          const double foundRoot = findUniqueFunctionRoot(aj, bj, func);
          result.insert(foundRoot);
          continue;
        }
        else
        {
          if (not lastSegment)
          {
            // bj - искомый корень
            result.insert(bj);
          }
          continue;
        }
      }
      else
      {
        // функция принимает значение 0 при значении аргумента aj
        // данная ситуация должна была быть рассмотренна на предыдущем шаге
        // либо проигнорирована, если это первый шаг, так как концы ОДЗ не рассматриваются
        continue;
      }
    }
    return result;
  }
}
