#include <stdexcept>
#include "math.h"

//-----------------------------------------------------------------------------

size_t math::factorial(size_t n)
{
  size_t result = 1;
  if (n == 0)
  {
    return 1;
  }
  for (size_t i = n; i > 1; --i)
  {
    result *= i;
  }

  return result;
}

//-----------------------------------------------------------------------------

size_t math::place(size_t n, size_t k)
{
  if (k > n)
  {
    throw std::logic_error("число размещений из n по k: k > n");
  }
  if (n == 0)
  {
    return 1;
  }
  size_t result = 1;
  for (size_t i = n; i > n - k; --i)
  {
    result *= i;
  }
  return result;
}

//-----------------------------------------------------------------------------

size_t math::comb(size_t n, size_t k)
{
  if (n == k)
  {
    return 1;
  }
  // для эффективности вычислений
  if (k > (n / 2))
  {
    k = n - k;
  }
  return place(n, k) / factorial(k);
}

//-----------------------------------------------------------------------------

double math::power(double value, size_t p)
{
  if (p == 0)
  {
    return 1.0;
  }
  else if (p == 1)
  {
    return value;
  }
  else
  {
    if (p % 2 == 0)
    {
      p /= 2;
      const double temp = power(value, p);
      return temp * temp;
    }
    else
    {
      p = (p - 1) / 2;
      const double temp = power(value, p);
      return value * temp * temp;
    }
  }
}
