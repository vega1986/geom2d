#include <stdexcept>

#include "avector.h"
#include "amath.h"

geom2d::vector::vector(const double theX, const double theY)
{
  x = theX;
  y = theY;
}

geom2d::vector::vector(const point p, const point q)
{
  x = q.x - p.x;
  y = q.y - p.y;
}

geom2d::vector::vector(const point p)
{
  x = p.x;
  y = p.y;
}

double geom2d::vector::length() const
{
  return std::sqrt((x * x) + (y * y));
}

void geom2d::vector::normalize()
{
  const auto len = length();
  if (len <= math::tolerance::tolPoint)
  {
    throw std::logic_error("length of vector is 0");
  }
  x /= len;
  y /= len;
}

geom2d::vector geom2d::operator+(const vector p, const vector q)
{
  return vector{ p.x + q.x, p.y + q.y };
}

geom2d::vector geom2d::operator-(const vector p, const vector q)
{
  return vector{ p.x - q.x, p.y - q.y };
}

geom2d::vector geom2d::operator*(const vector p, const double value)
{
  return vector{ p.x * value, p.y * value };
}

geom2d::vector geom2d::operator*(const double value, const vector p)
{
  return operator*(p, value);
}

geom2d::vector geom2d::operator/(const vector p, const double value)
{
  return vector{p.x / value, p.y / value};
}

double geom2d::operator,(const vector p, const vector q)
{
  return p.x * q.x + p.y * q.y;
}

double geom2d::operator^(const vector p, const vector q)
{
  return p.x * q.y - p.y * q.x;
}
