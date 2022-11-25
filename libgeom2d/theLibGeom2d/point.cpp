#include "point.h"
#include "amath.h"

double geom2d::point::distance(const point& p, const point& q)
{
  const double dx = p.x - q.x;
  const double dy = p.y - q.y;
  return std::sqrt(dx * dx + dy * dy);
}

bool geom2d::point::isSame(const point& p, const point& q)
{
  return (distance(p, q) <= math::tolerance::tolPoint);
}

geom2d::point geom2d::operator+(const point& p, const point& q)
{
  return point{ p.x + q.x, p.y + q.y };
}

geom2d::point geom2d::operator-(const point& p, const point& q)
{
  return point{ p.x - q.x, p.y - q.y };
}

geom2d::point geom2d::operator*(const point& p, const double value)
{
  return point{ p.x * value, p.y * value };
}

geom2d::point geom2d::operator*(const double value, const point & p)
{
  return operator*(p, value);
}

geom2d::point geom2d::operator/(const point& p, const double value)
{
  return point{ p.x / value, p.y / value };
}
