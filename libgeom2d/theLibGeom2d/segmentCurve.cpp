#include "segmentCurve.h"
#include "avector.h"
#include "amath.h"

#include "Eigen/Dense"

geom2d::segmentCurve::segmentCurve(const point first, const point second)
  :
  baseCurve(),
  m_startPoint(first),
  m_finishPoint(second)
{
  vector seg{ m_startPoint, m_finishPoint };
  if (seg.length() <= math::tolerance::tolPoint)
  {
    throw std::logic_error("segment length is too small");
  }
}

geom2d::point geom2d::segmentCurve::getPoint(double t) const
{
  const double alpha = 1.0 - t;
  const double beta = t;
  const point result = (m_startPoint * alpha) + (m_finishPoint * beta);
  return result;
}

geom2d::point geom2d::segmentCurve::getVelocity(double t) const
{
  const point result = m_finishPoint - m_startPoint;
  return result;
}

double geom2d::segmentCurve::parameterMin() const
{
  return 0.0;
}

double geom2d::segmentCurve::parameterMax() const
{
  return 1.0;
}

double geom2d::segmentCurve::nearestTo(const point p) const
{
  vector tau = { m_startPoint , m_finishPoint };
  const double tOfNearest = (vector{ m_startPoint, p }, tau) / (tau, tau);
  if (tOfNearest < 0.0)
  {
    return 0.0;
  }
  else if (tOfNearest > 1.0)
  {
    return 1.0;
  }
  else
  {
    return tOfNearest;
  }
}

geom2d::vector geom2d::segmentCurve::getTau() const
{
  vector unitTau{ m_startPoint, m_finishPoint };
  unitTau.normalize();
  return unitTau;
}

double
geom2d::segmentCurve::distanceTo(const point pnt) const
{
  const double tOfNearest = nearestTo(pnt);
  const point pntOnSegment = getPoint(tOfNearest);
  return point::distance(pnt, pntOnSegment);
}

const geom2d::point geom2d::segmentCurve::start() const
{
  return m_startPoint;
}

const geom2d::point geom2d::segmentCurve::finish() const
{
  return m_finishPoint;
}

std::optional<std::pair<double, double>>
  geom2d::segmentCurve::intersectionOf
  (
    segmentCurve seg1,
    segmentCurve seg2
  )
{
  const auto tau1 = seg1.getTau();
  const auto tau2 = seg2.getTau();

  const auto N = std::abs(tau1 ^ tau2); // z-компонента векторного произведения

  if (N <= 0.0)
  {
    // считаем, что отрезки параллельны
    vector betweenSegmentsDir{ seg1.start(), seg2.start() };
    betweenSegmentsDir.normalize();
    const auto M = tau1 ^ betweenSegmentsDir;
    if (M <= 0.0)
    {
      // считаем, что отрезки лежат на одной прямой

      if ((tau1, tau2) > 0.0)
      {
        const vector e1b2{ seg1.finish(), seg2.start() };
        const vector e2b1{ seg2.finish(), seg1.start() };

        const vector b1b2{ seg1.start(), seg2.start() };
        const vector e1e2{ seg1.finish(), seg2.finish() };

        // утилизируем граничные случаи
        if (e1b2.length() <= 0.0) return std::pair{ 1.0, 0.0 };
        if (e2b1.length() <= 0.0) return std::pair{ 0.0, 1.0 };

        if (b1b2.length() <= 0.0) return std::pair{ 0.0, 0.0 };
        if (e1e2.length() <= 0.0) return std::pair{ 1.0, 1.0 };

        const auto e1b2unit = e1b2.normalized();
        const auto e2b1unit = e2b1.normalized();
        const auto b1b2unit = b1b2.normalized();
        const auto e1e2unit = e1e2.normalized();

        // смотрим в одну сторону
        if (((e1b2, tau1) > 0.0) or ((e2b1, tau1) > 0.0))
        {
          return std::nullopt;
        }
        // тестируем точку b1
        {
          const bool b1_between_b2_and_e2 = ((b1b2unit, tau2) < 0.0) and ((e2b1unit, tau2) < 0.0);
          if (b1_between_b2_and_e2)
          {
            const double t1sol = 0.0;
            const double t2sol = b1b2.length() / seg2.length();
            return std::pair{ t1sol, t2sol };
          }
        }
        // тестируем точку e1
        {
          const bool e1_between_b2_and_e2 = ((e1b2unit, tau2) < 0.0) and ((e1e2unit, tau2) > 0.0);
          if (e1_between_b2_and_e2)
          {
            const double t1sol = 1.0;
            const double t2sol = e1b2.length() / seg2.length();
            return std::pair{ t1sol, t2sol };
          }
        }
        // тестируем точку b2
        {
          const bool b2_between_b1_and_e1 = ((b1b2unit, tau1) > 0.0) and ((e1b2unit, tau1) < 0.0);
          if (b2_between_b1_and_e1)
          {
            const double t2sol = 0.0;
            const double t1sol = b1b2.length() / seg1.length();
            return std::pair{ t1sol, t2sol };
          }
        }
        // тестируем точку e2
        {
          const bool e2_between_b1_and_e1 = ((e2b1unit, tau1) < 0.0) and ((e1e2unit, tau1) < 0.0);
          if (e2_between_b1_and_e1)
          {
            const double t2sol = 1.0;
            const double t1sol = e2b1.length() / seg1.length();
            return std::pair{ t1sol, t2sol };
          }
        }
        throw std::logic_error("impossible position in function geom2d::segmentCurve::intersectionOf");
      }
      else
      {
        const vector e1e2{ seg1.finish(), seg2.finish() };
        const vector b2b1{ seg2.start(), seg1.start() };

        const vector e1b2{ seg1.finish(), seg2.start() };
        const vector e2b1{ seg2.finish(), seg1.start() };

        // утилизируем граничные случаи
        if (e1e2.length() <= 0.0) return std::pair{ 1.0, 1.0 };
        if (b2b1.length() <= 0.0) return std::pair{ 0.0, 0.0 };

        if (e1b2.length() <= 0.0) return std::pair{ 1.0, 0.0 };
        if (e2b1.length() <= 0.0) return std::pair{ 0.0, 1.0 };

        const auto e1e2unit = e1e2.normalized();
        const auto b2b1unit = b2b1.normalized();
        const auto e1b2unit = e1b2.normalized();
        const auto e2b1unit = e2b1.normalized();

        // смотрим в разные стороны
        if (((e1e2, tau1) > 0.0) or ((b2b1, tau1) > 0.0))
        {
          return std::nullopt;
        }
        // тестируем точку b1
        {
          const bool b1_between_e2_and_b2 = ((b2b1unit, tau2) > 0.0) and ((e2b1unit, tau2) < 0.0);
          if (b1_between_e2_and_b2)
          {
            const double t1sol = 0.0;
            const double t2sol = b2b1.length() / seg2.length();
            return std::pair{ t1sol, t2sol };
          }
        }
        // тестируем точку e1
        {
          const bool e1_between_e2_and_b2 = ((e1b2unit, tau2) < 0.0) and ((e1e2unit, tau2) > 0.0);
          if (e1_between_e2_and_b2)
          {
            const double t1sol = 1.0;
            const double t2sol = e1b2.length() / seg2.length();
            return std::pair{ t1sol, t2sol };
          }
        }
        // тестируем точку b2
        {
          const bool b2_between_b1_and_e1 = ((b2b1unit, tau1) < 0.0) and ((e1b2unit, tau1) < 0.0);
          if (b2_between_b1_and_e1)
          {
            const double t2sol = 0.0;
            const double t1sol = b2b1.length() / seg1.length();
            return std::pair{ t1sol, t2sol };
          }
        }
        // тестируем точку e2
        {
          const bool e2_between_b1_and_e1 = ((e2b1unit, tau1) < 0.0) and ((e1e2unit, tau1) < 0.0);
          if (e2_between_b1_and_e1)
          {
            const double t2sol = 1.0;
            const double t1sol = e2b1.length() / seg1.length();
            return std::pair{ t1sol, t2sol };
          }
        }
        throw std::logic_error("impossible position in function geom2d::segmentCurve::intersectionOf");
      }
    }
    else
    {
      return std::nullopt;
    }
  }
  else
  {
    // решаем СЛАУ 2x2
    using namespace Eigen;
    Matrix2d A;
    Vector2d right_part;

    const auto b1 = seg1.start();
    const auto e1 = seg1.finish();

    const auto b2 = seg2.start();
    const auto e2 = seg2.finish();

    A << (e1 - b1).x, (b2 - e2).x, (e1 - b1).y, (b2 - e2).y;
    right_part << (b2 - b1).x, (b2 - b1).y;

    const Vector2d solution = A.fullPivLu().solve(right_part);
    const auto t1 = solution(0);
    const auto t2 = solution(1);
    return std::pair{ t1, t2 };
  }
  return std::nullopt;
}

std::tuple<double, double, double>
geom2d::segmentCurve::distanceBetween(segmentCurve seg1, segmentCurve seg2)
{
  const auto intersectionSolution = intersectionOf(seg1, seg2);
  if (intersectionSolution)
  {
    const auto [t1, t2] = intersectionSolution.value();
    if (t1 >= 0.0 and t1 <= 1.0 and t2 >= 0.0 and t2 <= 1.0)
    {
      return std::tuple{ 0.0, t1, t2 };
    }
  }

  // пересечения отрезков нет
  double dist = math::infinite::distance;
  double tex1 = 0.0;
  double tex2 = 0.0;
  {
    const auto t1 = seg1.nearestTo(seg2.start());
    const double curDist = seg1.distanceTo(seg2.start());
    if (curDist < dist)
    {
      dist = curDist;
      tex1 = t1;
      tex2 = seg2.parameterMin();
    }
  }
  {
    const auto t1 = seg1.nearestTo(seg2.finish());
    const double curDist = seg1.distanceTo(seg2.finish());
    if (curDist < dist)
    {
      dist = curDist;
      tex1 = t1;
      tex2 = seg2.parameterMax();
    }
  }
  {
    const auto t2 = seg2.nearestTo(seg1.start());
    const double curDist = seg2.distanceTo(seg1.start());
    if (curDist < dist)
    {
      dist = curDist;
      tex1 = seg1.parameterMin();
      tex2 = t2;
    }
  }
  {
    const auto t2 = seg2.nearestTo(seg1.finish());
    const double curDist = seg2.distanceTo(seg1.finish());
    if (curDist < dist)
    {
      dist = curDist;
      tex1 = seg1.parameterMax();
      tex2 = t2;
    }
  }
  return std::tuple{ dist, tex1, tex2 };
}

const double geom2d::segmentCurve::length() const
{
  return point::distance(m_startPoint, m_finishPoint);
}
