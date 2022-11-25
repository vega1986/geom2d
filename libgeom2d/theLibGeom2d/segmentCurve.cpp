#include "segmentCurve.h"
#include "avector.h"
#include "amath.h"

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

double geom2d::segmentCurve::distanceTo(const point pnt) const
{
  const double tOfNearest = nearestTo(pnt);
  const point pntOnSegment = getPoint(tOfNearest);
  return point::distance(pnt, pntOnSegment);
}
