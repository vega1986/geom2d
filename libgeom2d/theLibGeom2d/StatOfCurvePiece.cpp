#include <algorithm>

#include "StatOfCurvePiece.h"

//*************************************************************************************************

geom2d::StatOfCurvePiece::StatOfCurvePiece
(
  const double tp,
  const point p,
  const double tq,
  const point q
)
{
  m_tmin = std::min(tp, tq);
  m_pointOfTmin = (tp < tq) ? p : q;

  m_tmax = std::max(tp, tq);
  m_pointOfTmax = (tp < tq) ? q : p;

  m_tOfPointOfMinX = (p.x < q.x) ? tp : tq;
  m_pointOfMinX    = (p.x < q.x) ?  p :  q;

  m_tOfPointOfMaxX = (p.x < q.x) ? tq : tp;
  m_pointOfMaxX    = (p.x < q.x) ?  q :  p;

  m_tOfPointOfMinY = (p.y < q.y) ? tp : tq;
  m_pointOfMinY    = (p.y < q.y) ?  p :  q;

  m_tOfPointOfMaxY = (p.y < q.y) ? tq : tp;
  m_pointOfMaxY    = (p.y < q.y) ?  q :  p;
}

//*************************************************************************************************

const double geom2d::StatOfCurvePiece::tmin() const
{
  return m_tmin;
}

//*************************************************************************************************

const geom2d::point geom2d::StatOfCurvePiece::pointOftmin() const
{
  return m_pointOfTmin;
}

//*************************************************************************************************

const double geom2d::StatOfCurvePiece::tmax() const
{
  return m_tmax;
}

//*************************************************************************************************

const geom2d::point geom2d::StatOfCurvePiece::pointOftmax() const
{
  return m_pointOfTmax;
}

//*************************************************************************************************

const double geom2d::StatOfCurvePiece::tOfxmin() const
{
  return m_tOfPointOfMinX;
}

//*************************************************************************************************

const geom2d::point geom2d::StatOfCurvePiece::pointOfxmin() const
{
  return m_pointOfMinX;
}

//*************************************************************************************************

const double geom2d::StatOfCurvePiece::tOfxmax() const
{
  return m_tOfPointOfMaxX;
}

//*************************************************************************************************

const geom2d::point geom2d::StatOfCurvePiece::pointOfxmax() const
{
  return m_pointOfMaxX;
}

//*************************************************************************************************

const double geom2d::StatOfCurvePiece::tOfymin() const
{
  return m_tOfPointOfMinY;
}

//*************************************************************************************************

const geom2d::point geom2d::StatOfCurvePiece::pointOfymin() const
{
  return m_pointOfMinY;
}

//*************************************************************************************************

const double geom2d::StatOfCurvePiece::tOfymax() const
{
  return m_tOfPointOfMaxY;
}

//*************************************************************************************************

const geom2d::point geom2d::StatOfCurvePiece::pointOfymax() const
{
  return m_pointOfMaxY;
}

//*************************************************************************************************
