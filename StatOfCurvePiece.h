#pragma once
#include "point.h"

namespace geom2d
{

  class StatOfCurvePiece
  {
  public:
    StatOfCurvePiece(const double tp, const point p, const double tq, const point q);
    
    const double tmin() const;
    
    const point pointOftmin() const;
    
    const double tmax() const;
    
    const point pointOftmax() const;

    const double tOfxmin() const;

    const point pointOfxmin() const;

    const double tOfxmax() const;

    const point pointOfxmax() const;

    const double tOfymin() const;

    const point pointOfymin() const;

    const double tOfymax() const;

    const point pointOfymax() const;

  private:
    double m_tmin = 0.0;
    point m_pointOfTmin;

    double m_tmax = 0.0;
    point m_pointOfTmax;

    double m_tOfPointOfMinX = 0.0;
    point m_pointOfMinX;

    double m_tOfPointOfMaxX = 0.0;
    point m_pointOfMaxX;

    double m_tOfPointOfMinY = 0.0;
    point m_pointOfMinY;

    double m_tOfPointOfMaxY = 0.0;   
    point m_pointOfMaxY;
  };

}
