#pragma once
#include "point.h"
#include "theLibGeom2d.h"

namespace geom2d
{
  // Полезный класс помошник для анализа концов кривой на участке монотонности.
  class THELIBGEOM2D_API StatOfCurvePiece
  {
  public:
    StatOfCurvePiece(const double tp, const point p, const double tq, const point q);
    
    // минимальное значение параметра на участке монотонности
    const double tmin() const;
    
    // точка на кривой, которая соответствует минимальному значению параметра на участке монотонности
    const point pointOftmin() const;
    
    // максимальное значение параметра на участке монотонности
    const double tmax() const;
    
    // точка на кривой, которая соответствует максимальному значению параметра на участке монотонности
    const point pointOftmax() const;

    // значение параметра, которому соответствует точка на кривой, x-координата которой минимальна на участке монотонности
    const double tOfxmin() const;

    // точка на кривой, значение x-координаты которой минимально на заданном участке монотонности
    const point pointOfxmin() const;

    // значение параметра, которому соответствует точка на кривой, x-координата которой максимальна на участке монотонности
    const double tOfxmax() const;

    // точка на кривой, значение x-координаты которой максимальна на заданном участке монотонности
    const point pointOfxmax() const;

    // значение параметра, которому соответствует точка на кривой, y-координата которой минимальна на участке монотонности
    const double tOfymin() const;

    // точка на кривой, значение y-координаты которой минимально на заданном участке монотонности
    const point pointOfymin() const;

    // значение параметра, которому соответствует точка на кривой, y-координата которой максимальна на участке монотонности
    const double tOfymax() const;

    // точка на кривой, значение y-координаты которой максимальна на заданном участке монотонности
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
