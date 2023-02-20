#pragma once
#include <cmath>
#include <stdexcept>

#include "amath.h"
#include "theLibGeom2d.h"

namespace geom2d
{
  // Axis aligned boubding box
  class THELIBGEOM2D_API aabb
  {
    double m_xmin = 0.0;
    double m_xmax = 0.0;
    double m_ymin = 0.0;
    double m_ymax = 0.0;
  public:
    void clear()
    {
      m_xmin = 0.0;
      m_xmax = 0.0;
      m_ymin = 0.0;
      m_ymax = 0.0;
    }

    // Вычисляем границы AABB для произвольного непустого контейнера точек
    template <class pointsContainer>
    void reset(const pointsContainer& points)
    {
      if (points.end() == points.begin())
      {
        throw std::logic_error("void container of points in geom2d::aabb constructor");
      }
      m_xmin = points.begin()->x;
      m_ymin = points.begin()->y;

      m_xmax = points.begin()->x;
      m_ymax = points.begin()->y;

      for (const auto pnt : points)
      {
        // along X
        if (m_xmin > pnt.x)
        {
          m_xmin = pnt.x;
        }
        else if (m_xmax < pnt.x)
        {
          m_xmax = pnt.x;
        }
        else
        {
          // nothing to do
        }
        // along Y
        if (m_ymin > pnt.y)
        {
          m_ymin = pnt.y;
        }
        else if (m_ymax < pnt.y)
        {
          m_ymax = pnt.y;
        }
        else
        {
          // nothing to do
        }
      }
    }

    double xmin() const
    {
      return m_xmin;
    }

    double xmax() const
    {
      return m_xmax;
    }

    double ymin() const
    {
      return m_ymin;
    }

    double ymax() const
    {
      return m_ymax;
    }

    double lengthx() const
    {
      return m_xmax - m_xmin;
    }

    double lengthy() const
    {
      return m_ymax - m_ymin;
    }

    template <class pointsContainer>
    aabb(const pointsContainer & points)
    {
      reset(points);
    }
    
  };

}
