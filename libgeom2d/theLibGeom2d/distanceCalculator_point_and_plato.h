#pragma once
//#include "theLibGeom2d.h"
#include "CommonRangeHelper.h"
#include "baseCurve.h"

namespace geom2d
{
  namespace extrema
  {
    namespace point_and_plato
    {
      template <CurveDataGetter dataGetter>
      class distanceCalculator
      {
      public:
        
        distanceCalculator
        (
          const point & thePnt,
          const double theTofPlatoMin,
          const double theTofPlatoMax,
          const baseCurve& theCurvePlato
        ):pnt(thePnt), tofPlatoMin(theTofPlatoMin), tofPlatoMax(theTofPlatoMax), curvePlato(theCurvePlato) {}

        std::pair<double, double>
          execute() const
        {
          double dist = math::infinite::distance;
          double tofExtremaOnPlato = 0.0;

          dataGetter dg{ tofPlatoMin , tofPlatoMax , curvePlato };
          const auto pnt_abscissa = dataGetter::abscissaOf(pnt);
          if (pnt_abscissa <= dg.getMinCoord())
          {
            const auto curdist = point::distance(pnt, dg.getPointOfMin());
            if (curdist < dist)
            {
              dist = curdist;
              tofExtremaOnPlato = dg.getTofMin();
            }
          }
          else if (pnt_abscissa >= dg.getMaxCoord())
          {
            const auto curdist = point::distance(pnt, dg.getPointOfMax());
            if (curdist < dist)
            {
              dist = curdist;
              tofExtremaOnPlato = dg.getTofMax();
            }
          }
          else
          {
            const auto tonPlato = dg.getTofAbscissa(pnt_abscissa);
            const auto pntOnPlato = dg.getPoint(tonPlato);
            const auto curdist = point::distance(pnt, pntOnPlato);
            if (curdist < dist)
            {
              dist = curdist;
              tofExtremaOnPlato = tonPlato;
            }
          }
          return std::pair{ dist, tofExtremaOnPlato };
        }

      private:
        const point pnt;

        const double tofPlatoMin;
        const double tofPlatoMax;
        const baseCurve& curvePlato;
      };
    }
  }
}
