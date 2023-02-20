#pragma once
#include "CommonRangeHelper.h"
#include "baseCurve.h"
#include "avector.h"

#include <algorithm>

namespace geom2d
{
  namespace extrema
  {
    namespace two_parallel_plato
    {
      // класс для расчёта расстояния между кривыми (PlatoX & PlatoX) | (PlatoY & PlatoY).
      template <CurveDataGetter dataGetter>
      class distanceCalculator
      {
      public:
        distanceCalculator
        (
          const double theTmin1,
          const double theTmax1,
          const baseCurve& theCurve1,
          const double theTmin2,
          const double theTmax2,
          const baseCurve& theCurve2
        )
          :
          tmin1(theTmin1),
          tmax1(theTmax1),
          curve1(theCurve1),
          tmin2(theTmin2),
          tmax2(theTmax2),
          curve2(theCurve2)
        {}

        std::tuple<double, double, double>
          execute() const
        {
          double dist = math::infinite::distance;
          double textrema1 = 0.0;
          double textrema2 = 0.0;

          dataGetter dg1(tmin1, tmax1, curve1);
          dataGetter dg2(tmin2, tmax2, curve2);

          const auto P1 = dg1.getPointOfMin();
          const auto P1a = dataGetter::abscissaOf(P1);
          const auto Q1 = dg1.getPointOfMax();
          const auto Q1a = dataGetter::abscissaOf(Q1);

          const auto P2 = dg2.getPointOfMin();
          const auto P2a = dataGetter::abscissaOf(P2);
          const auto Q2 = dg2.getPointOfMax();
          const auto Q2a = dataGetter::abscissaOf(Q2);

          if (Q1a <= P2a)
          {
            dist = point::distance(Q1, P2);
            textrema1 = dg1.getTofMax();
            textrema2 = dg2.getTofMin();
            return std::tuple{ dist, textrema1, textrema2 };
          }
          else if (Q2a <= P1a)
          {
            dist = point::distance(P1, Q2);
            textrema1 = dg1.getTofMin();
            textrema2 = dg2.getTofMax();
            return std::tuple{ dist, textrema1, textrema2 };
          }
          else
          {
            if ((P1a > P2a) and (P1a < Q2a)) // test P1
            {
              const auto ton2 = dg2.getTofAbscissa(P1a);
              const auto pon2 = dg2.getPoint(ton2);

              dist = point::distance(P1, pon2);
              textrema1 = dg1.getTofMin();
              textrema2 = ton2;
             
              return std::tuple{dist, textrema1, textrema2};
            }
            else if ((Q1a > P2a) and (Q1a < Q2a)) // test Q1
            {
              const auto ton2 = dg2.getTofAbscissa(Q1a);
              const auto pon2 = dg2.getPoint(ton2);

              dist = point::distance(Q1, pon2);
              textrema1 = dg1.getTofMax();
              textrema2 = ton2;

              return std::tuple{ dist, textrema1, textrema2 };
            }
            else if ((P2a > P1a) and (P2a < Q1a)) // test P2
            {
              const auto ton1 = dg1.getTofAbscissa(P2a);
              const auto pon1 = dg1.getPoint(ton1);

              dist = point::distance(P2, pon1);
              textrema1 = ton1;
              textrema2 = dg2.getTofMin();

              return std::tuple{ dist, textrema1, textrema2 };
            }
            else if ((Q2a > P1a) and (Q2a < Q1a)) // test Q2
            {
              const auto ton1 = dg1.getTofAbscissa(Q2a);
              const auto pon1 = dg1.getPoint(ton1);

              dist = point::distance(Q2, pon1);
              textrema1 = ton1;
              textrema2 = dg2.getTofMax();

              return std::tuple{ dist, textrema1, textrema2 };
            }
            else
            {
              throw std::logic_error("impossible situation in geom2d::extrema::two_parallel_plato::distanceCalculator::execute");
            }
          }
        }

      private:
        const double tmin1;
        const double tmax1;
        const baseCurve& curve1;

        const double tmin2;
        const double tmax2;
        const baseCurve& curve2;
      };
    }
  }
}