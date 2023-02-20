#pragma once
#include "point.h"
#include "baseCurve.h"
#include "CommonRangeHelper.h"

#include <optional>

namespace geom2d
{
  namespace point_and_curve_alongaxis
  {
    // ����� ��� ������ ����� ����������� ������ Point � ������������ ������ (�� Point).
    // ������ ������ ���� ��������� ����� ���� �� ����� �� ���� ���������.
    // ����� �� ���������� �������� ����� � ���� ����� �� ������ ������� �� �� ��������.
    // ���� ���������� ����� ������� �� ����������� ���������, �� �������, ��� ����� ����� �� ������.
    template <CurveDataGetter dataGetter>
    class solver
    {
    public:
      solver(const point theP, const double theTmin, const double theTmax, const baseCurve& theCurve)
        :P{ theP }, tmin{ theTmin }, tmax{ theTmax }, curve{ theCurve } {}

      std::optional<double>
        execute() const
      {
        dataGetter theGetter{tmin, tmax, curve};
        if (dataGetter::abscissaOf(P) <= theGetter.getMinCoord())
        {
          if (point::isSame(P, theGetter.getPointOfMin()))
          {
            return theGetter.getTofMin();
          }
          else
          {
            return std::nullopt;
          }
        }
        else if (dataGetter::abscissaOf(P) >= theGetter.getMaxCoord())
        {
          if (point::isSame(P, theGetter.getPointOfMax()))
          {
            return theGetter.getTofMax();
          }
          else
          {
            return std::nullopt;
          }
        }
        else
        {
          const auto tofCurve = theGetter.getTofCoord(dataGetter::abscissaOf(P));
          const auto pofCurve = curve.getPoint(tofCurve);
          if (point::isSame(P, pofCurve))
          {
            return tofCurve;
          }
          else
          {
            return std::nullopt;
          }
        }
        return std::nullopt;
      }

    private:
      const point P;
      const double tmin;
      const double tmax;
      const baseCurve& curve;
    };
  }
}