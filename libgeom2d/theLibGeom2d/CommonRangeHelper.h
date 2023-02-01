#pragma once
#include <tuple>
#include <concepts>

#include "point.h"
#include "baseCurve.h"
#include "StatOfCurvePiece.h"
#include "findFunctionRoots.h"
#include "theLibGeom2d.h"

namespace geom2d
{
  // ��� ��� �������� ���������� ������ ��������� ��� �� ����� �� ���� ���������
  // ��� ����� (!) ������������ ���������� � ���������� ���������� ������ �� ���������,
  // � �� ��� � �������������� ������
  
  // 2 �������� - �������� ��������� ���������� �������������� ���������� ������ �1
  //              ��������� ��������� ��������.
  // 3 �������� - ����� �� ������ �1, �������������� ����� ���������
  // 4 �������� - �������� ��������� ���������� �������������� ���������� ������ �2
  //              ��������� ��������� ��������.
  // 5 �������� - ����� �� ������ �2, �������������� ����� ���������
  using CommonRangeOfTwoCurves =
    std::tuple
    <
      double,        // �������� ������ �������� ��� ��������� ��� �� X ��� Y ��� ���� ������.
      double,        // �������� ��������� ���������� �������������� ���������� ������ �1 ��������� ��������� ��������
      geom2d::point, // ����� �� ������ �1, �������������� ����� ���������
      double,        // �������� ��������� ���������� �������������� ���������� ������ �1 ��������� ��������� ��������
      geom2d::point  // ����� �� ������ �2, �������������� ����� ���������
    >;

  /////////////
  //         //
  // CONCEPT //
  //         //
  /////////////
  template <class T>
  concept CurveDataGetter = requires(const T a, const double t, const double theAny, const point & pnt)
  {
    { a.getPoint(t)     } -> std::convertible_to<geom2d::point>;
    { a.getCoord(t)     } -> std::convertible_to<double>;
    { a.getTmin()       } -> std::convertible_to<double>;
    { a.getTmax()       } -> std::convertible_to<double>;

    { a.getMinCoord()   } -> std::convertible_to<double>;
    { a.getPointOfMin() } -> std::convertible_to<geom2d::point>;
    { a.getTofMin()     } -> std::convertible_to<double>;

    { a.getMaxCoord()   } -> std::convertible_to<double>;
    { a.getPointOfMax() } -> std::convertible_to<geom2d::point>;
    { a.getTofMax()     } -> std::convertible_to<double>;

    { a.getTofCoord(theAny)  } -> std::convertible_to<double>; // the same as getTofAbscissa

    { a.getTofAbscissa(theAny) } -> std::convertible_to<double>;
    { a.getTofOrdinate(theAny) } -> std::convertible_to<double>;

    { T::abscissaOf(pnt) } -> std::convertible_to<double>;
    { T::ordinateOf(pnt) } -> std::convertible_to<double>;
  };
  
  ///////////////////
  //               //
  // DataGetterOfX //
  //               //
  ///////////////////
  // ����� ���� CurveDataGetter - ��� ��������� ���������� ����� ��� ���� ������ �� ��� X
  class DataGetterOfX
  {
  public:
    THELIBGEOM2D_API DataGetterOfX(const double theTmin, const double theTmax, const baseCurve& theCurve)
      :
      tmin{ theTmin },
      tmax{ theTmax },
      curve{ theCurve },
      scp{theTmin, theCurve.getPoint(theTmin), theTmax, theCurve.getPoint(theTmax)}{}

    THELIBGEOM2D_API const point getPoint(const double t) const { return curve.getPoint(t); }
    THELIBGEOM2D_API const double getCoord(const double t) const { return curve.getPoint(t).x; }
    THELIBGEOM2D_API const double getTmin() const { return tmin; }
    THELIBGEOM2D_API const double getTmax() const { return tmax; }

    
    THELIBGEOM2D_API const double getMinCoord() const { return scp.pointOfxmin().x; }
    THELIBGEOM2D_API const point getPointOfMin() const { return scp.pointOfxmin(); }
    THELIBGEOM2D_API const double getTofMin() const { return scp.tOfxmin(); }


    THELIBGEOM2D_API const double getMaxCoord() const { return scp.pointOfxmax().x; }
    THELIBGEOM2D_API const point getPointOfMax() const { return scp.pointOfxmax(); }
    THELIBGEOM2D_API const double getTofMax() const { return scp.tOfxmax(); }


    THELIBGEOM2D_API const double getTofCoord(const double theX) const { return curve.tofX(tmin, tmax, theX); }

    THELIBGEOM2D_API const double getTofAbscissa(const double theX) const { return curve.tofX(tmin, tmax, theX); }
    THELIBGEOM2D_API const double getTofOrdinate(const double theY) const { return curve.tofY(tmin, tmax, theY); }

    static inline const double abscissaOf(const point& pnt) { return pnt.x; }
    static inline const double ordinateOf(const point& pnt) { return pnt.y; }

  private:
    const double tmin;
    const double tmax;
    const baseCurve& curve;
    const StatOfCurvePiece scp;
  };

  ///////////////////
  //               //
  // DataGetterOfY //
  //               //
  ///////////////////
  // ����� ���� CurveDataGetter - ��� ��������� ���������� ����� ��� ���� ������ �� ��� Y
  class DataGetterOfY
  {
  public:
    THELIBGEOM2D_API DataGetterOfY(const double theTmin, const double theTmax, const baseCurve& theCurve)
      :
      tmin{ theTmin },
      tmax{ theTmax },
      curve{ theCurve },
      scp{ theTmin, theCurve.getPoint(theTmin), theTmax, theCurve.getPoint(theTmax) } {}

    THELIBGEOM2D_API const point getPoint(const double t) const { return curve.getPoint(t); }
    THELIBGEOM2D_API const double getCoord(const double t) const { return curve.getPoint(t).y; }
    THELIBGEOM2D_API const double getTmin() const { return tmin; }
    THELIBGEOM2D_API const double getTmax() const { return tmax; }


    THELIBGEOM2D_API const double getMinCoord() const { return scp.pointOfymin().y; }
    THELIBGEOM2D_API const point getPointOfMin() const { return scp.pointOfymin(); }
    THELIBGEOM2D_API const double getTofMin() const { return scp.tOfymin(); }


    THELIBGEOM2D_API const double getMaxCoord() const { return scp.pointOfymax().y; }
    THELIBGEOM2D_API const point getPointOfMax() const { return scp.pointOfymax(); }
    THELIBGEOM2D_API const double getTofMax() const { return scp.tOfymax(); }


    THELIBGEOM2D_API const double getTofCoord(const double theY) const { return curve.tofY(tmin, tmax, theY); }

    THELIBGEOM2D_API const double getTofAbscissa(const double theY) const { return curve.tofY(tmin, tmax, theY); }
    THELIBGEOM2D_API const double getTofOrdinate(const double theX) const { return curve.tofX(tmin, tmax, theX); }

    static inline const double abscissaOf(const point& pnt) { return pnt.y; }
    static inline const double ordinateOf(const point& pnt) { return pnt.x; }

  private:
    const double tmin;
    const double tmax;
    const baseCurve& curve;
    const StatOfCurvePiece scp;
  };

  ///////////////////////
  //                   //
  // CommonRangeHelper //
  //                   //
  ///////////////////////
  namespace CommonRangeHelper
  {
    // ���� ��������� ������ ����� ���
    template <CurveDataGetter dataGetter>
    CommonRangeOfTwoCurves
      ofLowest(const dataGetter& getter1, const dataGetter& getter2)
    {
      if (getter2.getMinCoord() < getter1.getMinCoord())
      {
        const auto commonMinCoord = getter1.getMinCoord();
        const auto tOfMinCoord1 = getter1.getTofMin();
        const auto pointOfMinCoord1 = getter1.getPointOfMin();

        auto func = [&getter2, value = commonMinCoord](const double t) -> double
        {
          return getter2.getCoord(t) - value;
        };
        // *
        const auto tOfMinCoord2 = math::findUniqueFunctionRoot(getter2.getTmin(), getter2.getTmax(), func);
        const auto pointOfMinCoord2 = getter2.getPoint(tOfMinCoord2);
        return std::tuple{ commonMinCoord, tOfMinCoord1, pointOfMinCoord1, tOfMinCoord2, pointOfMinCoord2 };
      }
      else if (getter1.getMinCoord() < getter2.getMinCoord())
      {
        const auto commonMinCoord = getter2.getMinCoord();
        
        const auto tOfMinCoord2 = getter2.getTofMin();
        const auto pointOfMinCoord2 = getter2.getPointOfMin();

        auto func = [&getter1, value = commonMinCoord](const double t) -> double
        {
          return getter1.getCoord(t) - value;
        };
        const auto tOfMinCoord1 = math::findUniqueFunctionRoot(getter1.getTmin(), getter1.getTmax(), func);
        const auto pointOfMinCoord1 = getter1.getPoint(tOfMinCoord1);
        return std::tuple{ commonMinCoord, tOfMinCoord1, pointOfMinCoord1, tOfMinCoord2, pointOfMinCoord2 };
      }
      else
      {
        const auto commonMinCoord = getter1.getMinCoord();
        
        const auto tOfMinCoord1 = getter1.getTofMin();
        const auto pointOfMinCoord1 = getter1.getPointOfMin();

        const auto tOfMinCoord2 = getter2.getTofMin();
        const auto pointOfMinCoord2 = getter2.getPointOfMin();
        return std::tuple{ commonMinCoord, tOfMinCoord1, pointOfMinCoord1, tOfMinCoord2, pointOfMinCoord2 };
      }
    }

    // ���� ��������� ������� ����� ���
    template <CurveDataGetter dataGetter>
    CommonRangeOfTwoCurves
      ofHighest(const dataGetter& getter1, const dataGetter& getter2)
    {
      if (getter2.getMaxCoord() < getter1.getMaxCoord())
      {
        const auto commonMaxCoord = getter2.getMaxCoord();
        const auto tOfMaxCoord2 = getter2.getTofMax();
        const auto pointOfMaxCoord2 = getter2.getPointOfMax();

        auto func = [&getter1, value = commonMaxCoord](const double t) -> double
        {
          return getter1.getCoord(t) - value;
        };
        // *
        const auto tOfMaxCoord1 = math::findUniqueFunctionRoot(getter1.getTmin(), getter1.getTmax(), func);
        const auto pointOfMaxCoord1 = getter1.getPoint(tOfMaxCoord1);
        return std::tuple{ commonMaxCoord, tOfMaxCoord1, pointOfMaxCoord1, tOfMaxCoord2, pointOfMaxCoord2 };
      }
      else if (getter1.getMaxCoord() < getter2.getMaxCoord())
      {
        const auto commonMaxCoord = getter1.getMaxCoord();

        const auto tOfMaxCoord1 = getter1.getTofMax();
        const auto pointOfMaxCoord1 = getter1.getPointOfMax();

        auto func = [&getter2, value = commonMaxCoord](const double t) -> double
        {
          return getter2.getCoord(t) - value;
        };
        const auto tOfMaxCoord2 = math::findUniqueFunctionRoot(getter2.getTmin(), getter2.getTmax(), func);
        const auto pointOfMaxCoord2 = getter2.getPoint(tOfMaxCoord2);
        return std::tuple{ commonMaxCoord, tOfMaxCoord1, pointOfMaxCoord1, tOfMaxCoord2, pointOfMaxCoord2 };
      }
      else
      {
        const auto commonMaxCoord = getter1.getMaxCoord();

        const auto tOfMaxCoord1 = getter1.getTofMax();
        const auto pointOfMaxCoord1 = getter1.getPointOfMax();

        const auto tOfMaxCoord2 = getter2.getTofMax();
        const auto pointOfMaxCoord2 = getter2.getPointOfMax();
        return std::tuple{ commonMaxCoord, tOfMaxCoord1, pointOfMaxCoord1, tOfMaxCoord2, pointOfMaxCoord2 };
      }
    }
  };

}