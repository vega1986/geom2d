#pragma once
#include <tuple>
#include <concepts>

#include "point.h"
#include "baseCurve.h"
#include "StatOfCurvePiece.h"
#include "findFunctionRoots.h"

namespace geom2d
{
  // Тип для хранения параметров общего интервала ОДЗ по одной из осей координат
  // ОДЗ здесь (!) определяется наименьшим и наибольшим параметром кривой на плоскости,
  // а не ОДЗ в алгебраическом смысле
  
  // 2 параметр - значение параметра прикотором соответствущая координата кривой №1
  //              достигает заданного значения.
  // 3 параметр - точка на кривой №1, соответствущая этому параметру
  // 4 параметр - значение параметра прикотором соответствущая координата кривой №2
  //              достигает заданного значения.
  // 5 параметр - точка на кривой №2, соответствущая этому параметру
  using CommonRangeOfTwoCurves =
    std::tuple
    <
      double,        // значение общего минимума или максимума ОДЗ по X или Y для двух кривых.
      double,        // значение параметра прикотором соответствущая координата кривой №1 достигает заданного значения
      geom2d::point, // точка на кривой №1, соответствущая этому параметру
      double,        // значение параметра прикотором соответствущая координата кривой №1 достигает заданного значения
      geom2d::point  // точка на кривой №2, соответствущая этому параметру
    >;

  /////////////
  //         //
  // CONCEPT //
  //         //
  /////////////
  template <class T>
  concept CurveDataGetter = requires(const T a, const double t)
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
  };
  
  ///////////////////
  //               //
  // DataGetterOfX //
  //               //
  ///////////////////
  // класс типа CurveDataGetter - для получения параметров конца ОДЗ двух кривых по оси X
  class DataGetterOfX
  {
  public:
    DataGetterOfX(const double theTmin, const double theTmax, const baseCurve& theCurve)
      :
      tmin{ theTmin },
      tmax{ theTmax },
      curve{ theCurve },
      scp{theTmin, theCurve.getPoint(theTmin), theTmax, theCurve.getPoint(theTmax)}{}

    const point getPoint(const double t) const { return curve.getPoint(t); }
    const double getCoord(const double t) const { return curve.getPoint(t).x; }
    const double getTmin() const { return tmin; }
    const double getTmax() const { return tmax; }

    
    const double getMinCoord() const { return scp.pointOfxmin().x; }
    const point getPointOfMin() const { return scp.pointOfxmin(); }
    const double getTofMin() const { return scp.tOfxmin(); }


    const double getMaxCoord() const { return scp.pointOfxmax().x; }
    const point getPointOfMax() const { return scp.pointOfxmax(); }
    const double getTofMax() const { return scp.tOfxmax(); }

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
  // класс типа CurveDataGetter - для получения параметров конца ОДЗ двух кривых по оси Y
  class DataGetterOfY
  {
  public:
    DataGetterOfY(const double theTmin, const double theTmax, const baseCurve& theCurve)
      :
      tmin{ theTmin },
      tmax{ theTmax },
      curve{ theCurve },
      scp{ theTmin, theCurve.getPoint(theTmin), theTmax, theCurve.getPoint(theTmax) } {}

    const point getPoint(const double t) const { return curve.getPoint(t); }
    const double getCoord(const double t) const { return curve.getPoint(t).y; }
    const double getTmin() const { return tmin; }
    const double getTmax() const { return tmax; }


    const double getMinCoord() const { return scp.pointOfymin().y; }
    const point getPointOfMin() const { return scp.pointOfymin(); }
    const double getTofMin() const { return scp.tOfymin(); }


    const double getMaxCoord() const { return scp.pointOfymax().y; }
    const point getPointOfMax() const { return scp.pointOfymax(); }
    const double getTofMax() const { return scp.tOfymax(); }

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
    // Ищем параметры левого конца ОДЗ
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

    // Ищем параметры правого конца ОДЗ
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