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
  // Концепт класса для получения параметров кривой на участке монотонности, где одна из осей считается осью абсцисс, а другая - осью ординат
  // Причём в качестве оси абсцисс может фигурировать как ось X так и ось Y, и наоборот для оси ординат, в зависимости от потребности алгоритма,
  // который использует класс, удовлетворяющий данному концепту
  template <class T>
  concept CurveDataGetter = requires(const T a, const double t, const double theAny, const point & pnt)
  {
    // получить точку кривой для заданного параметра t
    { a.getPoint(t)     } -> std::convertible_to<geom2d::point>;

    // получить значение абсциссы кривой для заданного параметра t
    { a.getCoord(t)     } -> std::convertible_to<double>;

    // получить минимальное значение параметра t на заданном участке монотонности
    { a.getTmin()       } -> std::convertible_to<double>;

    // получить максимальное значение параметра t на заданном учестке монотонности
    { a.getTmax()       } -> std::convertible_to<double>;

    // получить минимальное значение (крайнее левое) значение абсциссы точки кривой на участке монотонности
    { a.getMinCoord()   } -> std::convertible_to<double>;
    
    // получить точку кривой, в которой достигается минимум абсциссы точки кривой на участке монотонности
    { a.getPointOfMin() } -> std::convertible_to<geom2d::point>;

    // получить значение параметра t, при котором значение абсциссы точки кривой достигает минимального значения на участке монотонности
    { a.getTofMin()     } -> std::convertible_to<double>;

    // получить максимальное значение (крайнее правое) значение абсциссы точки кривой на участке монотонности
    { a.getMaxCoord()   } -> std::convertible_to<double>;

    // получить точку кривой, в которой достигается максимум абсциссы точки кривой на участке монотонности
    { a.getPointOfMax() } -> std::convertible_to<geom2d::point>;

    // получить значение параметра t, при котором значение абсциссы точки кривой достигает максимального значения на участке монотонности
    { a.getTofMax()     } -> std::convertible_to<double>;

    { a.getTofCoord(theAny)  } -> std::convertible_to<double>; // the same as getTofAbscissa

    // получить значение параметра t, при котором значение абсциссы точки кривой достигает значения theAny
    { a.getTofAbscissa(theAny) } -> std::convertible_to<double>;

    // получить значение параметра t, при котором значение ординаты точки кривой достигает значения theAny
    { a.getTofOrdinate(theAny) } -> std::convertible_to<double>;

    // Получить абсциссу произвольной точки
    { T::abscissaOf(pnt) } -> std::convertible_to<double>;

    // Получить ординату произвольной точки
    { T::ordinateOf(pnt) } -> std::convertible_to<double>;
  };
  
  ///////////////////
  //               //
  // DataGetterOfX //
  //               //
  ///////////////////
  // Данный класс удовлетворяет концепту CurveDataGetter. Здесь ось абсцисс - X, а ось ординат - Y.
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
  // Данный класс удовлетворяет концепту CurveDataGetter. Здесь ось абсцисс - Y, а ось ординат - X.
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
    // Вычисляем левый общий конец ОДЗ по оси абсцисс.
    // Концептуально возвращает max(xmin1, xmin2)
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

    // Вычисляем левый общий конец ОДЗ по оси абсцисс.
    // Концептуально возвращает min(xmax1, xmax2)
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