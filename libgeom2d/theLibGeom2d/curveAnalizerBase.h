#pragma once
#include "theLibGeom2d.h"
#include "baseCurve.h"

#include <set>

namespace geom2d
{

  class curveAnalizerBase
  {
  private:

    // кажда€ из 2 кривых на заданных [tmin1, tmax1] и [tmin2, tmax2] отрезке должна быть монотонна
    // анализируем пару кривых на заданных участках монотонности
    // определ€ем класс каждой кривой и в зависимости от сочетани€ этих классов вызываем
    // соответствующий метод.
    // всего есть 5 классов кривых на отрезках монотонности.
    // ѕоэтому всего 25 различных вариантов.
    // Ќесколько пар сочетаний двух классов кривых могут рассматриватьс€ в одной функции.
    //  роме того если классы кривых (сами кривые) помен€ть местами - суть алгоритма не мен€етс€,
    // поэтому в каждом классе наследнике присутствуют р€д статических методов, в которых передаютс€ сами кривые - в нужном пор€дке (1 и 2 или 2 и 1).
    void perform(const double tmin1, const double tmax1, const double tmin2, const double tmax2);
    
    // после выполнени€ основного алгоритма мы производим постпроцессинг - дл€ каждого Derived класса свой
    // по умолчанию мы не делаем ничего!
    virtual void postProcessing();

    // Ќиже идЄт р€д методов. ћетоды должны быть перезаписаны в конкретных классах анализаторах.
    // ѕринцип именовани€ следующий: perform + 'curveClass of 1 curve' + 1 + 'curve class of 2 curve' + 2
    // ≈сли в качестве 'curve class ...' стоит Any, то это значит что функци€ вызываетс€ и дл€ Screen и дл€ Normal
    // (подходит дл€ любой второй кривой, относ€щейс€ к одному из этих классов).

    virtual void performScreen1Screen2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performScreen1Normal2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performScreen1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performScreen1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performScreen1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performNormal1Screen2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performNormal1Normal2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performNormal1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performNormal1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performNormal1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoX1Any2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoX1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoX1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoX1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoY1Any2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoY1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoY1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPlatoY1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPoint1Any2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPoint1PlatoX2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPoint1PlatoY2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    virtual void performPoint1Point2(const double tmin1, const double tmax1, const double tmin2, const double tmax2) abstract;

    // ищем значени€ параметра при которых направление вдижени€ точки
    // мен€етс€ на противоположное вдоль одной из координатных осей
    static
      std::set<double>
        rootsOfCurveVelocity(const baseCurve& curve);

  protected:
    
    const baseCurve& m_curve1;
    const baseCurve& m_curve2;

    curveAnalizerBase(const baseCurve& curve1, const baseCurve& curve2)
      : m_curve1(curve1), m_curve2(curve2) {}

    ~curveAnalizerBase() = default; // деструктор не виртуальный, так как наследование от этого класса будет private!

  protected:
    // делим кривые на сегменты монотонности и дл€ каждой пары ќƒ« монотонности выполн€ем приватную абстрактную функцию
    // котора€ будет сво€ дл€ каждого Derived класса
    void perform();

  };

}