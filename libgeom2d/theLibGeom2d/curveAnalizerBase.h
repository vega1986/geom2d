#pragma once
#include "theLibGeom2d.h"
#include "baseCurve.h"

#include <set>

namespace geom2d
{

  class curveAnalizerBase
  {
  private:

    // ядро класса - исполняет своё предназначение для каждой пары отрезков монотонности
    void perform(const double tmin1, const double tmax1, const double tmin2, const double tmax2);
    
    // после выполнения основного алгоритма мы производим постпроцессинг - для каждого Derived класса свой
    // по умолчанию мы не делаем ничего!
    virtual void postProcessing();

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

    // ищем значения параметра при которых направление вдижения точки
    // меняется на противоположное вдоль одной из координатных осей
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
    // делим кривые на сегменты монотонности и для каждой пары ОДЗ монотонности выполняем пиватную абстрактную функцию
    // которая будет своя для каждого Derived класса
    void perform();

  };

}