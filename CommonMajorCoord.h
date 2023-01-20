#pragma once
#include <tuple>
#include "baseCurve.h"
#include "Axis.h"
#include "StatOfCurvePiece.h"

namespace geom2d
{
  class CurveMajorHelper
  {
  public:
    CurveMajorHelper(const double tlow, const double thigh, const baseCurve& crv)
      :tmin(tlow), tmax(thigh), curve(crv), scp{ tlow, crv.getPoint(tlow), thigh, crv.getPoint(thigh) } {}

    virtual const double getCoord(const double t) const abstract;

    virtual const double getCoordMin() const abstract;
    virtual const point getPointMin() const abstract;
    virtual const double getTofCoordMin() const abstract;

    virtual const double getCoordMax() const abstract;
    virtual const point getPointMax() const abstract;
    virtual const double getTofCoordMax() const abstract;
  protected:
    const baseCurve& curve;
    const double tmin;
    const double tmax;
    const StatOfCurvePiece scp;
  };

  class XofCurvePointGetter : public CurveMajorHelper
  {
  public:
    XofCurvePointGetter(const double tlow, const double thigh, const baseCurve& crv) :CurveMajorHelper(tlow, thigh, crv) {}

    virtual const double getCoord(const double t) const override
    {
      return curve.getPoint(t).x;
    }

    virtual const double getCoordMin() const override
    {
      return scp.pointOfxmin().x;
    }

    virtual const point getPointMin() const override
    {
      return scp.pointOfxmin();
    }

    virtual const double getTofCoordMin() const override
    {
      return scp.tOfxmin();
    }

    virtual const double getCoordMax() const override
    {
      return scp.pointOfxmax().x;
    }

    virtual const point getPointMax() const override
    {
      return scp.pointOfxmax();
    }

    virtual const double getTofCoordMax() const override
    {
      return scp.tOfxmax();
    }
  };

  class YofCurvePointGetter : public CurveMajorHelper
  {
  public:
    YofCurvePointGetter(const double tlow, const double thigh, const baseCurve& crv) :CurveMajorHelper(tlow, thigh, crv) {}

    virtual const double getCoord(const double t) const override
    {
      return curve.getPoint(t).y;
    }

    virtual const double getCoordMin() const override
    {
      return scp.pointOfymin().y;
    }

    virtual const point getPointMin() const override
    {
      return scp.pointOfymin();
    }

    virtual const double getTofCoordMin() const override
    {
      return scp.tOfymin();
    }

    virtual const double getCoordMax() const override
    {
      return scp.pointOfymax().y;
    }

    virtual const point getPointMax() const override
    {
      return scp.pointOfymax();
    }

    virtual const double getTofCoordMax() const override
    {
      return scp.tOfymax();
    }
  };

  // класс для вычисления минимальной и максимальной общей координаты
  class CommonMajorCoord
  {
  protected:
    CommonMajorCoord
    (
      const double tlow1,
      const double thigh1,
      const baseCurve& crv1,
      const double tlow2,
      const double thigh2,
      const baseCurve& crv2
    );

    const
      std::tuple<double, double, geom2d::point, double, geom2d::point>
        getMinCoords(const CurveMajorHelper& getter1, const CurveMajorHelper& getter2);
    const
      std::tuple<double, double, geom2d::point, double, geom2d::point>
        getMaxCoords(const CurveMajorHelper& getter1, const CurveMajorHelper& getter2);

  protected:

    const double tmin1;
    const double tmax1;
    const baseCurve& curve1;

    const double tmin2;
    const double tmax2;
    const baseCurve& curve2;
  };

  class XCommonMajorCoord : private CommonMajorCoord
  {
  public:
    XCommonMajorCoord
    (
      const double tlow1,
      const double thigh1,
      const baseCurve& crv1,
      const double tlow2,
      const double thigh2,
      const baseCurve& crv2
    ) :CommonMajorCoord(tlow1, thigh1, crv1, tlow2, thigh2, crv2)
    {
      setXmin();
      setXmax();
    }

  private:
    
    double xmin = 0.0;
    double tofxmin1 = 0.0;
    point pntofxmin1;
    double tofxmin2 = 0.0;
    point pntofxmin2;

    double xmax = 0.0;
    double tofxmax1 = 0.0;
    point pntofxmax1;
    double tofxmax2 = 0.0;
    point pntofxmax2;

    void setXmin();
    void setXmax();
  };

  class YCommonMajorCoord : private CommonMajorCoord
  {
  public:
    YCommonMajorCoord
    (
      const double tlow1,
      const double thigh1,
      const baseCurve& crv1,
      const double tlow2,
      const double thigh2,
      const baseCurve& crv2
    ) :CommonMajorCoord(tlow1, thigh1, crv1, tlow2, thigh2, crv2)
    {
      setYmin();
      setYmax();
    }

  private:

    double ymin = 0.0;
    double tofymin1 = 0.0;
    point pntofymin1;
    double tofymin2 = 0.0;
    point pntofymin2;

    double ymax = 0.0;
    double tofymax1 = 0.0;
    point pntofymax1;
    double tofymax2 = 0.0;
    point pntofymax2;

    void setYmin();
    void setYmax();
  };

}
