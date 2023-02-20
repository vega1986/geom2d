#pragma once
#include "baseCurve.h"
#include "IntersecctionSolutionType.h"
#include "StatOfCurvePiece.h"
#include "CommonRangeHelper.h"

namespace geom2d
{
  namespace uniqueIntersection
  {
    // Класс для поиска единственной точки пересечения между кривыми на участках монотонности.
    // В основе методики лежит метод деления ОДЗ other кривой пополам. ОДЗ кривой other обновляется, если на его концах
    // разность ординаты кривой other и reference имеет разный знак (что свидетельствует о пересечении).
    // В качестве оси абсцисс может использоваться как ось X так и ось Y - поэтому класс параметризован специальным классом геттером.
    template <typename DataGetter>
    class solver
    {
    public:
      solver(const double trmin, const double trmax, const baseCurve& reference, const double tomin, const double tomax, const baseCurve& other)
        :
        trefmin(trmin), trefmax(trmax), referenceCurve(reference),
        tothmin(tomin), tothmax(tomax), otherCurve(other)
      {}

      std::optional<IntersecctionSolutionType> execute() const
      {
        DataGetter referenceGetter{ trefmin , trefmax , referenceCurve };
        DataGetter otherGetter{ tothmin , tothmax , otherCurve };

        if (otherGetter.getMaxCoord() <= referenceGetter.getMinCoord())
        {
          const point point_oth = otherGetter.getPointOfMax();
          const double t_oth = otherGetter.getTofMax();

          const point point_ref = referenceGetter.getPointOfMin();
          const double t_ref = referenceGetter.getTofMin();

          if (point::isSame(point_oth, point_ref))
          {
            return IntersecctionSolutionType{ 0.5 * (point_oth + point_ref), t_ref, t_oth };
          }
          else
          {
            return std::nullopt;
          }
        }
        else if (otherGetter.getMinCoord() >= referenceGetter.getMaxCoord())
        {
          const point point_oth = otherGetter.getPointOfMin();
          const double t_oth = otherGetter.getTofMin();

          const point point_ref = referenceGetter.getPointOfMax();
          const double t_ref = referenceGetter.getTofMax();

          if (point::isSame(point_oth, point_ref))
          {
            return IntersecctionSolutionType{ 0.5 * (point_oth + point_ref), t_ref, t_oth };
          }
          else
          {
            return std::nullopt;
          }
        }
        else
        {
          ///////////////
          //           //
          // MIN COORD //
          //           //
          ///////////////
          auto [
            commonMin,
            trefOfCommonMin,
            pointRefOfCommonMin,
            tothOfCommonMin,
            pointOthOfCommonMin
          ] = CommonRangeHelper::ofLowest(referenceGetter, otherGetter);

          ///////////////
          //           //
          // MAX COORD //
          //           //
          ///////////////
          auto [
            commonMax,
            trefOfCommonMax,
            pointRefOfCommonMax,
            tothOfCommonMax,
            pointOthOfCommonMax
          ] = CommonRangeHelper::ofHighest(referenceGetter, otherGetter);
          
          // проверяем концы ОДЗ по X
          //  (*) проверяем левый конец
          if (point::isSame(pointOthOfCommonMin, pointRefOfCommonMin))
          {
            return IntersecctionSolutionType{ 0.5 * (pointOthOfCommonMin + pointRefOfCommonMin), trefOfCommonMin, tothOfCommonMin };
          }
          //  (*) проверяем правый конец
          if (point::isSame(pointOthOfCommonMax, pointRefOfCommonMax))
          {
            return IntersecctionSolutionType{ 0.5 * (pointOthOfCommonMax + pointRefOfCommonMax), trefOfCommonMax, tothOfCommonMax };
          }
          // лежит ли одна кривая полностью выше или ниже другой
          double fleft = DataGetter::ordinateOf(pointOthOfCommonMin) - DataGetter::ordinateOf(pointRefOfCommonMin);
          double fright = DataGetter::ordinateOf(pointOthOfCommonMax) - DataGetter::ordinateOf(pointRefOfCommonMax);
          const bool othAboveRef = (fleft > 0.0) and (fright > 0.0);
          const bool othUnderRef = (fleft < 0.0) and (fright < 0.0);
          if (othAboveRef or othUnderRef)
          {
            return std::nullopt;
          }
          while (true)
          {
            const double tothMiddle = 0.5 * (tothOfCommonMin + tothOfCommonMax);
            const point pointOthOnMiddle = otherCurve.getPoint(tothMiddle);
            
            double trefMiddle = 0.0;
            point pointRefOnMiddle;
            // обезопасим себя на случай абсолютно вертикальной кривой otherCurve вдоль нерефересной оси
            //if (DataGetter::abscissaOf(pointOthOnMiddle) <= DataGetter::abscissaOf(pointRefOfCommonMin))
            //{
            //  trefMiddle = trefOfCommonMin;
            //  pointRefOnMiddle = pointRefOfCommonMin;
            //}
            //else if (DataGetter::abscissaOf(pointOthOnMiddle) >= DataGetter::abscissaOf(pointRefOfCommonMax))
            //{
            //  trefMiddle = trefOfCommonMax;
            //  pointRefOnMiddle = pointRefOfCommonMax;
            //}
            //else
            //{
              auto func = [theRefCurve = &referenceCurve, abscissaToFind = DataGetter::abscissaOf(pointOthOnMiddle)](const double t) -> double
              {
                return DataGetter::abscissaOf(theRefCurve->getPoint(t)) - abscissaToFind;
              };
              const double trefMinCurrent = std::min(trefOfCommonMin, trefOfCommonMax);
              const double trefMaxCurrent = std::max(trefOfCommonMin, trefOfCommonMax);
              trefMiddle = math::findUniqueFunctionRoot(trefMinCurrent, trefMaxCurrent, func);
              pointRefOnMiddle = referenceCurve.getPoint(trefMiddle);
            //}
            // проверяем середину
            if (point::isSame(pointOthOnMiddle, pointRefOnMiddle) or std::abs(tothOfCommonMax - tothOfCommonMin) < math::tolerance::tolParameter)
            {
              return IntersecctionSolutionType{ 0.5 * (pointOthOnMiddle + pointRefOnMiddle), trefMiddle, tothMiddle };
            }
            // середина не является решением, тогда решение либо слева, либо справа
            fleft = DataGetter::ordinateOf(pointOthOfCommonMin) - DataGetter::ordinateOf(pointRefOfCommonMin);
            fright = DataGetter::ordinateOf(pointOthOfCommonMax) - DataGetter::ordinateOf(pointRefOfCommonMax);
            const double fmiddle = DataGetter::ordinateOf(pointOthOnMiddle) - DataGetter::ordinateOf(pointRefOnMiddle);
            if (((fleft > 0.0) and (fmiddle < 0.0)) or ((fleft < 0.0) and (fmiddle > 0.0)))
            {
              // середина становится правым концом, так как пересечение находится слева
              pointOthOfCommonMax = pointOthOnMiddle;
              tothOfCommonMax = tothMiddle;
              pointRefOfCommonMax = pointRefOnMiddle;
              trefOfCommonMax = trefMiddle;
            }
            else if (((fright > 0.0) and (fmiddle < 0.0)) or ((fright < 0.0) and (fmiddle > 0.0)))
            {
              // середина становится левым концом, так как пересечение находится справа
              pointOthOfCommonMin = pointOthOnMiddle;
              tothOfCommonMin = tothMiddle;
              pointRefOfCommonMin = pointRefOnMiddle;
              trefOfCommonMin = trefMiddle;
            }
            else
            {
              throw std::logic_error("impossible situation in functin geom2d::uniqueIntersection::solver::execute");
            }
          }
        }
        return std::nullopt;
      }

    private:

      const double trefmin;
      const double trefmax;
      const baseCurve& referenceCurve;

      const double tothmin;
      const double tothmax;
      const baseCurve& otherCurve;
    };


  }
}