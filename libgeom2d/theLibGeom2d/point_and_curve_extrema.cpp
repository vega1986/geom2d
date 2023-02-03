#include "point_and_curve_extrema.h"
#include "avector.h"
#include "findFunctionRoots.h"

std::pair<double, double>
  geom2d::extrema::point_and_curve::nearestPoint::execute() const
{
  constexpr size_t n = 100;

  double dist = math::infinite::distance;
  double tofExtrema = 0.0;

  const auto& crv = curve;
  const auto pnt = P;
  auto getScalarMultiplication = [pnt, &crv](const double t) -> const double
  {
    const auto q = crv.getPoint(t);
    const auto tau = vector{ crv.getVelocity(t) };
    const vector out{ q, pnt };
    if (point::isSame(q, pnt)) return 0.0;
    return (tau.normalized(), out);
  };

  for (size_t j = 0; j < n; ++j)
  {
    const double alfa_l = static_cast<double>(j      ) / static_cast<double>(n);
    const double beta_l = 1.0 - alfa_l;

    const double alfa_r = static_cast<double>(j + 1.0) / static_cast<double>(n);
    const double beta_r = 1.0 - alfa_r;

    const auto tl = tmin * beta_l + tmax * alfa_l;
      //static_cast<double>(j    ) * (tmax - tmin) / static_cast<double>(n);
    const auto scml = getScalarMultiplication(tl);

    const auto tr = tmin * beta_r + tmax * alfa_r;
      //static_cast<double>(j + 1) * (tmax - tmin) / static_cast<double>(n);
    const auto scmr = getScalarMultiplication(tr);

    if (scml == 0.0)
    {
      // экстрема найдена
      const auto pntOnCurve = curve.getPoint(tl);
      const auto curDist = point::distance(P, pntOnCurve);

      if (curDist < dist)
      {
        dist = curDist;
        tofExtrema = tl;
      }
    }
    else if (scmr == 0.0)
    {
      // экстрема найдена
      const auto pntOnCurve = curve.getPoint(tr);
      const auto curDist = point::distance(P, pntOnCurve);

      if (curDist < dist)
      {
        dist = curDist;
        tofExtrema = tr;
      }
    }
    else
    {
      if (((scml > 0.0) and (scmr > 0.0)) or ((scml < 0.0) and (scmr < 0.0)))
      {
        continue;
      }
      else
      {
        // решаем уравнение
        const auto tofnorm = math::findUniqueFunctionRoot(tl, tr, getScalarMultiplication);
        const auto ponCurve = curve.getPoint(tofnorm);
        const auto curDist = point::distance(P, ponCurve);
        if (curDist < dist)
        {
          dist = curDist;
          tofExtrema = tofnorm;
        }
      }
    }
  }
  // тестируем кончики) как обычно
  for (const auto tonCurve : std::initializer_list{tmin, tmax})
  {
    const auto pntOfVertex = curve.getPoint(tonCurve);
    const auto curDist = point::distance(pntOfVertex, P);
    if (curDist < dist)
    {
      dist = curDist;
      tofExtrema = tonCurve;
    }
  }
  return std::pair{ dist, tofExtrema };
}
