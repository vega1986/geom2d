#include "curveAnalizerBase.h"
#include "findFunctionRoots.h"

#include <vector>

// ********************************************************************************************************************

void geom2d::curveAnalizerBase::perform
(
  const double tmin1,
  const double tmax1,
  const double tmin2,
  const double tmax2
)
{
  //const auto pBegOfCurve1 = m_curve1.getPoint(tmin1);
  //const auto pEndOfCurve1 = m_curve1.getPoint(tmax1);

  //const auto pBegOfCurve2 = m_curve2.getPoint(tmin2);
  //const auto pEndOfCurve2 = m_curve2.getPoint(tmax2);

  //// возвращаем класс кривой
  //auto getCurveClass = [](const point pBegOfCurve, const point pEndOfCurve) -> curveClass
  //{
  //  const auto [xbeg, ybeg] = pBegOfCurve;
  //  const auto [xend, yend] = pEndOfCurve;
  //  constexpr auto tol2d = math::tolerance::tolPoint;
  //  //constexpr auto halfTol2d = tol2d / 2.0;

  //  const bool xInc = (xend - xbeg) > tol2d;
  //  const bool xDec = (xbeg - xend) > tol2d;

  //  const bool yInc = (yend - ybeg) > tol2d;
  //  const bool yDec = (ybeg - yend) > tol2d;

  //  if (xInc)
  //  {
  //    if (yDec)
  //    {
  //      return curveClass::Screen;
  //    }
  //    else if (yInc)
  //    {
  //      return curveClass::Normal;
  //    }
  //    else
  //    {
  //      return curveClass::PlatoY;
  //    }
  //  }
  //  else if (xDec)
  //  {
  //    if (yDec)
  //    {
  //      return curveClass::Normal;
  //    }
  //    else if (yInc)
  //    {
  //      return curveClass::Screen;
  //    }
  //    else
  //    {
  //      return curveClass::PlatoY;
  //    }
  //  }
  //  else
  //  {
  //    if (yDec)
  //    {
  //      return curveClass::PlatoX;
  //    }
  //    else if (yInc)
  //    {
  //      return curveClass::PlatoX;
  //    }
  //    else
  //    {
  //      return curveClass::Point;
  //    }
  //  }
  //};

  //const auto classOfCurve1 = getCurveClass(pBegOfCurve1, pEndOfCurve1);
  //const auto classOfCurve2 = getCurveClass(pBegOfCurve2, pEndOfCurve2);

  const auto classOfCurve1 = baseCurve::getCurveClass(tmin1, tmax1, m_curve1);
  const auto classOfCurve2 = baseCurve::getCurveClass(tmin2, tmax2, m_curve2);

  // решаем задачу для каждого возможного сочетания классов отдельно
  switch (classOfCurve1)
  {
  case curveClass::Screen:
    switch (classOfCurve2)
    {
    case curveClass::Screen:
    {
      performScreen1Screen2(tmin1, tmax1, tmin2, tmax2);
#if 0
      
#endif
    }
    break;
    case curveClass::Normal:
    {
      performScreen1Normal2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::PlatoX:
    {
      performScreen1PlatoX2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::PlatoY:
    {
      performScreen1PlatoY2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::Point:
    {
      performScreen1Point2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    }
    break;
  case curveClass::Normal:
    switch (classOfCurve2)
    {
    case curveClass::Screen:
    {
      performNormal1Screen2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::Normal:
    {
      performNormal1Normal2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::PlatoX:
    {
      performNormal1PlatoX2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::PlatoY:
    {
      performNormal1PlatoY2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::Point:
    {
      performNormal1Point2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    }
    break;
  case curveClass::PlatoX:
    switch (classOfCurve2)
    {
    case curveClass::Screen:
    case curveClass::Normal:
    {
      performPlatoX1Any2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::PlatoX:
    {
      performPlatoX1PlatoX2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::PlatoY:
    {
      performPlatoX1PlatoY2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::Point:
    {
      performPlatoX1Point2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    }
    break;
  case curveClass::PlatoY:
    switch (classOfCurve2)
    {
    case curveClass::Screen:
    case curveClass::Normal:
    {
      performPlatoY1Any2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::PlatoX:
    {
      performPlatoY1PlatoX2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::PlatoY:
    {
      performPlatoY1PlatoY2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::Point:
    {
      performPlatoY1Point2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    }
    break;
  case curveClass::Point:
    switch (classOfCurve2)
    {
    case curveClass::Screen:
    case curveClass::Normal:
    {
      performPoint1Any2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::PlatoX:
    {
      performPoint1PlatoX2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::PlatoY:
    {
      performPoint1PlatoY2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    case curveClass::Point:
    {
      performPoint1Point2(tmin1, tmax1, tmin2, tmax2);
    }
    break;
    }
    break;
  }
}

// ********************************************************************************************************************

void geom2d::curveAnalizerBase::postProcessing()
{
}

// ********************************************************************************************************************

std::set<double> geom2d::curveAnalizerBase::rootsOfCurveVelocity(const baseCurve& curve)
{
  const double tmin1 = curve.parameterMin();
  const double tmax1 = curve.parameterMax();

  // используем свой компаратор для исключения повторяющихся значений параметра
  auto curveParameterLess = [](const double lhs, const double rhs) -> bool
  {
    return (rhs - lhs) > math::tolerance::tolNumeric;
  };

  using parameterContainer = std::set<double, decltype(curveParameterLess)>;

  // I - ищем множество значений параметра, в которых кривая по x меняет знак монотонности
  // то есть, либо меняет уменьшение на увеличение, либо наоборот
  auto curve1u = [&curve](const double t) -> double
  {
    return curve.getVelocity(t).x;
  };
  const auto rootsOfCurve1u = math::findFunctionRoots(tmin1, tmax1, curve1u);

  // II - ищем множество значений параметра, в которых кривая по y меняет знак монотонности
  // то есть, либо меняет уменьшение на увеличение, либо наоборот
  auto curve1v = [&curve](const double t) -> double
  {
    return curve.getVelocity(t).y;
  };
  const auto rootsOfCurve1v = math::findFunctionRoots(tmin1, tmax1, curve1v);

  parameterContainer allRoots(curveParameterLess);
  // собираем корни по u
  for (const auto t : rootsOfCurve1u)
  {
    allRoots.insert(t);
  }
  // собираем корни по v
  for (const auto t : rootsOfCurve1v)
  {
    allRoots.insert(t);
  }
  // таким образом исключили все "дупликаты"
  // теперь собираем всё в обычный std::set
  std::set<double> result;
  for (const auto t : allRoots)
  {
    result.insert(t);
  }
  return result;
}

// ********************************************************************************************************************

void geom2d::curveAnalizerBase::perform()
{
  // Разбиваем кривую 1 на точки монотонности
  const auto tmanifold1 = rootsOfCurveVelocity(m_curve1);

  // Тоже для кривой 2
  const auto tmanifold2 = rootsOfCurveVelocity(m_curve2);

  // точки монотонности не включают начало и конец
  // для каждой пары значений интервала монотонности и каждой кривой ищем пересечение,
  // зная что кривые монотонны на выбранном участке
  // tmanifold - отсортированный массив t
  auto assembleVectorOft = [](const double tmin, const double tmax, const auto tmanifold)
  {
    std::vector<double> result;
    result.push_back(tmin);
    for (const auto t : tmanifold)
    {
      result.push_back(t);
    }
    result.push_back(tmax);
    return result;
  };
  const auto tman1 = assembleVectorOft(m_curve1.parameterMin(), m_curve1.parameterMax(), tmanifold1);
  const auto tman2 = assembleVectorOft(m_curve2.parameterMin(), m_curve2.parameterMax(), tmanifold2);
  for (size_t i = 0; i < tman1.size() - 1; ++i)
  {
    const double tmin1 = tman1[i];
    const double tmax1 = tman1[i + 1];
    for (size_t j = 0; j < tman2.size() - 1; ++j)
    {
      const double tmin2 = tman2[j];
      const double tmax2 = tman2[j + 1];
      perform(tmin1, tmax1, tmin2, tmax2);
    }
  }
  postProcessing();
}

// ********************************************************************************************************************

