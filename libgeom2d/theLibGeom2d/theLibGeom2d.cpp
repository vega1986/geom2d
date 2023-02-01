// theLibGeom2d.cpp : Определяет экспортируемые функции для DLL.
//

#include "pch.h"
#include "framework.h"
#include "theLibGeom2d.h"


// Пример экспортированной переменной
THELIBGEOM2D_API int ntheLibGeom2d=0;

// Пример экспортированной функции.
THELIBGEOM2D_API int fntheLibGeom2d(void)
{
    return 0;
}

// Конструктор для экспортированного класса.
CtheLibGeom2d::CtheLibGeom2d()
{
    return;
}
