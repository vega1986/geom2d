// Приведенный ниже блок ifdef — это стандартный метод создания макросов, упрощающий процедуру
// экспорта из библиотек DLL. Все файлы данной DLL скомпилированы с использованием символа THELIBGEOM2D_EXPORTS
// Символ, определенный в командной строке. Этот символ не должен быть определен в каком-либо проекте,
// использующем данную DLL. Благодаря этому любой другой проект, исходные файлы которого включают данный файл, видит
// функции THELIBGEOM2D_API как импортированные из DLL, тогда как данная DLL видит символы,
// определяемые данным макросом, как экспортированные.
#ifdef THELIBGEOM2D_EXPORTS
#define THELIBGEOM2D_API __declspec(dllexport)
#else
#define THELIBGEOM2D_API __declspec(dllimport)
#endif

#if 0
// Этот класс экспортирован из библиотеки DLL
class THELIBGEOM2D_API CtheLibGeom2d {
public:
	CtheLibGeom2d(void);
	// TODO: добавьте сюда свои методы.
};

extern THELIBGEOM2D_API int ntheLibGeom2d;

THELIBGEOM2D_API int fntheLibGeom2d(void);
#endif