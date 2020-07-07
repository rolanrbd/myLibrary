#pragma once
#ifndef GomokuApp_H
#define GomokuApp_H
#if defined(_DEBUG) && defined(_DLL)
#pragma comment (lib, "applib-mt-gd.lib")
#elif defined(_DEBUG) && !defined(_DLL)
#pragma comment (lib, "applib-mt-sgd.lib")
#elif !defined(_DEBUG) && defined(_DLL)
#pragma comment (lib, "applib-mt.lib")
#elif !defined(_DEBUG) && !defined(_DLL)
#pragma comment (lib, "applib-mt-s.lib")
#endif

namespace Gomoku
{
	class GomokuApp
	{
		static GomokuApp* thisApp;
		friend int main();
		static int main();
	protected:
		GomokuApp();
		virtual ~GomokuApp() = default;
		virtual int execute() = 0;
	};
}
#endif