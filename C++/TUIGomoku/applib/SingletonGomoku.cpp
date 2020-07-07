#include <SingletonGomoku.h>
#include <cassert>

Console::SingletonGomoku* Console::SingletonGomoku::thisSingleton = nullptr;
Console::AdapterConsole* Console::SingletonGomoku::window = nullptr;
Console::SingletonGomoku::SingletonGomoku()
{
	assert(!thisSingleton && "Error: Already initialized!");
	thisSingleton = this;
	window = new AdapterConsole(new WinConsole("Gomoku Gang", 89, 34));
}


Console::SingletonGomoku* Console::SingletonGomoku::getInstance()
{
	if (thisSingleton == 0)
	{
		thisSingleton = new SingletonGomoku();
	}

	return thisSingleton;
}

Console::AdapterConsole * Console::SingletonGomoku::getAdapter()
{
	return window;
}

