#include <GomokuApp.h>

GomokuApp* GomokuApp::thisApp = nullptr;

GomokuApp::GomokuApp()
{
	thisApp = this;
}

int GomokuApp::main()
{
	return thisApp->execute();
}

int main() {
	return GomokuApp::main();
}
