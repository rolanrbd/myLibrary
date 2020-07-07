#pragma comment (lib, "applib-mt.lib")

#include "GomokuApp.h"

/* Created by: Claudia Gonzalez
*/
using namespace Gomoku;
class MyGomokuApp : public GomokuApp {
	int execute() override {

		GomokuController* game = new GomokuController();
		game->ExecuteInitCommand();

		return 0;
	}
} MyGomokuApp;


