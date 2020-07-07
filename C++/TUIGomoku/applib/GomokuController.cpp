

#include "GomokuController.h"


Gomoku::GomokuController::GomokuController() :SingletonGomoku()
{

	getAdapter()->Initialize();

	modelSubject = new GomokuModel(getAdapter());

	modelSubject->Attach(new ScoreView(1, 1));
	modelSubject->Attach(new BoardView());
	modelSubject->Attach(new MovesView(1, 1));

}

void Gomoku::GomokuController::ExecuteInitCommand()
{
	InitCommand *c = new InitCommand(modelSubject);
	c->Execute();
}

Gomoku::GomokuController::~GomokuController() {}
