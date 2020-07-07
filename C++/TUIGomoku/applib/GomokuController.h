#pragma once
#ifndef GomokuController_H
#define GomokuController_H

#include "WinConsole.h"
#include "IAdapterConsole.h"
#include "AdapterConsole.h"

#include "SingletonGomoku.h"
#include "GomokuModel.h"
#include "GomokuViews.h"
#include "GomokuCommand.h"

/* Created by: Claudia Gonzalez
Client
*/


namespace Gomoku
{
	class View;
	class ScoreView;
	class BoardView;
	class MovesView;

	class GomokuController : public SingletonGomoku {
		Model* modelSubject;

	public:
		GomokuController();
		void ExecuteInitCommand();
		~GomokuController();
	};
}
#endif
