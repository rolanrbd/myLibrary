#pragma once

#ifndef Gomoku_H
#define Gomoku_H

#include "WinConsole.h"
#include "IAdapterConsole.h"

namespace Gomoku
{
	WORD const FOREGROUND_BLACK = 0;
	WORD const BACKGROUND_WHITE = BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_BLUE;

	enum class State {
		NO_PLAYER,
		BLUE_PLAYER,
		RED_PLAYER
	};

	enum class ViewType {
		BOARD_VIEW,
		MOVES_VIEW,
		SCORE_VIEW
	};
	
}
#endif