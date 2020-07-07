#include "GomokuModel.h"
#include "GomokuViews.h"


Gomoku::Model::Model(AdapterConsole * a) { adapt = a; }

AdapterConsole * Gomoku::Model::getAdapter() { return adapt; }

bool Gomoku::Model::getGameState() { return gameOver; }

void Gomoku::Model::setGameState(BOOL b) { gameOver = b; }

ClickCoords * Gomoku::Model::getLastClick() { return adapt->getLastClick(); }

void Gomoku::Model::setLastClick(ClickCoords * click) { lastClick = click; notify(); }

Gomoku::Model::~Model() {}


void Gomoku::GomokuModel::StartPlaying()
{
}

Gomoku::GomokuModel::GomokuModel(AdapterConsole * a) : Model(a)
{
	//fills out the movesMatrix
	for (int i = 0; i < 15; i++)
		for (int j = 0; j < 15; j++)
		{
			movesMatrix[i][j] = State::NO_PLAYER;
		}

	//StartPlaying();
}

void Gomoku::GomokuModel::Attach(View * obsvr) { obsvrs.push_back(obsvr); }

void Gomoku::GomokuModel::Detach(View * obsvr) { obsvrs.pop_back(); }

BOOL Gomoku::GomokuModel::getPlayerInTurn() { return isRedTurn; }

void Gomoku::GomokuModel::setPlayerInTurn(bool st) { isRedTurn = st; }

Gomoku::State Gomoku::GomokuModel::getState(int i, int j) { return movesMatrix[i][j]; }

void Gomoku::GomokuModel::setState(int i, int j, State st) { movesMatrix[i][j] = st; }

bool Gomoku::GomokuModel::isWinMove(int xStart, int yStart, State playerInTurn)
{
	bool result = false;
	bool searchNext = true;
	int x = xStart, y = yStart;

	int countInLine = 1;
	//look up
	while (searchNext)
	{   // [x][y-1]
		if ((y - 1) >= 0 && movesMatrix[x][y - 1] == playerInTurn)
		{
			countInLine++;
			if (countInLine >= 5)
				return true;
			y--;
		}
		else
		{
			searchNext = false;
			y = yStart;
		}
	}

	//look down
	searchNext = true;
	while (searchNext)
	{   // g_State[x][y+1]
		if ((y + 1) <= 14 && movesMatrix[x][y + 1] == playerInTurn)
		{
			countInLine++;
			if (countInLine >= 5)
				return true;
			y++;
		}
		else
		{
			searchNext = false;
			y = yStart;
		}
	}

	//look left
	countInLine = 1;
	//look left and up
	searchNext = true;
	while (searchNext)
	{   // g_State[x-1][y]
		if ((x - 1) >= 0 && movesMatrix[x - 1][y] == playerInTurn)
		{
			countInLine++;
			if (countInLine >= 5)
				return true;
			x--;
		}
		else
		{
			searchNext = false;
			x = xStart;
		}
	}
	//look right
	searchNext = true;
	while (searchNext)
	{   // g_State[x+1][y]
		if ((x + 1) <= 14 && movesMatrix[x + 1][y] == playerInTurn)
		{
			countInLine++;
			if (countInLine >= 5)
				return true;
			x++;
		}
		else
		{
			searchNext = false;
			x = xStart;
		}
	}
	//look diagonal left up
	countInLine = 1;
	searchNext = true;
	while (searchNext)
	{   // g_State[x-1][y-1]
		if (((x - 1) >= 0 && (y - 1 >= 0)) && movesMatrix[x - 1][y - 1] == playerInTurn)
		{
			countInLine++;
			if (countInLine >= 5)
				return true;
			x--; y--;
		}
		else
		{
			searchNext = false;
			x = xStart; y = yStart;
		}
	}
	//look diagonal left down
	searchNext = true;
	while (searchNext)
	{   // g_State[x+1][y+1]
		if (((x + 1) <= 14 && (y + 1 <= 14)) && movesMatrix[x + 1][y + 1] == playerInTurn)
		{
			countInLine++;
			if (countInLine >= 5)
				return true;
			x++; y++;
		}
		else
		{
			searchNext = false;
			x = xStart; y = yStart;
		}
	}
	//look diagonal right up
	countInLine = 1;
	searchNext = true;

	while (searchNext)
	{   // g_State[x+1][y-1]
		if (((x + 1) <= 14 && (y - 1 >= 0)) && movesMatrix[x + 1][y - 1] == playerInTurn)
		{
			countInLine++;
			if (countInLine >= 5)
				return true;
			x++; y--;
		}
		else
		{
			searchNext = false;
			x = xStart; y = yStart;
		}
	}
	//look diagonal right down
	searchNext = true;
	while (searchNext)
	{   // g_State[x-1][y+1]
		if (((x - 1) >= 0 && (y + 1 <= 14)) && movesMatrix[x - 1][y + 1] == playerInTurn)
		{
			countInLine++;
			if (countInLine >= 5)
				return true;
			x--; y++;
		}
		else
		{
			searchNext = false;
			x = xStart; y = yStart;
		}
	}
	return result;
}

void Gomoku::GomokuModel::notify()
{
	for (size_t i = 0; i < obsvrs.size(); i++)
	{
		switch (obsvrs[i]->getType())
		{
		case ViewType::BOARD_VIEW:
			obsvrs[i]->update(new DrawPlayCommand(this));
			break;
		case ViewType::MOVES_VIEW:
			obsvrs[i]->update(new DrawMovesCommand(this));
			break;
		case ViewType::SCORE_VIEW:
			obsvrs[i]->update(new DrawScoreCommand(this));
			break;
		}
	}


}

Gomoku::GomokuModel::~GomokuModel() {}




