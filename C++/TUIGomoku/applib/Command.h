#pragma once
/*#include <GomokuModel.h>


class Command
{
public:
	Command() {}
	virtual ~Command() {}
	virtual void Execute() = 0;
};

class InitCommand : public Command
{
protected:
	Model * modelSubject;
public:
	InitCommand(Model *mdlSbjct) : Command() { modelSubject = mdlSbjct; }
	void Execute()
	{
		ClickCoords* c;
		while (modelSubject->getGameState())
		{
			c = modelSubject->getAdapter->getLastClick();
			if (c != nullptr)
			{
				modelSubject->setLastClick(c);
			}
		}
	}
};

/*class DrawPlayCommand : public Command
{
protected:
	GomokuModel * mdlSbjct;
public:
	DrawPlayCommand(GomokuModel *_mdlSbjct) : Command() { mdlSbjct = _mdlSbjct; }
	void Execute()
	{
		ClickCoords* c = mdlSbjct->getLastClick();

		// Determine which tile was clicked
		int i = c->x / mdlSbjct->getAdapter()->getCellWidth(); //3 is the width of each cell
		int j = c->y / mdlSbjct->getAdapter()->getCellHeight(); //2 is the height of each cell

		if (c->y >= 4 && c->x <= 44)
		{
			if (mdlSbjct->getState(i, j) == State::NO_PLAYER)
			{
				//check who's turn is
				if (mdlSbjct->getPlayerInTurn())
				{
					mdlSbjct->setState(i, j, State::RED_PLAYER);
					mdlSbjct->getAdapter()->PaintScreen(i, j, BACKGROUND_RED); //c->x, c->y
					mdlSbjct->getAdapter()->setLastClick(nullptr);

					if (mdlSbjct->isWinMove(i, j, State::RED_PLAYER))
						mdlSbjct->getAdapter()->DrawWin("red");
					else
						mdlSbjct->setPlayerInTurn(false);
				}
				else
				{
					mdlSbjct->setState(i, j, State::BLUE_PLAYER);
					//paint cell with color blue
					mdlSbjct->getAdapter()->PaintScreen(i, j, BACKGROUND_BLUE);//c->x, c->y
					mdlSbjct->getAdapter()->setLastClick(nullptr);

					if (mdlSbjct->isWinMove(i, j, State::BLUE_PLAYER))
						mdlSbjct->getAdapter()->DrawWin("blue");
					else
						mdlSbjct->setPlayerInTurn(true);
				}
			}
		}
		mdlSbjct->getAdapter()->Initialize();
	}
};

class DrawScoreCommand : public Command
{
protected:
	GomokuModel * modelSubject;
public:
	DrawScoreCommand(GomokuModel *mdlSbjct) : Command() { modelSubject = mdlSbjct; }
	void Execute() {
	
		modelSubject()->DrawScores(2, 1, buildScoresMsg("RED", redCounter++));
	
	}
};

class DrawMovesCommand : public Command
{
	GomokuModel* modelSubject;
public:
	DrawMovesCommand(GomokuModel *mdlSbjct) : Command() { modelSubject = mdlSbjct; }
	void Execute() {}
};*/