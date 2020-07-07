#include "GomokuViews.h"

Gomoku::BoardView::BoardView() {}

Gomoku::ViewType Gomoku::BoardView::getType() { return ViewType::BOARD_VIEW; }

void Gomoku::BoardView::update(Command * c) { c->Execute(); }

Gomoku::ScoreView::ScoreView(int _x, int _y) { x = _x; y = _y; longestRun = 0; }

Gomoku::ViewType Gomoku::ScoreView::getType() { return ViewType::SCORE_VIEW; }

void Gomoku::ScoreView::update(Command * c)
{/*
 //pintar el score view
 auto m = game.getMove();
 int lRun = game->getLongestRunOfLastMove();
 if (lRun > longestRun &&  m.getPlayerColor() != colorLastLongestRun)
 {
 longestRun = lRun;
 colorLastLongestRun = m.getPlayerColor();
 string msg = "color tal con tal cantidad";
 s->getAdapter()->DrawScores(x, y, msg)
 }
 else if (lRun > longestRun &&  m.getPlayerColor() == colorLastLongestRun)
 {
 longestRun = lRun;
 string msg = "color tal con tal cantidad";
 s->getAdapter()->DrawScores(x, y, msg)
 }*/
}

Gomoku::MovesView::MovesView(int _x, int _y) { x = _x; y = _y; }

Gomoku::ViewType Gomoku::MovesView::getType() { return ViewType::MOVES_VIEW; }

void Gomoku::MovesView::update(Command * c)
{
	/*
	//pintar el moves list
	auto m = game.getMove();
	string msg = "";
	WORD color;
	if (m.getPlayerColor() == "red")
	{
	msg = "red";
	color = RED;
	}
	else
	{
	msg = "blue";
	color = BLUE;
	}
	y++;
	s->getAdapter()->DrawMovesList(x, y, color, msg);*/
}

Gomoku::View::~View() {}
