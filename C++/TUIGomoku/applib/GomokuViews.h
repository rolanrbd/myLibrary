#pragma once
#ifndef GomokuView_H
#define GomokuView_H
#include "Gomoku.h"
#include "GomokuCommand.h"

using namespace std;

namespace Gomoku
{
	

	class Command;

	//my observer
	class View {
	public:
		virtual void update(Command* c) = 0;
		virtual ViewType getType() = 0;
		virtual ~View();;
	};

	class BoardView : public View {
	public:
		BoardView();
		ViewType getType() override;
		void update(Command* c) override;
	};

	class ScoreView : public View
	{
		int x, y;
		int longestRun;
	public:
		ScoreView(int _x, int _y);
		ViewType getType() override;
		void update(Command* c) override;
	};

	class MovesView : public View {
		int x, y;
	public:
		MovesView(int _x, int _y);
		ViewType getType() override;
		void update(Command* c) override;

	};

}
#endif