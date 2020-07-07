#pragma once
#ifndef GomokuCommand_H
#define GomokuCommand_H

#include "GomokuModel.h"


namespace Gomoku
{
	class Model;
	class GomokuModel;

	class Command
	{
	public:
		virtual ~Command();
		virtual void Execute() = 0;
	};

	class InitCommand : public Command
	{
	protected:
		Model * modelSubject;
	public:
		InitCommand(Model *mdlSbjct);
		void Execute();
	};

	class DrawPlayCommand : public Command
	{
	protected:
		GomokuModel * mdlSbjct;
	public:
		DrawPlayCommand(GomokuModel *_mdlSbjct);
		void Execute();
	};

	class DrawScoreCommand : public Command
	{
	protected:
		GomokuModel * modelSubject;
	public:
		DrawScoreCommand(GomokuModel *mdlSbjct);
		void Execute();
	};

	class DrawMovesCommand : public Command
	{
		GomokuModel* modelSubject;
	public:
		DrawMovesCommand(GomokuModel *mdlSbjct);
		void Execute();
	};
}

#endif