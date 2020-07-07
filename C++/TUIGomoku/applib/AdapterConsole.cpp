#include <AdapterConsole.h>

Console::AdapterConsole::AdapterConsole(IWinConsole* winConsole) {
	window = winConsole;
	lastClick = nullptr;
	window->ResizeWindow();
}

/*void AdapterConsole::SetConsoleSize(SHOR) {
	WINDOW_HEIGHT = wh.X;
	WINDOW_WIDTH = wh.Y;
}*/

void Console::AdapterConsole::SaveConsoleState() {
	window->SaveConsoleState();
}

void Console::AdapterConsole::ResizeWindow() {
	window->ResizeWindow();
}

void Console::AdapterConsole::HideCursor() {
	window->HideCursor();
}

void Console::AdapterConsole::PaintScreen(int w, int h, WORD color) {
	window->Draw(w, h, color);
}
void Console::AdapterConsole::DrawMovesList(int w, int h, WORD color, string msg) {
	window->DrawMovesList(w, h, color, msg);
}
void Console::AdapterConsole::DrawScores(int xCoord, int yCoord, string msg) {
	window->DrawScores(xCoord, yCoord, msg);
}

void Console::AdapterConsole::DrawMatrix() {
	window->DrawMatrix();
}

void Console::AdapterConsole::DrawWin(string winnerColor) {
	window->DrawWin(winnerColor);
}

void Console::AdapterConsole::RestoreWindow() {
	window->RestoreWindow();
}


void Console::AdapterConsole::ClearWindow() {
	window->ClearWindow();
}

INPUT_RECORD Console::AdapterConsole::GetInput() {
	return window->ReadInputConsole();
}

void Console::AdapterConsole::Initialize() {
	//window->ResizeWindow();
	INPUT_RECORD inputRecord;
	while (running) {
		inputRecord = GetInput();
		switch (inputRecord.EventType) {
		case WINDOW_BUFFER_SIZE_EVENT:
			break;
		case KEY_EVENT:
			ProcessKeyEvent(inputRecord.Event.KeyEvent);
			break;
		case MOUSE_EVENT:
			if(!MouseEventProc(inputRecord.Event.MouseEvent))
				return;
			break;
		}
	}
}

BOOL  Console::AdapterConsole::CtrlHandler(DWORD ctrlType) {
	switch (ctrlType) {
	case CTRL_C_EVENT:
		cout << "Ctrl-C pressed\n" << endl;
		running = false;
		return TRUE;
	}
	return FALSE;
}

void  Console::AdapterConsole::ProcessKeyEvent(KEY_EVENT_RECORD const& ker) {
	if (ker.uChar.AsciiChar == 'r')
		RestoreWindow();
}

//WORK AROUND for building the matrix
int clickCount = 0;

bool Console::AdapterConsole::MouseEventProc(MOUSE_EVENT_RECORD const& mer) {
	auto mask = mer.dwButtonState;
	//left clicks alternates between colors
	if (mask&FROM_LEFT_1ST_BUTTON_PRESSED)
	{
		clickCount++;
		if (clickCount == 1)
		{
			DrawMatrix();			
		}
		else
		{
			//save the console state before click in case the user wants to undo the move
			SaveConsoleState();

			lastClick = new ClickCoords();
			lastClick->x = mer.dwMousePosition.X;
			lastClick->y = mer.dwMousePosition.Y;
			
			return false;
		}
		
	}
	//undo last move
	else if (mask&RIGHTMOST_BUTTON_PRESSED) {
		RestoreWindow();
		return false;
	}
	return true;
}

int Console::AdapterConsole::getCellHeight()
{
	//**static for now
	return 2;
}

int Console::AdapterConsole::getCellWidth()
{
	//**static for now
	return 3; //6
}
