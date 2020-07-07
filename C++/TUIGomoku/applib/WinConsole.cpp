#include <WinConsole.h>


Console::WinConsole::WinConsole()
{
	SetConsoleMode(hConsoleInput, consoleInitialMode);
	SetConsoleTitleA("Gomoku Claudia Gonzalez");
}

Console::WinConsole::WinConsole(string consoleTitle, unsigned int width, unsigned int height)
{
	SetConsoleMode(hConsoleInput, consoleInitialMode);
	SetConsoleTitleA(consoleTitle.c_str());
	WINDOW_WIDTH = width;
	WINDOW_HEIGHT = height;

	//DrawMatrix();
}

Console::WinConsole::~WinConsole()
{
	FlushConsoleInputBuffer(hConsoleInput);
}

INPUT_RECORD Console::WinConsole::ReadInputConsole()
{
	vector<INPUT_RECORD> inBuffer(128);
	DWORD numEvents;

	if (!ReadConsoleInput(hConsoleInput, inBuffer.data(), (DWORD)inBuffer.size(), &numEvents)) {
		cerr << "Failed to console input\n";
	}
	
	for (size_t iEvent = 0; iEvent < numEvents; ++iEvent) {
		return inBuffer[iEvent];
	}
}

HANDLE Console::WinConsole::GetConsoleInput() {
	return hConsoleInput;
}

HANDLE Console::WinConsole::GetConsoleOutput() {
	return hConsoleOutput;
}

COORD Console::WinConsole::GetConsoleSize() {
	return COORD{ (short)WINDOW_WIDTH, (short)WINDOW_HEIGHT };
}

void Console::WinConsole::SaveConsoleState()
{
	// Get the old window/buffer size
	GetConsoleScreenBufferInfo(hConsoleOutput, &originalCSBI);

	// Save the desktop
	originalBuffer.resize(originalCSBI.dwSize.X*originalCSBI.dwSize.Y);
	originalBufferCoord = COORD{ 0 };
	bufferRect = { 0 };
	bufferRect.Right = originalCSBI.dwSize.X - 1;
	bufferRect.Bottom = originalCSBI.dwSize.Y - 1;
	ReadConsoleOutputA(hConsoleOutput, originalBuffer.data(), originalCSBI.dwSize, originalBufferCoord, &bufferRect);

	// Save the cursor
	GetConsoleCursorInfo(hConsoleOutput, &originalCCI);

	HideCursor();
}

void Console::WinConsole::ResizeWindow() {
	SMALL_RECT sr{ 0 };
	SetConsoleWindowInfo(hConsoleOutput, TRUE, &sr);

	SetConsoleScreenBufferSize(hConsoleOutput, GetConsoleSize());

	CONSOLE_SCREEN_BUFFER_INFO sbi;
	GetConsoleScreenBufferInfo(hConsoleOutput, &sbi);

	sr.Top = sr.Left = 0;
	WINDOW_WIDTH = std::min((SHORT)WINDOW_WIDTH, sbi.dwMaximumWindowSize.X);
	WINDOW_HEIGHT = std::min((SHORT)WINDOW_HEIGHT, sbi.dwMaximumWindowSize.Y);

	sr.Right = WINDOW_WIDTH - 1;
	sr.Bottom = WINDOW_HEIGHT - 1;

	SetConsoleWindowInfo(hConsoleOutput, TRUE, &sr);
	currentConsoleWidth = sr.Right - sr.Left + 1;

	HideCursor();
}

void Console::WinConsole::HideCursor() {
	auto newCCI = originalCCI;
	newCCI.bVisible = FALSE;
	SetConsoleCursorInfo(hConsoleOutput, &newCCI);
}

void Console::WinConsole::RestoreWindow() {
	//Restore the original settings/size
	SMALL_RECT sr{ 0 };
	SetConsoleWindowInfo(hConsoleOutput, TRUE, &sr);
	SetConsoleScreenBufferSize(hConsoleOutput, originalCSBI.dwSize);
	SetConsoleWindowInfo(hConsoleOutput, TRUE, &originalCSBI.srWindow);

	//restore the desktop contents
	bufferRect = SMALL_RECT{ 0 };
	bufferRect.Right = originalCSBI.dwSize.X - 1;
	bufferRect.Bottom = originalCSBI.dwSize.Y - 1;

	WriteConsoleOutputA(hConsoleOutput, originalBuffer.data(), originalCSBI.dwSize, originalBufferCoord, &bufferRect);
	SetConsoleTitleA(originalTitle.data());

	//Restore the cursor
	SetConsoleCursorInfo(hConsoleOutput, &originalCCI);
	SetConsoleCursorPosition(hConsoleOutput, originalCSBI.dwCursorPosition);

	HideCursor();
}

void Console::WinConsole::ClearWindow() {
	GetConsoleScreenBufferInfo(hConsoleOutput, &csbi);
	dwConSize = csbi.dwSize.X * csbi.dwSize.Y;
	FillConsoleOutputCharacter(hConsoleOutput, TEXT(' '), dwConSize, screenCoords, &cCharsWritten);
	GetConsoleScreenBufferInfo(hConsoleOutput, &csbi);
	FillConsoleOutputAttribute(hConsoleOutput, csbi.wAttributes, dwConSize, screenCoords, &cCharsWritten);
	SetConsoleCursorPosition(hConsoleOutput, screenCoords);
}

void Console::WinConsole::DrawMatrix() {
	DrawTopBar();
	DrawRightSide();
	SHORT yCoord = 4;
	int xCount = 0;
	SHORT xCoord = 0;
	bool blankCharacter = true;
	for (int k = 0; k < 15; k++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int i = 0; i < 15; i++)
			{
				// Fill the title section, at the top
				dwConSize = 2;
				screenCoords = COORD{ xCoord, yCoord };
				CHAR character = blankCharacter ? ' ' : '_';
				FillConsoleOutputCharacterA(hConsoleOutput, character, dwConSize, screenCoords, &cCharsWritten);
				FillConsoleOutputAttribute(hConsoleOutput, MATRIX_ATTR, dwConSize, screenCoords, &cCharsWritten);
				xCount = xCount + 2;
				xCoord = xCount;
				dwConSize = 1;
				screenCoords = COORD{ xCoord, yCoord };
				FillConsoleOutputCharacterA(hConsoleOutput, '|', dwConSize, screenCoords, &cCharsWritten);
				FillConsoleOutputAttribute(hConsoleOutput, MATRIX_ATTR, dwConSize, screenCoords, &cCharsWritten);
				xCount = xCount + 1;
				xCoord = xCount;
			}
			yCoord = yCoord + 1;
			xCount = 0;
			xCoord = 0;
			blankCharacter = k == 14 ? true : false;
		}
		blankCharacter = true;
	}
}

void Console::WinConsole::DrawWin(string winnerColor) {

	// Get the number of character cells in the current buffer
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	GetConsoleScreenBufferInfo(hConsoleOutput, &csbi);

	// Fill the entire screen area green
	DWORD consoleSize = csbi.dwSize.X * csbi.dwSize.Y;
	COORD cursorHomeCoord{ 0, 0 };
	FillConsoleOutputCharacterA(hConsoleOutput, ' ', consoleSize, cursorHomeCoord, &cCharsWritten);
	FillConsoleOutputAttribute(hConsoleOutput, BACKGROUND_GREEN, consoleSize, cursorHomeCoord, &cCharsWritten);
	
	// Print message in the middle of the screen
	string msg = "";
	if (winnerColor == "red")
		msg = "!!!!!!!!!!!!!!! CONGRATULATIONS RED, YOU WIN !!!!!!!!!!!!!!!";
	else
		msg = "!!!!!!!!!!!!!!! CONGRATULATIONS BLUE, YOU WIN !!!!!!!!!!!!!!!";
	
	vector<WORD> attr{
		FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE
	};	

	COORD loc;
	loc.X = (WINDOW_WIDTH - (SHORT)msg.size()) / 2;
	loc.Y = WINDOW_HEIGHT / 2;
	DWORD nCharsWritten;
	WriteConsoleOutputCharacterA(hConsoleOutput, msg.c_str(), (DWORD)msg.size(), loc, &nCharsWritten);
	WriteConsoleOutputAttribute(hConsoleOutput, attr.data(), (DWORD)attr.size(), loc, &nCharsWritten);
}

void Console::WinConsole::Draw(int w, int h, WORD color) {
	COORD coord;
	w = w * 3;
	h = h * 2;

	coord.X = (SHORT)w;
	coord.Y = (SHORT)h;

	FillConsoleOutputAttribute(hConsoleOutput, color, 2, coord, &cCharsWritten);
	h++;
	coord.X = (SHORT)w;
	coord.Y = (SHORT)h;

	FillConsoleOutputAttribute(hConsoleOutput, color, 2, coord, &cCharsWritten);
}

void Console::WinConsole::DrawMovesList(int xCoord, int yCoord, WORD color, string msg) {
	vector<WORD> attr{
		color
	};
	screenCoords = COORD{ (SHORT)xCoord, (SHORT)yCoord };
	WriteConsoleOutputCharacterA(hConsoleOutput, msg.c_str(), 5, screenCoords, &cCharsWritten);
//	WriteConsoleOutputAttribute(hConsoleOutput, attr.data(), 5, screenCoords, &cCharsWritten);
	
}

void Console::WinConsole::DrawScores(int xCoord, int yCoord, string msg)
{
	screenCoords = COORD{ (SHORT)xCoord, (SHORT)yCoord };
	WriteConsoleOutputCharacterA(hConsoleOutput, msg.c_str(), (DWORD)msg.size(), screenCoords, &cCharsWritten);

}

void Console::WinConsole::DrawTopBar()
{
	SHORT yCoord = 0;
	SHORT xCoord = 0;
	for (int k = 0; k < 4; k++)
	{
		dwConSize = WINDOW_WIDTH;
		screenCoords = COORD{ xCoord, yCoord };
		CHAR character = ' ';
		FillConsoleOutputCharacterA(hConsoleOutput, character, dwConSize, screenCoords, &cCharsWritten);
		FillConsoleOutputAttribute(hConsoleOutput, BACKGROUND_BLUE|FOREGROUND_WHITE, dwConSize, screenCoords, &cCharsWritten);
		yCoord = yCoord + 1;
	}
}

void Console::WinConsole::DrawRightSide()
{
	SHORT yCoord = 4;
	SHORT xCoord = 44;
	for (int k = 4; k < WINDOW_HEIGHT; k++)
	{
		dwConSize = WINDOW_WIDTH;
		screenCoords = COORD{ xCoord, yCoord };
		CHAR character = ' ';
		//	ScrollWindow(hConsoleOutput, 0, -screenCoords.Y, screenCoords, &cCharsWritten);
		FillConsoleOutputCharacterA(hConsoleOutput, character, dwConSize, screenCoords, &cCharsWritten);
		FillConsoleOutputAttribute(hConsoleOutput, BACKGROUND_WHITE, dwConSize, screenCoords, &cCharsWritten);
		yCoord = yCoord + 1;
	}
}