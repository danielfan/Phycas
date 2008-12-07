#include <iostream>
#include <string>
#include <windows.h>

// Prompts user with string s, then returns the first character typed. If any problems arise 
// (e.g. cannot obtain handle to console input buffer), bails out by returning the null 
// character (i.e. '\0').
char AskUser(std::string s)
	{
	HANDLE h; 
	DWORD num_chars_read, new_console_mode, prev_console_mode; 
	char char_buffer[1]; 

	// Output the prompt string
	std::cerr << s << std::endl;

	// Get handle to console input buffer
	h = GetStdHandle(STD_INPUT_HANDLE); 
	if (h == INVALID_HANDLE_VALUE) 
		return '\0'; 

	// Save the current input mode (will restore it before we leave this function) 
	if (!GetConsoleMode(h, &prev_console_mode) ) 
		return '\0'; 

	// Set new console mode. There are five mode flags defined in wincon.h (ENABLE_LINE_INPUT, ENABLE_ECHO_INPUT, 
	// ENABLE_PROCESSED_INPUT, ENABLE_WINDOW_INPUT and ENABLE_MOUSE_INPUT), only ENABLE_PROCESSED_INPUT is useful
	// to us, and we specifically want to avoid ENABLE_LINE_INPUT because it requires the user to press the enter
	// key before ReadConsole returns (much better to have this function return the instant the user presses any
	// key).
	new_console_mode = ENABLE_PROCESSED_INPUT;
	if (!SetConsoleMode(h, new_console_mode)) 
		return '\0';

	// Read 1 character and place it in char_buffer. num_chars_read should be 1 afterwards. Note that 
	// the last argument is reserved and must be NULL.
	if (!ReadConsole(h, char_buffer, 1, &num_chars_read, NULL))
		return '\0';

	// Be nice and return console mode to its previous value
	if (!SetConsoleMode(h, prev_console_mode)) 
		return '\0';

	return char_buffer[0];
	}
