#include <Python.h>
#include <string.h>

int
main(int argc, char *argv[])
{
	char cmd[1000];
	Py_Initialize();
	sprintf(cmd, "execfile(\"%s\")\n", argv[1]);
	PyRun_SimpleString(cmd);
	Py_Finalize();
	return 0;
}
