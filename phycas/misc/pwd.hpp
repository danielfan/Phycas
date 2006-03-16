#if !defined(PHO_PWD_H)
#define PHO_PWD_H

class NxsCommandManager;

class PWD
	{
	public:
		void SetupPWD(NxsCommandManager *cmdMgr);
		CmdResult HandlePWD();
	private:
		PWD();
		friend PWDCreator;
	};

#endif

