#!/usr/bin/python
from read_xcode import *
import cStringIO
fileStart = '''<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>
<?codewarrior exportversion="1.0.1" ideversion="5.0" ?>

<!DOCTYPE PROJECT [
<!ELEMENT PROJECT (TARGETLIST, TARGETORDER, GROUPLIST, DESIGNLIST?)>
<!ELEMENT TARGETLIST (TARGET+)>
<!ELEMENT TARGET (NAME, SETTINGLIST, FILELIST?, LINKORDER?, SEGMENTLIST?, OVERLAYGROUPLIST?, SUBTARGETLIST?, SUBPROJECTLIST?, FRAMEWORKLIST?)>
<!ELEMENT NAME (#PCDATA)>
<!ELEMENT USERSOURCETREETYPE (#PCDATA)>
<!ELEMENT PATH (#PCDATA)>
<!ELEMENT FILELIST (FILE*)>
<!ELEMENT FILE (PATHTYPE, PATHROOT?, ACCESSPATH?, PATH, PATHFORMAT?, ROOTFILEREF?, FILEKIND?, FILEFLAGS?)>
<!ELEMENT PATHTYPE (#PCDATA)>
<!ELEMENT PATHROOT (#PCDATA)>
<!ELEMENT ACCESSPATH (#PCDATA)>
<!ELEMENT PATHFORMAT (#PCDATA)>
<!ELEMENT ROOTFILEREF (PATHTYPE, PATHROOT?, ACCESSPATH?, PATH, PATHFORMAT?)>
<!ELEMENT FILEKIND (#PCDATA)>
<!ELEMENT FILEFLAGS (#PCDATA)>
<!ELEMENT FILEREF (TARGETNAME?, PATHTYPE, PATHROOT?, ACCESSPATH?, PATH, PATHFORMAT?)>
<!ELEMENT TARGETNAME (#PCDATA)>
<!ELEMENT SETTINGLIST ((SETTING|PANELDATA)+)>
<!ELEMENT SETTING (NAME?, (VALUE|(SETTING+)))>
<!ELEMENT PANELDATA (NAME, VALUE)>
<!ELEMENT VALUE (#PCDATA)>
<!ELEMENT LINKORDER (FILEREF*)>
<!ELEMENT SEGMENTLIST (SEGMENT+)>
<!ELEMENT SEGMENT (NAME, ATTRIBUTES?, FILEREF*)>
<!ELEMENT ATTRIBUTES (#PCDATA)>
<!ELEMENT OVERLAYGROUPLIST (OVERLAYGROUP+)>
<!ELEMENT OVERLAYGROUP (NAME, BASEADDRESS, OVERLAY*)>
<!ELEMENT BASEADDRESS (#PCDATA)>
<!ELEMENT OVERLAY (NAME, FILEREF*)>
<!ELEMENT SUBTARGETLIST (SUBTARGET+)>
<!ELEMENT SUBTARGET (TARGETNAME, ATTRIBUTES?, FILEREF?)>
<!ELEMENT SUBPROJECTLIST (SUBPROJECT+)>
<!ELEMENT SUBPROJECT (FILEREF, SUBPROJECTTARGETLIST)>
<!ELEMENT SUBPROJECTTARGETLIST (SUBPROJECTTARGET*)>
<!ELEMENT SUBPROJECTTARGET (TARGETNAME, ATTRIBUTES?, FILEREF?)>
<!ELEMENT FRAMEWORKLIST (FRAMEWORK+)>
<!ELEMENT FRAMEWORK (FILEREF, LIBRARYFILE?, VERSION?)>
<!ELEMENT LIBRARYFILE (FILEREF)>
<!ELEMENT VERSION (#PCDATA)>
<!ELEMENT TARGETORDER (ORDEREDTARGET|ORDEREDDESIGN)*>
<!ELEMENT ORDEREDTARGET (NAME)>
<!ELEMENT ORDEREDDESIGN (NAME, ORDEREDTARGET+)>
<!ELEMENT GROUPLIST (GROUP|FILEREF)*>
<!ELEMENT GROUP (NAME, (GROUP|FILEREF)*)>
<!ELEMENT DESIGNLIST (DESIGN+)>
<!ELEMENT DESIGN (NAME, DESIGNDATA)>
<!ELEMENT DESIGNDATA (#PCDATA)>
]>

<PROJECT>\n\t<TARGETLIST>\n'''

windowsSearchPathSettings = '''				<SETTING>
					<SETTING><NAME>SearchPath</NAME>
						<SETTING><NAME>Path</NAME><VALUE>:MSL:</VALUE></SETTING>
						<SETTING><NAME>PathFormat</NAME><VALUE>MacOS</VALUE></SETTING>
						<SETTING><NAME>PathRoot</NAME><VALUE>CodeWarrior</VALUE></SETTING>
					</SETTING>
					<SETTING><NAME>Recursive</NAME><VALUE>true</VALUE></SETTING>
					<SETTING><NAME>FrameworkPath</NAME><VALUE>false</VALUE></SETTING>
					<SETTING><NAME>HostFlags</NAME><VALUE>All</VALUE></SETTING>
				</SETTING>
				<SETTING>
					<SETTING><NAME>SearchPath</NAME>
						<SETTING><NAME>Path</NAME><VALUE>:Win32-x86 Support:</VALUE></SETTING>
						<SETTING><NAME>PathFormat</NAME><VALUE>MacOS</VALUE></SETTING>
						<SETTING><NAME>PathRoot</NAME><VALUE>CodeWarrior</VALUE></SETTING>
					</SETTING>
					<SETTING><NAME>Recursive</NAME><VALUE>true</VALUE></SETTING>
					<SETTING><NAME>FrameworkPath</NAME><VALUE>false</VALUE></SETTING>
					<SETTING><NAME>HostFlags</NAME><VALUE>All</VALUE></SETTING>
				</SETTING>
'''
macSearchPathSettings = '''					<SETTING>
                        <SETTING><NAME>SearchPath</NAME>
                            <SETTING><NAME>Path</NAME><VALUE>MSL/MSL_C</VALUE></SETTING>
                            <SETTING><NAME>PathFormat</NAME><VALUE>Unix</VALUE></SETTING>
                            <SETTING><NAME>PathRoot</NAME><VALUE>CodeWarrior</VALUE></SETTING>
                        </SETTING>
                        <SETTING><NAME>Recursive</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>FrameworkPath</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>HostFlags</NAME><VALUE>All</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>SearchPath</NAME>
                            <SETTING><NAME>Path</NAME><VALUE>MSL/MSL_C++</VALUE></SETTING>
                            <SETTING><NAME>PathFormat</NAME><VALUE>Unix</VALUE></SETTING>
                            <SETTING><NAME>PathRoot</NAME><VALUE>CodeWarrior</VALUE></SETTING>
                        </SETTING>
                        <SETTING><NAME>Recursive</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>FrameworkPath</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>HostFlags</NAME><VALUE>All</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>SearchPath</NAME>
                            <SETTING><NAME>Path</NAME><VALUE>usr/include</VALUE></SETTING>
                            <SETTING><NAME>PathFormat</NAME><VALUE>Unix</VALUE></SETTING>
                            <SETTING><NAME>PathRoot</NAME><VALUE>OS X Volume</VALUE></SETTING>
                        </SETTING>
                        <SETTING><NAME>Recursive</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>FrameworkPath</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>HostFlags</NAME><VALUE>All</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>SearchPath</NAME>
                            <SETTING><NAME>Path</NAME><VALUE>usr/lib</VALUE></SETTING>
                            <SETTING><NAME>PathFormat</NAME><VALUE>Unix</VALUE></SETTING>
                            <SETTING><NAME>PathRoot</NAME><VALUE>OS X Volume</VALUE></SETTING>
                        </SETTING>
                        <SETTING><NAME>Recursive</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>FrameworkPath</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>HostFlags</NAME><VALUE>All</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>SearchPath</NAME>
                            <SETTING><NAME>Path</NAME><VALUE>MacOS X Support</VALUE></SETTING>
                            <SETTING><NAME>PathFormat</NAME><VALUE>Unix</VALUE></SETTING>
                            <SETTING><NAME>PathRoot</NAME><VALUE>CodeWarrior</VALUE></SETTING>
                        </SETTING>
                        <SETTING><NAME>Recursive</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>FrameworkPath</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>HostFlags</NAME><VALUE>All</VALUE></SETTING>
                    </SETTING>
'''
debuggerRuntimeSettings = '''				<!-- Settings for "Debugger Runtime" panel -->
				<SETTING><NAME>MWRuntimeSettings_WorkingDirectory</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWRuntimeSettings_CommandLine</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWRuntimeSettings_HostApplication</NAME>
					<SETTING><NAME>Path</NAME><VALUE></VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>Generic</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Absolute</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>MWRuntimeSettings_EnvVars</NAME><VALUE></VALUE></SETTING>
'''

def getLibFileGroups(m,w):
	return '''			<GROUP><NAME>Macintosh libraries</NAME>
				<GROUP><NAME>MacOS libraries</NAME>
					<GROUP><NAME>Mach-o</NAME>
						<FILEREF>
							<TARGETNAME>%s</TARGETNAME>
							<PATHTYPE>Name</PATHTYPE>
							<PATH>mwcrt1.o</PATH>
							<PATHFORMAT>Unix</PATHFORMAT>
						</FILEREF>
					</GROUP>
				</GROUP>
				<GROUP><NAME>MSL</NAME>
					<GROUP><NAME>Debug builds</NAME>
						<GROUP><NAME>Mach-o</NAME>
							<FILEREF>
								<TARGETNAME>%s</TARGETNAME>
								<PATHTYPE>Name</PATHTYPE>
								<PATH>console_OS_X.c</PATH>
								<PATHFORMAT>Unix</PATHFORMAT>
							</FILEREF>
							<FILEREF>
								<TARGETNAME>%s</TARGETNAME>
								<PATHTYPE>Name</PATHTYPE>
								<PATH>MSL_All_Mach-O_D.lib</PATH>
								<PATHFORMAT>Unix</PATHFORMAT>
							</FILEREF>
						</GROUP>
					</GROUP>
				</GROUP>
			</GROUP>
			<GROUP><NAME>Win32 libraries</NAME>
				<GROUP><NAME>MSL ANSI libraries</NAME>
					<FILEREF>
						<TARGETNAME>%s</TARGETNAME>
						<PATHTYPE>PathRelative</PATHTYPE>
						<PATHROOT>CodeWarrior</PATHROOT>
						<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
						<PATH>Libraries/Runtime/Libs/MSL_All_x86.lib</PATH>
						<PATHFORMAT>Unix</PATHFORMAT>
					</FILEREF>
				</GROUP>
				<GROUP><NAME>Win32 SDK libraries</NAME>
					<FILEREF>
						<TARGETNAME>%s</TARGETNAME>
						<PATHTYPE>PathRelative</PATHTYPE>
						<PATHROOT>CodeWarrior</PATHROOT>
						<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
						<PATH>Libraries/Win32 SDK/ADVAPI32.LIB</PATH>
						<PATHFORMAT>Unix</PATHFORMAT>
					</FILEREF>
					<FILEREF>
						<TARGETNAME>%s</TARGETNAME>
						<PATHTYPE>PathRelative</PATHTYPE>
						<PATHROOT>CodeWarrior</PATHROOT>
						<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
						<PATH>Libraries/Win32 SDK/COMCTL32.LIB</PATH>
						<PATHFORMAT>Unix</PATHFORMAT>
					</FILEREF>
					<FILEREF>
						<TARGETNAME>%s</TARGETNAME>
						<PATHTYPE>PathRelative</PATHTYPE>
						<PATHROOT>CodeWarrior</PATHROOT>
						<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
						<PATH>Libraries/Win32 SDK/COMDLG32.LIB</PATH>
						<PATHFORMAT>Unix</PATHFORMAT>
					</FILEREF>
					<FILEREF>
						<TARGETNAME>%s</TARGETNAME>
						<PATHTYPE>PathRelative</PATHTYPE>
						<PATHROOT>CodeWarrior</PATHROOT>
						<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
						<PATH>Libraries/Win32 SDK/GDI32.LIB</PATH>
						<PATHFORMAT>Unix</PATHFORMAT>
					</FILEREF>
					<FILEREF>
						<TARGETNAME>%s</TARGETNAME>
						<PATHTYPE>PathRelative</PATHTYPE>
						<PATHROOT>CodeWarrior</PATHROOT>
						<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
						<PATH>Libraries/Win32 SDK/KERNEL32.LIB</PATH>
						<PATHFORMAT>Unix</PATHFORMAT>
					</FILEREF>
					<FILEREF>
						<TARGETNAME>%s</TARGETNAME>
						<PATHTYPE>PathRelative</PATHTYPE>
						<PATHROOT>CodeWarrior</PATHROOT>
						<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
						<PATH>Libraries/Win32 SDK/SHELL32.LIB</PATH>
						<PATHFORMAT>Unix</PATHFORMAT>
					</FILEREF>
					<FILEREF>
						<TARGETNAME>%s</TARGETNAME>
						<PATHTYPE>PathRelative</PATHTYPE>
						<PATHROOT>CodeWarrior</PATHROOT>
						<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
						<PATH>Libraries/Win32 SDK/USER32.LIB</PATH>
						<PATHFORMAT>Unix</PATHFORMAT>
					</FILEREF>
					<FILEREF>
						<TARGETNAME>%s</TARGETNAME>
						<PATHTYPE>PathRelative</PATHTYPE>
						<PATHROOT>CodeWarrior</PATHROOT>
						<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
						<PATH>Libraries/Win32 SDK/WINMM.LIB</PATH>
						<PATHFORMAT>Unix</PATHFORMAT>
					</FILEREF>
					<FILEREF>
						<TARGETNAME>%s</TARGETNAME>
						<PATHTYPE>PathRelative</PATHTYPE>
						<PATHROOT>CodeWarrior</PATHROOT>
						<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
						<PATH>Libraries/Win32 SDK/wsock32.lib</PATH>
						<PATHFORMAT>Unix</PATHFORMAT>
					</FILEREF>
				</GROUP>
			</GROUP>''' % (m, m, m, w, w, w, w, w, w, w, w, w, w)
def getTargetSettingsPanel(targetWindows, name):
	l = targetWindows and 'Win32 x86 Linker' or 'MacOS X PPC Linker'
	return '''				
				<!-- Settings for "Target Settings" panel -->
				<SETTING><NAME>Linker</NAME><VALUE>%s</VALUE></SETTING>
				<SETTING><NAME>PreLinker</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>PostLinker</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>Targetname</NAME><VALUE>%s</VALUE></SETTING>
				<SETTING><NAME>OutputDirectory</NAME>
					<SETTING><NAME>Path</NAME><VALUE>build/mwerks_proj</VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>Unix</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>PHYCAS_ROOT</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>SaveEntriesUsingRelativePaths</NAME><VALUE>true</VALUE></SETTING>
''' % (l, name)
def getFileMappingsPanel(targetWindows):
	return targetWindows and '''				
				<!-- Settings for "File Mappings" panel -->
				<SETTING><NAME>FileMappings</NAME>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.c</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ x86</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.c++</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ x86</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.cc</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ x86</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.cp</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ x86</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.cpp</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ x86</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.h</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ x86</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>true</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.pch</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ x86</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>true</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.pch++</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ x86</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>true</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.rc</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>MW WinRC</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.res</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>WinRes Import</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.txt</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>true</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>W8BN</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.word</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>true</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>WDBN</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.word</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>true</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileType</NAME><VALUE>WDBN</VALUE></SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.word</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>true</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.doc</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>true</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.lib</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>Lib Import x86</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
					<SETTING>
						<SETTING><NAME>FileExtension</NAME><VALUE>.obj</VALUE></SETTING>
						<SETTING><NAME>Compiler</NAME><VALUE>Obj Import x86</VALUE></SETTING>
						<SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
						<SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
						<SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
					</SETTING>
				</SETTING>
				''' or '''                <!-- Settings for "File Mappings" panel -->
                <SETTING><NAME>FileMappings</NAME>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>APPL</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>Appl</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>MMLB</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>Lib Import PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>MPLF</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>Lib Import PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>MWCD</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>RSRC</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.arr</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.bh</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>Balloon Help</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.c</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE>C/C++</VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.c++</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE>C/C++</VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.cc</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE>C/C++</VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.cp</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE>C/C++</VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.cpp</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE>C/C++</VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.exp</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.h</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE>C/C++</VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>true</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.pch</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE>C/C++</VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.pch++</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>MW C/C++ PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE>C/C++</VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.r</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>Rez</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE>Rez</VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>TEXT</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.s</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>PPCAsm</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>XCOF</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>XCOFF Import PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>docu</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>rsrc</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>shlb</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>PEF Import PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileType</NAME><VALUE>stub</VALUE></SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>PEF Import PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.doc</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>true</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.o</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE>XCOFF Import PPC</VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.ppob</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                    <SETTING>
                        <SETTING><NAME>FileExtension</NAME><VALUE>.rsrc</VALUE></SETTING>
                        <SETTING><NAME>Compiler</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>EditLanguage</NAME><VALUE></VALUE></SETTING>
                        <SETTING><NAME>Precompile</NAME><VALUE>false</VALUE></SETTING>
                        <SETTING><NAME>Launchable</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>ResourceFile</NAME><VALUE>true</VALUE></SETTING>
                        <SETTING><NAME>IgnoredByMake</NAME><VALUE>false</VALUE></SETTING>
                    </SETTING>
                </SETTING>
'''

def getConstTargetSettings(includeFile):
	return '''				<!-- Settings for "Build Extras" panel -->
				<SETTING><NAME>CacheModDates</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>DumpBrowserInfo</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>CacheSubprojects</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>UseThirdPartyDebugger</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>BrowserGenerator</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>DebuggerAppPath</NAME>
					<SETTING><NAME>Path</NAME><VALUE></VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>Generic</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Absolute</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>DebuggerCmdLineArgs</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>DebuggerWorkingDir</NAME>
					<SETTING><NAME>Path</NAME><VALUE></VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>Generic</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Absolute</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>CodeCompletionPrefixFileName</NAME><VALUE>Win32Headers.pch</VALUE></SETTING>
				<SETTING><NAME>CodeCompletionMacroFileName</NAME><VALUE>Win32_C++_Macros.h</VALUE></SETTING>
				
				<!-- Settings for "Debugger Target" panel -->
				<SETTING><NAME>ConsoleEncoding</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>LogSystemMessages</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>AutoTargetDLLs</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>StopAtWatchpoints</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>PauseWhileRunning</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>PauseInterval</NAME><VALUE>5</VALUE></SETTING>
				<SETTING><NAME>PauseUIFlags</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>AltExePath</NAME>
					<SETTING><NAME>Path</NAME><VALUE></VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>Generic</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Absolute</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>StopAtTempBPOnLaunch</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>CacheSymbolics</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>TempBPFunctionName</NAME><VALUE>main</VALUE></SETTING>
				<SETTING><NAME>TempBPType</NAME><VALUE>0</VALUE></SETTING>
				
				<!-- Settings for "Remote Debug" panel -->
				<SETTING><NAME>Enabled</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>ConnectionName</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>DownloadPath</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>LaunchRemoteApp</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>RemoteAppPath</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>CoreID</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>JTAGClockSpeed</NAME><VALUE>8000</VALUE></SETTING>
				<SETTING><NAME>IsMultiCore</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>OSDownload</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>UseGlobalOSDownload</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>OSDownloadConnectionName</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>OSDownloadPath</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>AltDownload</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>AltDownloadConnectionName</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "Auto-target" panel -->
				<SETTING><NAME>OtherExecutables</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "Custom Keywords" panel -->
				<SETTING><NAME>CustomColor1</NAME>
					<SETTING><NAME>Red</NAME><VALUE>0</VALUE></SETTING>
					<SETTING><NAME>Green</NAME><VALUE>39321</VALUE></SETTING>
					<SETTING><NAME>Blue</NAME><VALUE>0</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>CustomColor2</NAME>
					<SETTING><NAME>Red</NAME><VALUE>0</VALUE></SETTING>
					<SETTING><NAME>Green</NAME><VALUE>32767</VALUE></SETTING>
					<SETTING><NAME>Blue</NAME><VALUE>0</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>CustomColor3</NAME>
					<SETTING><NAME>Red</NAME><VALUE>0</VALUE></SETTING>
					<SETTING><NAME>Green</NAME><VALUE>32767</VALUE></SETTING>
					<SETTING><NAME>Blue</NAME><VALUE>0</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>CustomColor4</NAME>
					<SETTING><NAME>Red</NAME><VALUE>0</VALUE></SETTING>
					<SETTING><NAME>Green</NAME><VALUE>32767</VALUE></SETTING>
					<SETTING><NAME>Blue</NAME><VALUE>0</VALUE></SETTING>
				</SETTING>
				
				<!-- Settings for "C/C++ Compiler" panel -->
				<SETTING><NAME>MWFrontEnd_C_cplusplus</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_checkprotos</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_arm</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_trigraphs</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_onlystdkeywords</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_enumsalwaysint</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_mpwpointerstyle</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_prefixname</NAME><VALUE>%s</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_ansistrict</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_mpwcnewline</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_wchar_type</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_enableexceptions</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_dontreusestrings</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_poolstrings</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_dontinline</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_useRTTI</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_multibyteaware</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_unsignedchars</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_autoinline</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_booltruefalse</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_inlinelevel</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_ecplusplus</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_objective_c</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_defer_codegen</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_templateparser</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_c99</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWFrontEnd_C_bottomupinline</NAME><VALUE>1</VALUE></SETTING>
				
				<!-- Settings for "C/C++ Warnings" panel -->
				<SETTING><NAME>MWWarning_C_warn_illpragma</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_warn_emptydecl</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_warn_possunwant</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_warn_unusedvar</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_warn_unusedarg</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_warn_extracomma</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_pedantic</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_warningerrors</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_warn_hidevirtual</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_warn_implicitconv</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_warn_notinlined</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWWarning_C_warn_structclass</NAME><VALUE>1</VALUE></SETTING>
				
				<!-- Settings for "FTP Panel" panel -->
				<SETTING><NAME>MWFTP_Post_hostName</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_username</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_password</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_remoteDir</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_ftp_PathVersion</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_ftp_PathType</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_ftp_PathFormat</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_ftp_tree</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_uploadDir</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_ftp_port</NAME><VALUE>21</VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_SendBin</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWFTP_Post_ShouldLog</NAME><VALUE>1</VALUE></SETTING>
				
				<!-- Settings for "Java Command Line" panel -->
				<SETTING><NAME>MWCommandLine_Java_clsName</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWCommandLine_Java_args</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "PJavaDebugging" panel -->
				<SETTING><NAME>MWVJavaDebugging_Protocol</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWVJavaDebugging_JDKVersion</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWVJavaDebugging_TimeOut</NAME><VALUE>15</VALUE></SETTING>
				<SETTING><NAME>MWVJavaDebugging_SupportSlowDevices</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWVJavaDebugging_UseRemoteLaunchAgent</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWVJavaDebugging_LaunchVMasServer</NAME><VALUE>false</VALUE></SETTING>
				
				<!-- Settings for "Java Language" panel -->
				<SETTING><NAME>MWJava_Language_optimize</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_warnDeprecated</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_emitMap</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_strictFileNames</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_strictFileHierarchy</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_1_1_Compatible</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_emitHeaders</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_headerType</NAME><VALUE>JNINativeHeaders</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_packageFilter</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_genComments</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_genHeaders</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_enableAsserts</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_Language_targetVM</NAME><VALUE>1.1</VALUE></SETTING>
				
				<!-- Settings for "Java Manifest-JAD Setting Info" panel -->
				<SETTING><NAME>Manifest-JAD Attributes</NAME>
					<SETTING>
						<SETTING><NAME>Attribute</NAME><VALUE>Main-Class</VALUE></SETTING>
						<SETTING><NAME>Value</NAME><VALUE>Auto-Generated</VALUE></SETTING>
					</SETTING>
				</SETTING>
				
				<!-- Settings for "Java MRJAppBuilder" panel -->
				<SETTING><NAME>MWJava_MRJAppBuilder_outFile</NAME><VALUE>MRJApplication</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_merge</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_quitMenu</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_grow</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdoutType</NAME><VALUE>Console</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stderrType</NAME><VALUE>Console</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdinType</NAME><VALUE>Console</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_appIconPVersion</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_appIconPType</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_appIconPFormat</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_appIconPTree</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_appIconFile</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_splashScreenPVersion</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_splashScreenPType</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_splashScreenPFormat</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_splashScreenPTree</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_splashScreenPICTFile</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_aboutName</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdoutPVersion</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdoutPType</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdoutPFormat</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdoutPTree</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdoutFile</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdoutAppend</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stderrPType</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stderrPFormat</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stderrPTree</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stderrFile</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stderrAppend</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdinPType</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdinPFormat</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdinPTree</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_MRJAppBuilder_stdinFile</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "Java Output" panel -->
				<SETTING><NAME>MWJava_Output_outputtype</NAME><VALUE>JarFile</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_outfile</NAME><VALUE>JavaClasses.jar</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_ftype</NAME><VALUE>1514754080</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_fcreator</NAME><VALUE>1297570384</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_compress</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_genManifest</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_trunctype</NAME><VALUE>Front</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_deleteClasses</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_consoleApp</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_preverify</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_genJad</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_Output_obfuscate</NAME><VALUE>0</VALUE></SETTING>
				
				<!-- Settings for "Java Project" panel -->
				<SETTING><NAME>MWJava_Proj_projtype</NAME><VALUE>Applet</VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_mainClassName</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_HTMLAppCreator</NAME><VALUE>1463898714</VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_HTMLAppName</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_PathVersion</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_PathType</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_PathFormat</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_tree</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_HTMLAppWin32Name</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_compress</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_simulator</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_useVM</NAME><VALUE>macosx</VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_vmarguments</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_vmName</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJava_Proj_simPropFile</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "JavaDoc Project" panel -->
				<SETTING><NAME>MWJavaDoc_Proj_Version</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_Deprecated</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_Author</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_Index</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_Tree</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_SunResolveToSame</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_Shortnames</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_Folder</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_GenerateAPILinks</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_scope</NAME><VALUE>Public</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_fcreator</NAME><VALUE>1297303877</VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_encodingName</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_decodingName</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWJavaDoc_Proj_javaPackagePath</NAME><VALUE>http://java.sun.com/products/jdk/1.1/docs/api/</VALUE></SETTING>
				
				<!-- Settings for "MacOS Merge Panel" panel -->
				<SETTING><NAME>MWMerge_MacOS_projectType</NAME><VALUE>Application</VALUE></SETTING>
				<SETTING><NAME>MWMerge_MacOS_outputName</NAME><VALUE>Merge Out</VALUE></SETTING>
				<SETTING><NAME>MWMerge_MacOS_outputCreator</NAME><VALUE>????</VALUE></SETTING>
				<SETTING><NAME>MWMerge_MacOS_outputType</NAME><VALUE>APPL</VALUE></SETTING>
				<SETTING><NAME>MWMerge_MacOS_suppressWarning</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWMerge_MacOS_copyFragments</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWMerge_MacOS_copyResources</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWMerge_MacOS_flattenResource</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWMerge_MacOS_flatFileName</NAME><VALUE>a.rsrc</VALUE></SETTING>
				<SETTING><NAME>MWMerge_MacOS_flatFileOutputPath</NAME>
					<SETTING><NAME>Path</NAME><VALUE>:</VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>MacOS</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Project</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>MWMerge_MacOS_skipResources</NAME>
					<SETTING><VALUE>DLGX</VALUE></SETTING>
					<SETTING><VALUE>ckid</VALUE></SETTING>
					<SETTING><VALUE>Proj</VALUE></SETTING>
					<SETTING><VALUE>WSPC</VALUE></SETTING>
				</SETTING>
				
				<!-- Settings for "Output Flags" panel -->
				<SETTING><NAME>FileLocked</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>ResourcesMapIsReadOnly</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>PrinterDriverIsMultiFinderCompatible</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>Invisible</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>HasBundle</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>NameLocked</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>Stationery</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>HasCustomIcon</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>Shared</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>HasBeenInited</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>Label</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>Comments</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "Packager Panel" panel -->
				<SETTING><NAME>MWMacOSPackager_UsePackager</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWMacOSPackager_FolderToPackage</NAME>
					<SETTING><NAME>Path</NAME><VALUE></VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>MacOS</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Absolute</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>MWMacOSPackager_CreateClassicAlias</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWMacOSPackager_ClassicAliasMethod</NAME><VALUE>UseTargetOutput</VALUE></SETTING>
				<SETTING><NAME>MWMacOSPackager_ClassicAliasPath</NAME>
					<SETTING><NAME>Path</NAME><VALUE></VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>MacOS</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Absolute</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>MWMacOSPackager_CreatePkgInfo</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWMacOSPackager_PkgCreatorType</NAME><VALUE>????</VALUE></SETTING>
				<SETTING><NAME>MWMacOSPackager_PkgFileType</NAME><VALUE>APPL</VALUE></SETTING>
				
				<!-- Settings for "PPC CodeGen" panel -->
				<SETTING><NAME>MWCodeGen_PPC_structalignment</NAME><VALUE>PPC_mw</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_tracebacktables</NAME><VALUE>Inline</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_processor</NAME><VALUE>P601</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_function_align</NAME><VALUE>4</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_tocdata</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_largetoc</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_profiler</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_vectortocdata</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_poolconst</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_peephole</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_readonlystrings</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_linkerpoolsstrings</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_volatileasm</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_schedule</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_altivec</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_altivec_move_block</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_strictIEEEfp</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_fpcontract</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_genfsel</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_PPC_orderedfpcmp</NAME><VALUE>0</VALUE></SETTING>
				
				<!-- Settings for "PPC CodeGen Mach-O" panel -->
				<SETTING><NAME>MWCodeGen_MachO_structalignment</NAME><VALUE>PPC_mw</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_profiler_enum</NAME><VALUE>Off</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_processor</NAME><VALUE>Generic</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_function_align</NAME><VALUE>4</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_common</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_peephole</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_readonlystrings</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_linkerpoolsstrings</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_volatileasm</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_schedule</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_altivec</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_vecmove</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_fp_ieee_strict</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_fpcontract</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_genfsel</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_MachO_fp_cmps_ordered</NAME><VALUE>0</VALUE></SETTING>
				
				<!-- Settings for "PPC Disassembler" panel -->
				<SETTING><NAME>MWDisassembler_PPC_showcode</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWDisassembler_PPC_extended</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWDisassembler_PPC_mix</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWDisassembler_PPC_nohex</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWDisassembler_PPC_showdata</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWDisassembler_PPC_showexceptions</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWDisassembler_PPC_showsym</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWDisassembler_PPC_shownames</NAME><VALUE>1</VALUE></SETTING>
				
				<!-- Settings for "PPC Global Optimizer" panel -->
				<SETTING><NAME>GlobalOptimizer_PPC_optimizationlevel</NAME><VALUE>Level0</VALUE></SETTING>
				<SETTING><NAME>GlobalOptimizer_PPC_optfor</NAME><VALUE>Speed</VALUE></SETTING>
				
				<!-- Settings for "PPC Linker" panel -->
				<SETTING><NAME>MWLinker_PPC_linksym</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWLinker_PPC_symfullpath</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWLinker_PPC_linkmap</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWLinker_PPC_nolinkwarnings</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWLinker_PPC_dontdeadstripinitcode</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_PPC_permitmultdefs</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_PPC_linkmode</NAME><VALUE>Fast</VALUE></SETTING>
				<SETTING><NAME>MWLinker_PPC_initname</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWLinker_PPC_mainname</NAME><VALUE>__start</VALUE></SETTING>
				<SETTING><NAME>MWLinker_PPC_termname</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "PPC Mac OS X Linker" panel -->
				<SETTING><NAME>MWLinker_MacOSX_linksym</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_symfullpath</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_nolinkwarnings</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_linkmap</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_dontdeadstripinitcode</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_permitmultdefs</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_use_objectivec_semantics</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_strip_debug_symbols</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_split_segs</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_report_msl_overloads</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_objects_follow_linkorder</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_linkmode</NAME><VALUE>Fast</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_exports</NAME><VALUE>ReferencedGlobals</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_sortcode</NAME><VALUE>None</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MacOSX_mainname</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "PPC Mac OS X Project" panel -->
				<SETTING><NAME>MWProject_MacOSX_type</NAME><VALUE>Executable</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_outfile</NAME><VALUE>PAUP 4.0d?? (PPC, test)</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_filecreator</NAME><VALUE>PAUP</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_filetype</NAME><VALUE>APPL</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_vmaddress</NAME><VALUE>4096</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_usedefaultvmaddr</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_flatrsrc</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_flatrsrcfilename</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_flatrsrcoutputdir</NAME>
					<SETTING><NAME>Path</NAME><VALUE>:</VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>MacOS</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Project</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>MWProject_MacOSX_installpath</NAME><VALUE>./</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_dont_prebind</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_flat_namespace</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_frameworkversion</NAME><VALUE>A</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_currentversion</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_MacOSX_flat_oldimpversion</NAME><VALUE>0</VALUE></SETTING>
				
				<!-- Settings for "PPC Mach-O Linker" panel -->
				<SETTING><NAME>MWLinker_MachO_exports</NAME><VALUE>None</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_mainname</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_currentversion</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_compatibleversion</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_symfullpath</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_supresswarnings</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_multisymerror</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_prebind</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_deadstrip</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_objectivecsemantics</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_whichfileloaded</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_whyfileloaded</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_readonlyrelocs</NAME><VALUE>Errors</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_undefinedsymbols</NAME><VALUE>Errors</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_twolevelnamespace</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWLinker_MachO_stripdebugsymbols</NAME><VALUE>0</VALUE></SETTING>
				
				<!-- Settings for "PPC Mach-O Target" panel -->
				<SETTING><NAME>MWProject_MachO_type</NAME><VALUE>Executable</VALUE></SETTING>
				<SETTING><NAME>MWProject_MachO_outfile</NAME><VALUE>a.exe</VALUE></SETTING>
				<SETTING><NAME>MWProject_MachO_filecreator</NAME><VALUE>????</VALUE></SETTING>
				<SETTING><NAME>MWProject_MachO_filetype</NAME><VALUE>MEXE</VALUE></SETTING>
				<SETTING><NAME>MWProject_MachO_vmaddress</NAME><VALUE>4096</VALUE></SETTING>
				<SETTING><NAME>MWProject_MachO_flatrsrc</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWProject_MachO_flatrsrcfilename</NAME><VALUE>a.rsrc</VALUE></SETTING>
				<SETTING><NAME>MWProject_MachO_flatrsrcoutputdir</NAME>
					<SETTING><NAME>Path</NAME><VALUE>:</VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>MacOS</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Project</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>MWProject_MachO_installpath</NAME><VALUE>./</VALUE></SETTING>
				<SETTING><NAME>MWProject_MachO_frameworkversion</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "PPC PEF" panel -->
				<SETTING><NAME>MWPEF_exports</NAME><VALUE>None</VALUE></SETTING>
				<SETTING><NAME>MWPEF_libfolder</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWPEF_sortcode</NAME><VALUE>None</VALUE></SETTING>
				<SETTING><NAME>MWPEF_expandbss</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWPEF_sharedata</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWPEF_olddefversion</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWPEF_oldimpversion</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWPEF_currentversion</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWPEF_fragmentname</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWPEF_collapsereloads</NAME><VALUE>0</VALUE></SETTING>
				
				<!-- Settings for "PPC Project" panel -->
				<SETTING><NAME>MWProject_PPC_type</NAME><VALUE>Application</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_outfile</NAME><VALUE>PAUP 4.0d?? (PPC, test)</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_filecreator</NAME><VALUE>PAUP</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_filetype</NAME><VALUE>APPL</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_size</NAME><VALUE>10000</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_minsize</NAME><VALUE>1024</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_stacksize</NAME><VALUE>64</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_flags</NAME><VALUE>22720</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_symfilename</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_rsrcname</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_rsrcheader</NAME><VALUE>Native</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_rsrctype</NAME><VALUE>????</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_rsrcid</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_rsrcflags</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_rsrcstore</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_rsrcmerge</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_flatrsrc</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_PPC_flatrsrcoutputdir</NAME>
					<SETTING><NAME>Path</NAME><VALUE>:</VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>MacOS</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Project</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>MWProject_PPC_flatrsrcfilename</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "PPCAsm Panel" panel -->
				<SETTING><NAME>MWAssembler_PPC_auxheader</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWAssembler_PPC_symmode</NAME><VALUE>Mac</VALUE></SETTING>
				<SETTING><NAME>MWAssembler_PPC_dialect</NAME><VALUE>PPC</VALUE></SETTING>
				<SETTING><NAME>MWAssembler_PPC_prefixfile</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWAssembler_PPC_typecheck</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWAssembler_PPC_warnings</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWAssembler_PPC_casesensitive</NAME><VALUE>0</VALUE></SETTING>
				
				<!-- Settings for "Property List" panel -->
				<SETTING><NAME>PList_OutputType</NAME><VALUE>File</VALUE></SETTING>
				<SETTING><NAME>PList_OutputEncoding</NAME><VALUE>UTF-8</VALUE></SETTING>
				<SETTING><NAME>PList_Prefix</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>PList_FileFilename</NAME><VALUE>Info.plist</VALUE></SETTING>
				<SETTING><NAME>PList_FileDirectory</NAME>
					<SETTING><NAME>Path</NAME><VALUE>:</VALUE></SETTING>
					<SETTING><NAME>PathFormat</NAME><VALUE>MacOS</VALUE></SETTING>
					<SETTING><NAME>PathRoot</NAME><VALUE>Project</VALUE></SETTING>
				</SETTING>
				<SETTING><NAME>PList_ResourceType</NAME><VALUE>plst</VALUE></SETTING>
				<SETTING><NAME>PList_ResourceID</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>PList_ResourceName</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "Rez Compiler" panel -->
				<SETTING><NAME>MWRez_Language_maxwidth</NAME><VALUE>80</VALUE></SETTING>
				<SETTING><NAME>MWRez_Language_script</NAME><VALUE>Roman</VALUE></SETTING>
				<SETTING><NAME>MWRez_Language_alignment</NAME><VALUE>Align1</VALUE></SETTING>
				<SETTING><NAME>MWRez_Language_filtermode</NAME><VALUE>FilterSkip</VALUE></SETTING>
				<SETTING><NAME>MWRez_Language_suppresswarnings</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWRez_Language_escapecontrolchars</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWRez_Language_prefixname</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWRez_Language_filteredtypes</NAME><VALUE>'CODE' 'DATA' 'PICT'</VALUE></SETTING>
				
				<!-- Settings for "WinRC Compiler" panel -->
				<SETTING><NAME>MWWinRC_prefixname</NAME><VALUE></VALUE></SETTING>
				
				<!-- Settings for "x86 CodeGen" panel -->
				<SETTING><NAME>MWCodeGen_X86_processor</NAME><VALUE>PentiumII</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_alignment</NAME><VALUE>bytes8</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_exceptions</NAME><VALUE>ZeroOverhead</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_name_mangling</NAME><VALUE>MWWin32</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_use_extinst</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_extinst_mmx</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_extinst_3dnow</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_use_mmx_3dnow_convention</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_extinst_cmov</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_extinst_sse</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_extinst_sse2</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_intrinsics</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_optimizeasm</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_disableopts</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_profile</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_runtime</NAME><VALUE>Custom</VALUE></SETTING>
				<SETTING><NAME>MWCodeGen_X86_readonlystrings</NAME><VALUE>0</VALUE></SETTING>
				
				<!-- Settings for "x86 COFF" panel -->
				<SETTING><NAME>MWLinker_X86_subsysmajorid</NAME><VALUE>4</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_subsysminorid</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCOFF_X86_opsysmajorid</NAME><VALUE>4</VALUE></SETTING>
				<SETTING><NAME>MWCOFF_X86_opsysminorid</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_usrmajorid</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_usrminorid</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_maxstacksize</NAME><VALUE>1024</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_minstacksize</NAME><VALUE>4</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_size</NAME><VALUE>1024</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_minsize</NAME><VALUE>4</VALUE></SETTING>
				<SETTING><NAME>MWCOFF_X86_coff_flags</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWCOFF_X86_dll_flags</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_baseaddress</NAME><VALUE>4194304</VALUE></SETTING>
				<SETTING><NAME>MWCOFF_X86_filealign</NAME><VALUE>512</VALUE></SETTING>
				<SETTING><NAME>MWCOFF_X86_sectionalign</NAME><VALUE>4096</VALUE></SETTING>
				
				<!-- Settings for "x86 Disassembler" panel -->
				<SETTING><NAME>PDisasmX86_showHeaders</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showSectHeaders</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showSymTab</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showCode</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showData</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showDebug</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showExceptions</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showRaw</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showAllRaw</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showSource</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showRelocation</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showHex</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showComments</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_showSymDefs</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_unmangle</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_verbose</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_resolveRelocs</NAME><VALUE>true</VALUE></SETTING>
				<SETTING><NAME>PDisasmX86_resolveLocals</NAME><VALUE>false</VALUE></SETTING>
				
				<!-- Settings for "x86 Exceptions Panel" panel -->
				<SETTING><NAME>MWDebugger_X86_Exceptions</NAME>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
					<SETTING><VALUE>0</VALUE></SETTING>
				</SETTING>
				
				<!-- Settings for "x86 Global Optimizer" panel -->
				<SETTING><NAME>GlobalOptimizer_X86_optimizationlevel</NAME><VALUE>Level4</VALUE></SETTING>
				<SETTING><NAME>GlobalOptimizer_X86_optfor</NAME><VALUE>Speed</VALUE></SETTING>
				
				<!-- Settings for "x86 Linker" panel -->
				<SETTING><NAME>MWLinker_X86_entrypointusage</NAME><VALUE>Default</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_entrypoint</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_subsystem</NAME><VALUE>WinCUI</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_commandfile</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_generatemap</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_linksym</NAME><VALUE>0</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_linkCV</NAME><VALUE>1</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_symfullpath</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_linkdebug</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_checksum</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_zero_init_bss</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_mergedata</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_usedefaultlibs</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_adddefaultlibs</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_nowarnings</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWLinker_X86_verbose</NAME><VALUE>false</VALUE></SETTING>
				
				<!-- Settings for "x86 Project" panel -->
				<SETTING><NAME>MWProject_X86_type</NAME><VALUE>Application</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_outfile</NAME><VALUE>d_PhycasMW-console.EXE</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_importlib</NAME><VALUE></VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_setimportlibdir</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_dontgenerateimportlib</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_oldformatlib</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_replaceobjextension</NAME><VALUE>false</VALUE></SETTING>
				<SETTING><NAME>MWProject_X86_copyallfiles</NAME><VALUE>false</VALUE></SETTING>
''' % includeFile
def getLibFiles(targetWindows):
	if targetWindows: return '''				<FILE>
					<PATHTYPE>PathRelative</PATHTYPE>
					<PATHROOT>CodeWarrior</PATHROOT>
					<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
					<PATH>Libraries/Win32 SDK/KERNEL32.LIB</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Unknown</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
				<FILE>
					<PATHTYPE>PathRelative</PATHTYPE>
					<PATHROOT>CodeWarrior</PATHROOT>
					<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
					<PATH>Libraries/Win32 SDK/USER32.LIB</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Unknown</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
				<FILE>
					<PATHTYPE>PathRelative</PATHTYPE>
					<PATHROOT>CodeWarrior</PATHROOT>
					<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
					<PATH>Libraries/Win32 SDK/COMCTL32.LIB</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Unknown</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
				<FILE>
                    <PATHTYPE>PathRelative</PATHTYPE>
                    <PATHROOT>CodeWarrior</PATHROOT>
                    <ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
                    <PATH>Libraries/Win32 SDK/wsock32.lib</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                    <FILEKIND>Library</FILEKIND>
                    <FILEFLAGS>Debug</FILEFLAGS>
                </FILE>
                <FILE>
					<PATHTYPE>PathRelative</PATHTYPE>
					<PATHROOT>CodeWarrior</PATHROOT>
					<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
					<PATH>Libraries/Win32 SDK/GDI32.LIB</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Unknown</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
				<FILE>
					<PATHTYPE>PathRelative</PATHTYPE>
					<PATHROOT>CodeWarrior</PATHROOT>
					<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
					<PATH>Libraries/Win32 SDK/COMDLG32.LIB</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Unknown</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
				<FILE>
					<PATHTYPE>PathRelative</PATHTYPE>
					<PATHROOT>CodeWarrior</PATHROOT>
					<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
					<PATH>Libraries/Win32 SDK/SHELL32.LIB</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Unknown</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
				<FILE>
					<PATHTYPE>PathRelative</PATHTYPE>
					<PATHROOT>CodeWarrior</PATHROOT>
					<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
					<PATH>Libraries/Win32 SDK/ADVAPI32.LIB</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Unknown</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
				<FILE>
					<PATHTYPE>PathRelative</PATHTYPE>
					<PATHROOT>CodeWarrior</PATHROOT>
					<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
					<PATH>Libraries/Win32 SDK/WINMM.LIB</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Unknown</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
				<FILE>
					<PATHTYPE>PathRelative</PATHTYPE>
					<PATHROOT>CodeWarrior</PATHROOT>
					<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
					<PATH>Libraries/Runtime/Libs/MSL_All_x86.lib</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Unknown</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
'''
	else: return '''				<FILE>
					<PATHTYPE>Name</PATHTYPE>
					<PATH>MSL_All_Mach-O_D.lib</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Library</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
				<FILE>
					<PATHTYPE>Name</PATHTYPE>
					<PATH>console_OS_X.c</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Text</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
				<FILE>
					<PATHTYPE>Name</PATHTYPE>
					<PATH>mwcrt1.o</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Library</FILEKIND>
					<FILEFLAGS>Debug</FILEFLAGS>
				</FILE>
'''
def getLibLinkFileRef(targetWindows):
	if targetWindows: return '''                <FILEREF>
                    <PATHTYPE>PathRelative</PATHTYPE>
                    <PATHROOT>CodeWarrior</PATHROOT>
                    <ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
                    <PATH>Libraries/Win32 SDK/KERNEL32.LIB</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
                <FILEREF>
                    <PATHTYPE>PathRelative</PATHTYPE>
                    <PATHROOT>CodeWarrior</PATHROOT>
                    <ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
                    <PATH>Libraries/Win32 SDK/USER32.LIB</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
                <FILEREF>
                    <PATHTYPE>PathRelative</PATHTYPE>
                    <PATHROOT>CodeWarrior</PATHROOT>
                    <ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
                    <PATH>Libraries/Win32 SDK/COMCTL32.LIB</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
                <FILEREF>
                    <PATHTYPE>PathRelative</PATHTYPE>
                    <PATHROOT>CodeWarrior</PATHROOT>
                    <ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
                    <PATH>Libraries/Win32 SDK/GDI32.LIB</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
                <FILEREF>
                    <PATHTYPE>PathRelative</PATHTYPE>
                    <PATHROOT>CodeWarrior</PATHROOT>
                    <ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
                    <PATH>Libraries/Win32 SDK/COMDLG32.LIB</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
                <FILEREF>
                    <PATHTYPE>PathRelative</PATHTYPE>
                    <PATHROOT>CodeWarrior</PATHROOT>
                    <ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
                    <PATH>Libraries/Win32 SDK/SHELL32.LIB</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
                <FILEREF>
                    <PATHTYPE>PathRelative</PATHTYPE>
                    <PATHROOT>CodeWarrior</PATHROOT>
                    <ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
                    <PATH>Libraries/Win32 SDK/ADVAPI32.LIB</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
                <FILEREF>
                    <PATHTYPE>PathRelative</PATHTYPE>
                    <PATHROOT>CodeWarrior</PATHROOT>
                    <ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
                    <PATH>Libraries/Win32 SDK/WINMM.LIB</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
                <FILEREF>
                    <PATHTYPE>PathRelative</PATHTYPE>
                    <PATHROOT>CodeWarrior</PATHROOT>
                    <ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
                    <PATH>Libraries/Runtime/Libs/MSL_All_x86.lib</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
				<FILEREF>
					<TARGETNAME>%s</TARGETNAME>
					<PATHTYPE>PathRelative</PATHTYPE>
					<PATHROOT>CodeWarrior</PATHROOT>
					<ACCESSPATH>:Win32-x86 Support:</ACCESSPATH>
					<PATH>Libraries/Win32 SDK/wsock32.lib</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
				</FILEREF>
'''
	else: return '''                <FILEREF>
                    <PATHTYPE>Name</PATHTYPE>
                    <PATH>MSL_All_Mach-O_D.lib</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
                <FILEREF>
                    <PATHTYPE>Name</PATHTYPE>
                    <PATH>mwcrt1.o</PATH>
                    <PATHFORMAT>Unix</PATHFORMAT>
                </FILEREF>
'''
def getRelativeFilePathElement(sourceTree, relPath):
	pathType = 'PathRelative' # relPath.rfind('/') < 0 and 'Name' or 'PathRelative'
	if relPath.endswith('.cpp') or relPath.endswith('.c'):
		dbStr = 'Debug'
	else:
		dbStr = ''
	return '''				<FILE>
					<PATHTYPE>%s</PATHTYPE>
					<PATHROOT>%s</PATHROOT>
					<ACCESSPATH></ACCESSPATH>
					<PATH>%s</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
					<FILEKIND>Text</FILEKIND>
					<FILEFLAGS>%s</FILEFLAGS>
				</FILE>
''' % (pathType, sourceTree, relPath, dbStr)
	
def getLinkOrderFileRefElement(sourceTree, relPath):
	pathType = 'PathRelative' # relPath.rfind('/') < 0 and 'Name' or 'PathRelative'
	return '''				<FILEREF>
					<PATHTYPE>%s</PATHTYPE>
					<PATHROOT>%s</PATHROOT>
					<ACCESSPATH></ACCESSPATH>
					<PATH>%s</PATH>
					<PATHFORMAT>Unix</PATHFORMAT>
				</FILEREF>
''' % (pathType, sourceTree, relPath)
	
def composeFileElement(p):
	splitPath = p.split('/')
	return getRelativeFilePathElement(splitPath[0], '/'.join(splitPath[1:]))
def composeLinkOrderFileRef(p):
	splitPath = p.split('/')
	return getLinkOrderFileRefElement(splitPath[0], '/'.join(splitPath[1:]))
def mwerksSetting(settingName, val = ''):
	return xmlWrap('SETTING', xmlWrap('NAME', settingName) + xmlWrap('VALUE', val))
def mwerksStartSetting(settingName):
	return '<SETTING>' + xmlWrap('NAME', settingName)
def mwerksSearchPath(name, pathRoot, pathFmt):
	tabstr = '\t\t\t\t\t'
	return ('<SETTING>\n\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n\t%s</SETTING>\n\t%s\n\t%s\n\t%s\n%s</SETTING>' % (
		tabstr + mwerksStartSetting('SearchPath'),
		tabstr + mwerksSetting('Path', name), 
		tabstr + mwerksSetting('PathFormat', pathFmt), 
		tabstr + mwerksSetting('PathRoot', pathRoot), 
		tabstr,
		tabstr + mwerksSetting('Recursive', 'true'), 
		tabstr + mwerksSetting('FrameworkPath', 'false'),
		tabstr + mwerksSetting('HostFlags', 'All'),
		tabstr))


class MWerksSourceFileTree(SourceFileTree):
		if False:
			cpi = len(commonPrefix)
			fileChildren = []
			dirChildren = {}
			for f in listOfFiles:
				p = f[cpi:].split('/')
				if len(p) == 1: 
					fileChildren.append(f)
				elif dirChildren.get(p[0]) == None:
					dirChildren[p[0]] = [f]
				else:
					dirChildren[p[0]].append(f)
			for f in fileChildren:
				self.writeFlatVCFile(out, f, targets, indentLevel)
			if createAllDirs or len(dirChildren) > 1 or len(fileChildren) > 1:
				newIndentLevel = indentLevel + '\t'
				keys = [k for k in dirChildren.iterkeys()]
				keys.sort()
				for k in keys:
					out.write('%s<Filter\n%sName="%s"\n%sFilter="">\n' % (indentLevel, newIndentLevel, k, newIndentLevel))
					self.writeHierarchichalVCFiles(out, targets, dirChildren[k], commonPrefix + k + '/', newIndentLevel, True)
					out.write('%s</Filter>\n' % indentLevel)

class MWerksTarget(PBXNativeTarget):
	def getMWerksTargetName(self, targetWindows):
		return (targetWindows and 'win_' or 'mac_') + self.buildSettings.PRODUCT_NAME + self.buildSettings.buildStyleName
	def writeMWerksTarget(self, out, targetWindows):
		self.mwname = self.getMWerksTargetName(targetWindows)
		out.write('\t\t<TARGET>\n')
		out.write('\t\t\t%s\n' % xmlWrap('NAME', self.mwname))
		self.writeMWerksSettings(out, targetWindows)
		self.writeMWerksFileList(out, targetWindows)
		self.writeMWerksLinkOrder(out, targetWindows)
		out.write('\t\t</TARGET>\n')
	def writeMWerksSettings(self, out, targetWindows):
		out.write('\t\t\t<SETTINGLIST>\n\n')
		out.write('\t\t\t\t<!-- Settings for "Source Trees" panel -->\n')
		out.write('\t\t\t\t%s\n\n' % mwerksSetting('UserSourceTrees'))
		out.write('\t\t\t\t<!-- Settings for "Access Paths" panel -->\n')
		out.write('\t\t\t\t%s\n' % mwerksSetting('AlwaysSearchUserPaths', 'true'))
		out.write('\t\t\t\t%s\n' % mwerksSetting('InterpretDOSAndUnixPaths', 'true'))
		out.write('\t\t\t\t%s\n' % mwerksSetting('RequireFrameworkStyleIncludes', 'false'))
		out.write('\t\t\t\t%s\n' % mwerksStartSetting('UserSearchPaths'))
		out.write('\t\t\t\t\t%s\n' % mwerksSearchPath('', 'PHYCAS_ROOT', 'UNIX'))
		out.write('\t\t\t\t\t%s\n' % mwerksSearchPath('', 'NCL_ROOT', 'UNIX'))
		out.write('\t\t\t\t</SETTING>\n')
		out.write('\t\t\t\t%s\n' % mwerksStartSetting('SystemSearchPaths'))
		if targetWindows: out.write(windowsSearchPathSettings)
		else:	out.write(macSearchPathSettings)
		out.write('\t\t\t\t\t%s\n' % mwerksSearchPath('', 'BOOST_ROOT', 'UNIX'))
		out.write('\t\t\t\t</SETTING>\n')
		out.write(debuggerRuntimeSettings)
		out.write(getTargetSettingsPanel(targetWindows, self.mwname))
		out.write(getFileMappingsPanel(targetWindows))
		
		out.write(getConstTargetSettings(self.buildSettings.GCC_PREFIX_HEADER))
		out.write('\t\t\t</SETTINGLIST>\n')
	def writeMWerksFileList(self, out, targetWindows):
		out.write('\t\t\t<FILELIST>\n')
		out.write(getLibFiles(targetWindows))
		fp =  self.getFilePaths()
		for f in fp:
			out.write(composeFileElement(f))
		out.write('\t\t\t</FILELIST>\n')
	def writeMWerksLinkOrder(self, out, targetWindows):
		out.write('\t\t\t<LINKORDER>\n')
		out.write(getLibLinkFileRef(targetWindows))
		sourceFiles = self.getSources()
		for f in sourceFiles:
			out.write(composeLinkOrderFileRef(f))
		out.write('\t\t\t</LINKORDER>\n')
		
class MWerksProject(PBXProject):
	def writeMWerksXMLProj(self, out, allFiles):
		global fileStart
		out.write(fileStart)
		for t in self.targets:
			t.writeMWerksTarget(out, True)
		for t in self.targets:
			t.writeMWerksTarget(out, False)
		out.write('\t\t</TARGETLIST>\n\t\t<TARGETORDER>\n')
		for t in self.targets: out.write('\t\t\t%s\n' % xmlWrap('ORDEREDTARGET', xmlWrap('NAME', t.getMWerksTargetName(True))))
		for t in self.targets: out.write('\t\t\t%s\n' % xmlWrap('ORDEREDTARGET', xmlWrap('NAME', t.getMWerksTargetName(False))))
		out.write('\t\t</TARGETORDER>\n')
		aMacName, aWinName = (self.targets[0].getMWerksTargetName(False), self.targets[0].getMWerksTargetName(True))
		out.write('\t\t<GROUPLIST>\n')
		for st in allFiles:
			st.writeMWerksFileGroups(out, self.targets)
		out.write('\t\t\t%s\n' % getLibFileGroups(aMacName, aWinName))
		out.write('\t\t</GROUPLIST>\n\n</PROJECT>\n')

if __name__ == '__main__':
	if len(sys.argv) < 2:	print 'Usage: %s xcode_project {Makefile_name}' % sys.argv[0]; sys.exit(1)
	xcodeFileName = sys.argv[1]
	typeToTypeMap = { PBXProject : MWerksProject, PBXNativeTarget : MWerksTarget}
	project, allFiles = readXCodeProject(xcodeFileName, [], {}, typeToTypeMap)
	out = cStringIO.StringIO()
	project.writeMWerksXMLProj(out, allFiles)
	print out.getvalue()
			
