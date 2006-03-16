PhycasExecutor - a utility to run the Phycas main process and user interface.

Setup
-----
**Initially:

1) Download the Java JDK from java.sun.com (avoid getting the Enterprise Edition, or EE, version
   which contains much more than you need)

2) Download the boost library from http://www.boost.org/ and unpack it in the parent directory
   of phycasdev (so that boost and phycasdev are siblings)

3) Download ant from ant.apache.org (put in directory with no spaces in name), set environment variables 

ANT_HOME={ant install dir}
JAVA_HOME={jdk or jre dir} *may already be set

and add ANT_HOME/bin to path.

4) Download swixml from www.swixml.org and copy swixml.jar from {swixml install dir}/build and jdom.jar from {swixml install dir}/lib to {phycasdev path}/gui directory.

5) Download XMLBeans from xmlbeans.apache.org  (put in directory with no spaces in name), set environment variables

XMLBEANS_HOME={xmlbeans install dir}

and add XMLBEANS_HOME/bin to path. Then, copy xbean.jar from {xmlbeans install dir}/lib to 
{phycasdev path}/gui directory.

6) Download jfreechart & jcommon from www.jfree.org/jfreechart, copy jfreechart_[version#].jar to {phycasdev path}/gui directory then rename jfreechart.jar, and copy jcommon_[version#].jar to {phycasdev path}/gui directory then rename jcommon.jar.

7) Download JavaHelp from java.sun.com/products/javahelp/download_binary.html (put in directory with no spaces in name), set environment variable

JAVAHELP_HOME={javaHelp install dir}

and add JAVAHELP_HOME/javahelp/bin to path.  Then, copy jh.jar to {phycasdev path}/gui directory.

**Initially, and again when the schema changes:

Note: the following steps have been largely automated for Windows by creation of the 
renew.bat batch file.

1) Run the following commands to create the XMLBeans jar files:

cd {phycasdev path}/gui directory
scomp -out commandLanguage.jar ..\command_archive\xml\phycas_command_language.xsd
scomp -out commandState.jar ..\command_archive\xml\command_state.xsd

Under Windows (and perhaps other platforms), if you get a CreateProcess exception, be 
sure the JDK bin directory comes first in PATH (in System Variables) because scomp is 
trying to locate the java compiler in the Java runtime directory instead of the JDK bin
directory. 

Build & Run
-----------
1) cd to {phycasdev path} dir 
2) exec "python transform_cmd_xml.py"
3) Build the cpp program and move .exe to {phycasdev path}/executor dir and rename PhycasSocket.exe
4) cd to {phycasdev path}/gui/help
5) exec "jhindexer master"
6) cd to {phycasdev path}/executor
7) exec "ant -f {phycasdev path}/build/buildGui.xml"
8) exec "ant -f {phycasdev path}/build/buildExec.xml"
9) exec "java -jar PhycasExec.jar" or double-click on icon

*a properties file called "phycas.properties" will be stored in the user's home directory 
(C:\Documents and Settings\<userid>\phycas directory in Windows) upon first running 
the executor containing the phycas specific properties:

phycas.home={phycasdev/executor dir full path}
phycas.jre_home={jre dir full path}

