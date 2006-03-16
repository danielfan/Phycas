#! /usr/bin/python2.2
from sax_constructible import SAXConstructible, ChildrenToHandle, TextOnlyElement
import sys, cStringIO, re, string, os
from command_param_ncl_transform import *
from sets import Set
from phycas_sax_util import *
from os import path
wsToUnderscoreTable = string.maketrans(' ', '_')
dotToUnderscoreTable = string.maketrans('.', '_')
allchars = string.maketrans('','')

def writeCInclude(out, s):
	if len(s) > 0:
		if s[0] == '<':
			out.write('#include %s\n' % s)
		else:
			out.write('#include "%s"\n' % s)


class PhycasCmdSetupInfo(PhycasCallback): pass

class CmdPhycasImpl(SAXConstructible):
	def getSetupRequiredHeader(self):
		if self.setup:
			return self.setup.getSetupRequiredHeader()
		return []
	def getSetupFunctionSig(self):
		return 'void %s::%s(NxsCommandManager *cmdMgr)' % (self.setup.getReceiverClass(), self.setup.getFuncName())
	def getSetupReceiverClass(self):
		return self.setup.getReceiverClass()
	def getSetupFuncName(self):
		return self.setup.getFuncName()
	def getHandlerReceiverClass(self):
		return self.handler.getReceiverClass()
	def getHandlerFuncName(self):
		return self.handler.getFuncName()

class CmdParamGroup(SAXConstructible): pass

class Command(SAXConstructible):
	def getVarAddress(self, paramLabel):
		for cp in self.cmd_param:
			if cp.getLabel() == paramLabel: return '&settingsStruct->' + cp.getCPPVarLabel()
		s = 'Parameter Label %s was not found in Command.getVarAddress()' % paramLabel
		raise ValueError, s
	def getCommandName(self):
		return self.label # automaticallyConvertToString = fale -> self.label.getRawChars()
	def endSelfElement(self, name):
		self.cmd_param = []
		for cmdParamGroup in self.cmd_param_group:
			self.cmd_param.extend(cmdParamGroup.cmd_param)
		for cmdParam in self.cmd_param:
			cmdParam.setCommandParent(self)
	def getSettingsStructFilename(self):
		global wsToUnderscoreTable
		return self.getCommandName().translate(wsToUnderscoreTable).lower() + '_settings.hpp'
	def getSetupFilename(self): 
		global wsToUnderscoreTable
		return self.getCommandName().translate(wsToUnderscoreTable).lower() + '_setup.cpp'
	def writeOutput(self, tmpDir):
		if len(self.cmd_param) > 0:
			self.writeSettingsHeader(tmpDir)
		self.writeSetupFile(tmpDir)
	def getSetupFunctionSig(self):
		return self.phycas_impl.getSetupFunctionSig()
	def getSetupCmdInstancePtrName(self):
		return 'v' +  self.getCommandName().translate(allchars, '_') + 'Command'
	def writeTests(self, out):
		tests = self.availability + self.requirement
		for t in tests:
			t.writeTestDefinition(out)
		for a in self.availability:
			out.write('\t%s->AddCommandAvailabilityReq(%s);\n' % (self.getSetupCmdInstancePtrName(), a.getNCLTestInstanceName()))
		for r in self.requirement:
			out.write('\t%s->AddExecutionReq(%s);\n' % (self.getSetupCmdInstancePtrName(), r.getNCLTestInstanceName()))
	def hasTests(self):
		return len(self.availability) > 0 or len(self.requirement) > 0
	def composeSetupFile(self, out): 
		requiredHeaders = Set(self.phycas_impl.getSetupRequiredHeader())
		if self.hasTests():
			requiredHeaders.add('ncl/misc/nxs_test_impl.hpp')
			requiredHeaders.add('phycas/taxa/taxa_manager.hpp') # //@ TEMP
		for cmdParam in self.cmd_param:
			requiredHeaders.update(cmdParam.getNCLCmdParamRequiredHeaders())
		settingsStructName = self.getSettingStructName()
		underscoredCmdName = self.getCommandName().translate(allchars, '_')
		callbackTypedef = 'Phycas' + settingsStructName+ 'CmdCallback'
		settingsPtrTypedef = settingsStructName + 'Ptr'
		cmdTypedef = 'Phycas' + underscoredCmdName + 'Command'
		cmdPtrTypedef = cmdTypedef + 'Ptr'
		vCmdPtr = self.getSetupCmdInstancePtrName()
		vCmdCallback = 'v' + underscoredCmdName + 'CmdCallback'
		#	out = cStringIO.StringIO()
		#	out.write('#if !defined (OLD_XML_CMDS_PROCESSING)\n')
		out.write('#include "phycas/force_include.h"\n')
		if self.compilation_target != 'ALL':
			cond_comp_statement =  self.compilation_target[0] == '!' and ('if !defined(%s)\n' % self.compilation_target[1:]) or ('if defined(%s)\n' % self.compilation_target)
			out.write('#%s\n' % cond_comp_statement)
		out.write('#include "ncl/nxs_defs.hpp"\n#include "ncl/command/nxs_auto_command.hpp"\n')
		if len(self.cmd_param) > 0:
			out.write('#include "phycas/command/%s"\n' % self.getSettingsStructFilename())
		for i in requiredHeaders:
			writeCInclude(out, i)
		out.write('\n%s\n' % self.getSetupFunctionSig())
		out.write('\t{\n')
		if len(self.cmd_param) > 0:
			out.write('\ttypedef NxsExecuteCommandCallback<%s, %s> %s;\n' % (self.phycas_impl.getHandlerReceiverClass(), settingsStructName, callbackTypedef))
			out.write('\ttypedef boost::shared_ptr<%s> %s;\n' % (settingsStructName, settingsPtrTypedef))
			out.write('\t%s settingsStruct = %s(new %s);\n' %(settingsPtrTypedef, settingsPtrTypedef, settingsStructName))
			out.write('\t%s %s = %s(this, &%s::%s, settingsStruct);\n' % (callbackTypedef, vCmdCallback, callbackTypedef, self.phycas_impl.getHandlerReceiverClass(), self.phycas_impl.getHandlerFuncName()))
		else:
			out.write('\ttypedef NxsVoidExecuteCommandCallback<%s> %s;\n' % (self.phycas_impl.getHandlerReceiverClass(), callbackTypedef))
			out.write('\t%s v%sCmdCallback = %s(this, &%s::%s);\n' % (callbackTypedef, underscoredCmdName, callbackTypedef, self.phycas_impl.getHandlerReceiverClass(), self.phycas_impl.getHandlerFuncName()))
		out.write('\n\ttypedef NxsAutoCommand<%s> %s;\n' %(callbackTypedef, cmdTypedef))
		out.write('\ttypedef boost::shared_ptr<%s> %s;\n' %(cmdTypedef, cmdPtrTypedef))
		out.write('\n\t%s %s = %s(new %s("%s", %s));\n' % (cmdPtrTypedef, vCmdPtr, cmdPtrTypedef, cmdTypedef, underscoredCmdName, vCmdCallback))
		out.write('\t%s->SetDescription("%s");\n' % (vCmdPtr, str(self.description)))
		for cmdParam in self.cmd_param:
			cmdParam.writeInstantiateNCLParamObj(out)
		for cmdParam in self.cmd_param:
			n = cmdParam.getSetupInstanceName()
			funcToCall = len(cmdParam.getLabel()) > 0 and 'AddKeyword' or 'AddUnnamedSetting'
			out.write('\t%s->%s(NxsCmdOptionShPtr(%s), false);\n' %(vCmdPtr, funcToCall,  n))
			if funcToCall == 'AddUnnamedSetting' and cmdParam.followWithEquals(): # hack to accomodate definition commands in NCL (e.g. TaxSet)
				out.write('\t%s->ExpectEqualsSign();\n' % vCmdPtr)
		for cmdParam in self.cmd_param:
			d = cmdParam.getDescription()
			if len(d) > 0:
				n = cmdParam.getSetupInstanceName()
				out.write('\t%s->SetDescription("%s");\n' % (n, d))
		self.writeTests(out)
		out.write('\tcmdMgr->AddCommand(%s);\n' % vCmdPtr)
		out.write('\t}\n')
		if self.compilation_target != 'ALL':
			out.write('#endif //%s\n' % cond_comp_statement)
		#	out.write('#endif //if !defined (OLD_XML_CMDS_PROCESSING)\n')
	def getSettingStructName(self):
		return self.getCommandName().translate(allchars, '_') + 'Settings'
	def writeSettingsHeader(self, tmpDir):
		fileName = self.getSettingsStructFilename()
		filePath = os.path.join(tmpDir,  fileName)
		out = open(filePath, 'w')
			#	CmdParamOutfileTypeInfo.usingNXSOutputCmdParam = False; out.write('#if !defined(IMPLEMENTING_GENERIC_OUTPUT_CMD_PARAM)\n')
		self.composeSettingsHeader(out, fileName)
			#	CmdParamOutfileTypeInfo.usingNXSOutputCmdParam = True; out.write('#else //if defined(IMPLEMENTING_GENERIC_OUTPUT_CMD_PARAM)\n')
			#	self.composeSettingsHeader(out, fileName)
			#	out.write('#endif //defined(IMPLEMENTING_GENERIC_OUTPUT_CMD_PARAM)\n')
		out.close()
	def writeSetupFile(self, tmpDir):
		fileName = self.getSetupFilename()
		filePath = os.path.join(tmpDir,  fileName)
		out = open(filePath, 'w')
			#	CmdParamOutfileTypeInfo.usingNXSOutputCmdParam = False; out.write('#if !defined(IMPLEMENTING_GENERIC_OUTPUT_CMD_PARAM)\n')
		self.composeSetupFile(out)
			#	CmdParamOutfileTypeInfo.usingNXSOutputCmdParam = True; out.write('#else //if defined(IMPLEMENTING_GENERIC_OUTPUT_CMD_PARAM)\n')
			#	self.composeSetupFile(out)
			#	out.write('#endif //defined(IMPLEMENTING_GENERIC_OUTPUT_CMD_PARAM)\n')
		out.close()
	def composeSettingsHeader(self, out, fileName):
		sentinelString = 'PHYC_' + fileName.translate(dotToUnderscoreTable).upper()
		paramTypeInstanceRequiredHeader = Set()
		for cmdParam in self.cmd_param:
			paramTypeInstanceRequiredHeader.update(cmdParam.getVarInstanceHeaders())
		className = self.getSettingStructName()
		#out = cStringIO.StringIO()
		#out.write('#if defined (OLD_XML_CMDS_PROCESSING)\n#error including wrong autogenerated code\n#endif //defined (OLD_XML_CMDS_PROCESSING)\n')
		out.write('/*########## %s ##########*/\n' % fileName)
		out.write('#if !defined (%s)\n#define %s\n' % (sentinelString, sentinelString) )
		for i in paramTypeInstanceRequiredHeader:
			out.write('#include "%s"\n'%i)
		out.write('class %s\n\t{\n\tpublic:\n\t\t' % className)
		for cmdParam in self.cmd_param:
			cmdParam.writeVarDeclaration(out)
			out.write('\n\t\t')
		out.write('\n\t\t%s()\n\t\t\t:' % className)
		first = True
		for cmdParam in self.cmd_param:
			if not first:
				out.write(',')
				out.write('\n\t\t\t')
			first = False
			cmdParam.writeVarDefinition(out)
		out.write('\n\t\t\t')
		out.write('{}\n')
		out.write('\t};\n')
		out.write('#endif // if !defined (%s)\n' % sentinelString )
		out.write('/*%%%%%%%%%%%% /%s %%%%%%%%%%%%*/\n' % fileName)
		

class CommandSet(SAXConstructible):
	def writeOutput(self, tmpDir):
		for command in self.command:
			command.writeOutput(tmpDir)

class CommandList(SAXConstructible):
	def writeOutput(self, tmpDir):
		for commandSet in self.command_set:
			commandSet.writeOutput(tmpDir)

if __name__ == '__main__':
	TextOnlyElement.automaticallyConvertToString = True
	SAXConstructible.initializeFromClassStatic = True
	writingNCL = True
	if len(sys.argv) > 1:	pathToPhycasRoot = sys.argv[1]
	else:	pathToPhycasRoot = os.path.join('..', '..')
	if len(sys.argv) > 2:	relativePathToDest = sys.argv[2]
	else:	relativePathToDest = os.path.join('tmp', 'phycas', 'command')
	destDir = os.path.join(pathToPhycasRoot,  relativePathToDest)
	xmlFilename = os.path.join(pathToPhycasRoot,  'gui', 'xml', 'all_commands.xml')
	print 'reading ', xmlFilename, ' and writing to ', destDir
	
	initStaticsFile = os.path.join(pathToPhycasRoot , 'python', 'phycas_sax_init_util.py')
	print ' about to execute:  ', initStaticsFile
	execfile(initStaticsFile)
   	cmdList = CommandList('phyc:command_list', {})
	if not cmdList.parseFile(xmlFilename):
		print xmlFilename, ' not parsed'
		sys.exit(1)
	cmdList.writeOutput(destDir)

