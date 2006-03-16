#! /usr/bin/python2.2
import sys, cStringIO, re, os
if os.environ.get('PHYCAS_ROOT') is not None:
	print 'appending'
	sys.path.append(os.path.join(os.environ['PHYCAS_ROOT'], 'python'))
from sax_constructible import SAXConstructible, ChildrenToHandle, TextOnlyElement, IgnoreContentElement
def uncapitalize(s):
	if len(s) > 1: return str(s[0]).lower() + s[1:]
	return s.lower()
class CommandRequireOrAvail(SAXConstructible):
		#this class seems like an unnecessary wrapper - we should make sure that we are going to meaningfully extend it or delete it
	def writeTestDefinition(self, out):
		self.test_info.writeTestDefinition(out)
	def getNCLTestInstanceName(self): 
		return self.test_info.getNCLTestInstanceName()
class CommandTestInfo(SAXConstructible):
	def getName(self):
		return self.comparison.name
	name = property(getName)
	def getNCLTestInstanceName(self): return 'test_' + self.name
	def supportedTest(self):
		return bool(self.comparison)
	def writeTestDefinition(self, out):
		'''currently only support comparisons'''
		if self.supportedTest():
			self.message.writeVarDefinition(out, self.name)
		if self.comparison:
			self.comparison.writeTestDefinition(out, self.name)
		return None
class CommandTest(SAXConstructible):
	xmlToNclOp = {'and' : ('kBoth', 'logical'), 
		'xor' : ('kXor', 'logical'), 
		'or' : ('kEither', 'logical'), 
		'if_fir_sec' : ('kIfFirstSecond', 'logical'), 
		'not_and' : ('kNotBoth', 'logical'),
		'not_xor': ('kNotOne', 'logical'), 
		#'not_or' : ('kBothOrNeither', 'logical'),   not implemented in ncl
		'greater_than' : ('kGreaterThan', 'comparison'), 
		'less_than' : ('kLessThan', 'comparison'), 
		'equals' : ('kEquals', 'comparison'), 
		'not_equal' : ('kNotEqual', 'comparison'), 
		'less_or_equal' : ('kLessOrEq', 'comparison'), 
		'greater_or_equal' : ('kGreaterOrEq', 'comparison'), }
	def getCPPType(self):	
		return 'UInt' # \\@ TEMP
	type = property(getCPPType)
	def toNCLOperator(op):
			#kNeither, and kFirstOnly are not handled in our xml schema
		return 'NxsTest::' + str(CommandTest.xmlToNclOp[op][0])
	toNCLOperator = staticmethod(toNCLOperator)
	def writeTestDefinition(self, out, name):
		self.left_operand.writeVarDefinition(out, 'left_' + name, self.type)
		self.right_operand.writeVarDefinition(out, 'right_' + name, self.type)
		operatorString = CommandTest.toNCLOperator(self.operator)
		varType = self.type + self.left_operand.nclTag + 'To' + self.right_operand.nclTag + 'Test'
		s = '\tNxsTestShPtr test_' + name + ' = NxsTestShPtr(new ' + varType + '(left_' + name + ', ' + operatorString + ', right_' + name + ', sms_' + name + '));\n'
		out.write(s)

class CommandTestOperand(SAXConstructible):
	def getNCLOperandType(self, type):
		s = self.labile and 'Callback' or self.constant_val and 'Value' or 'Var'
		return 'Nxs%s%sWrapper' % (type, s)
	def getNCLTag(self):
		return self.labile and 'Func' or self.constant_val and 'Const' or 'Var'
	nclTag = property(getNCLTag)
	def writeVarDefinition(self, out, name, type):
		pref = self.getNCLOperandType(type) + ' ' + name + '('
		if  self.labile:
			s = self.labile.phycas_impl.getBoostFuncDef()
		elif self.constant_val:
			s = self.constant_val
		else:
			raise ValueError, 'unimplemented'
		fullDef = '\t' + pref + s + ');\n'
		out.write(fullDef)
class CommandTestMessage(SAXConstructible):
	def writeVarDefinition(self, out, name):
		fm = self.failure and '"%s"' % self.failure or 'std::string()'
		sm = self.success and '"%s"' % self.success or 'std::string()'
		s = '\t' + 'StoredMsgSource sms_%s(%s, %s);\n' % (name, fm, sm)
		out.write(s)
	

class CmdParam(SAXConstructible):
	_typeInfoMatchPattern = re.compile(r'.+_type_info$')
	
	def setCommandParent(self, parent):
		self._commandParent = parent
	def getLabel(self): return self.type_info.decorateParamLabel(self.label)
	def getLabelOrVar(self): 
		s = self.getLabel()
		if len(s) == 0:
			return self.getCPPVarLabel()
		return s
	def getType(self): return self.type_info.getType()
	def childDone(self, obj, name):
		SAXConstructible.childDone(self, obj, name)
		if CmdParam._typeInfoMatchPattern.match(name):
			self.type_info = obj
	def endSelfElement(self, name): self.type_info.setBase(self)	
	def getVarInstanceHeaders(self): return self.type_info.getVarInstanceHeaders()
	def getCPPVarLabel(self):
		if not self.phycas_impl: return ''
		s = self.phycas_impl.getManipulatedVariableName()
		return uncapitalize(s)
	def writeVarDefinition(self, out): self.type_info.writeVarDefinition(out,  self.getCPPVarLabel())
	def writeVarDeclaration(self, out): self.type_info.writeVarDeclaration(out,  self.getCPPVarLabel())
	def getBoostRequiredHeaders(self): return ['<boost/shared_ptr.hpp>', '<boost/bind.hpp>']
	def getNCLCmdParamRequiredHeaders(self): return self.type_info.getNCLCmdParamInstanceHeaders() + self.getBoostRequiredHeaders()
	def writeInstantiateNCLParamObj(self, out):  self.type_info.writeInstantiateNCLParamObj(out)
	def getSetupInstanceName(self): return self.type_info.getSetupInstanceName()
	def getDescription(self): return self.description
	def appendMixed(self, obj):
		p = self.__dict__.get('_mixed', [])
		p.append(obj)
		self._mixed = p
	def followWithEquals(self): return self.type_info.followWithEquals()
	
class CmdParamTypeInfo(SAXConstructible):
	def setBase(self, b): self.base = b
	def getVarInstanceHeaders(self): return []
	def getCPPType(self): return 'C++Type'
	def getDefaultCPPValue(self):
		if self.__dict__.get('default'): return self.default
		return ''
	def writeVarDefinition(self, out, label): 
		initArg = self.getDefaultCPPValue()
		#if initArg == '""': initArg = ''
		out.write('%s(%s)' % (label, initArg))
	def writeVarDeclaration(self, out, label): out.write('%-16s\t%s;' % (self.getCPPType(), label))
	def getNCLCmdParamInstanceHeaders(self): return []
	def getSetupInstanceName(self): return 'v' + self.base.getLabelOrVar() + 'CmdOpt'
	def getStartNCLOptCtorArgs(self): 
		lab = '"' + self.base.getLabel() + '"'
		if lab == '""': lab = 'std::string()'
		return lab + ', &settingsStruct->' + self.base.getCPPVarLabel()
	def getMiddleNCLOptCtorArgs(self): return self.getDefaultAsString()
	def getDefaultAsString(self): return str(self.default)
	def getFinalNCLOptCtorArgs(self):
		return self.getPersistenceAndUser()
	def getPersistenceAndUser(self):
		p = str(bool(self.base.persistent)).lower()
		if self.base.user_level == 'basic_user': l = 'kCmdPermBasicUser'
		elif self.base.user_level == 'advanced_user': l = 'kCmdPermAdvancedUser'
		elif self.base.user_level == 'developer': l = 'kCmdPermDeveloper'
		else: raise TypeError, 'Unexpected user_level: %s' % self.base.user_level
		return p + ', ' + l
	def writeInstantiateNCLParamObj(self, out):
		nclOptType = self.getNCLOptType()
		setupInstanceName = self.getSetupInstanceName()
		nclOptCtorArg = self.getStartNCLOptCtorArgs() + ', ' +  self.getMiddleNCLOptCtorArgs() + ', ' + self.getFinalNCLOptCtorArgs()
		out.write('\t%s * %s = new %s(%s);\n' % (nclOptType, setupInstanceName, nclOptType, nclOptCtorArg))
	def decorateParamLabel(self, label): return label
	def followWithEquals(self): return False

class CmdParamBoolTypeInfo(CmdParamTypeInfo):
	def getCPPType(self): return 'bool'
	def getNCLCmdParamInstanceHeaders(self): return ['ncl/command/nxs_primitive_cmd_param.hpp']
	def getNCLOptType(self): return 'BoolCmdOption'
	
class CmdParamNumberTypeInfo(CmdParamTypeInfo):
	def getVariableValueProvider(self, gve):
		wrapperName = 'NamedNumberSource<%s>(' % self.getCPPType()
		if gve.label != None:	args = '"%s", '% gve.label + self.base._commandParent.getVarAddress(gve.label)
		else: args = '"", ' + gve.labile.phycas_impl.getBoostFuncDef()
		return 	wrapperName + args + ')'
	def getMinValAsString(self):
		if self.min_val.isVariable: return self.getVariableValueProvider(self.min_val)
		elif self.min_val.constant_val == '': return self.getMinOfCPPType()
		else: return self.min_val.constant_val
	def getMaxValAsString(self):
		if self.max_val.isVariable: return self.getVariableValueProvider(self.max_val)
		elif self.max_val.constant_val == '': return self.getMaxOfCPPType()
		else: return self.max_val.constant_val
	def getNCLCmdParamInstanceHeaders(self): return ['ncl/command/nxs_primitive_cmd_param.hpp']
	def getMiddleNCLOptCtorArgs(self):
		return self.getDefaultAsString() + ', '  + self.getMinValAsString() + ', ' + self.getMaxValAsString()
	
class CmdParamStringTextField(CmdParamTypeInfo):
	def getVarInstanceHeaders(self): return []
	def getCPPType(self): return 'std::string'
	def getDefaultCPPValue(self):
		s = str(self.default)
		if len(s) > 0:
			return '"%s"' % s
		return ''
	def getNCLCmdParamInstanceHeaders(self): return ['ncl/command/nxs_cmd_param.hpp']
	def getNCLOptType(self): return 'NxsStringCmdOption'
	def getDefaultAsString(self): 
		defStr = str(self.default)
		if len(defStr) > 0: return '"' + str(self.default) + '"'
		return 'std::string()'

class CmdParamNxsStringTextField(CmdParamStringTextField):
	def getMiddleNCLOptCtorArgs(self):
		return self.getDefaultAsString() + ', false'
	
class CmdParamRestrictedStringTypeInfo(CmdParamStringTextField):
	def getNCLCmdParamInstanceHeaders(self): return ['ncl/command/nxs_restricted_string_cmd_param.hpp']
	def getNCLOptType(self): return 'RestrictNameCmdOption'
	def getMiddleNCLOptCtorArgs(self): return self.getDefaultAsString() + ', '  + self.disallowed_values.getValueAsString()

class CmdParamNameTypeInfo(CmdParamRestrictedStringTypeInfo):
	def followWithEquals(self): return True

class CmdParamChoiceTypeInfo(CmdParamStringTextField):
	def getNCLCmdParamInstanceHeaders(self): return ['ncl/command/nxs_choice_cmd_param.hpp']
	def writeVarDeclaration(self, out, label): out.write('%-16s\t%s;\n\t\t%-16s\t%s;' % ('std::string', label, 'unsigned', self.phycas_impl.manipulated_var))
	def getNCLOptType(self): return 'NxsChoiceCmdOption'	
	def getStartNCLOptCtorArgs(self): 
		lab = '"' + self.base.getLabel() + '"'
		if lab == '""': lab = 'std::string()'
		return lab + ', &settingsStruct->' + self.phycas_impl.manipulated_var + ', &settingsStruct->' + self.base.getCPPVarLabel()
	def getMiddleNCLOptCtorArgs(self): 
		s = self.getDefaultAsString() + ', '  + self.choices.getValueAsString() + ', '
		if self.choices.isVariable:	return s + 'false'
		return s + 'true'

class CmdParamInfileTypeInfo(CmdParamStringTextField):
	def getVarInstanceHeaders(self): return ['ncl/misc/nxs_file_path.hpp']
	def getCPPType(self): return 'NxsInFilePath'
	def getNCLCmdParamInstanceHeaders(self): return ['ncl/command/nxs_file_cmd_param.hpp']
	def getNCLOptType(self): return 'NxsInFileCmdOption'
	def getMiddleNCLOptCtorArgs(self): return 'NxsInFilePath(' + self.getDefaultAsString() + ')'

class CmdParamOutfileTypeInfo(CmdParamStringTextField):
	def getFinalNCLOptCtorArgs(self): return CmdParamTypeInfo.getFinalNCLOptCtorArgs(self)
	def getDefaultCPPValue(self): return CmdParamStringTextField.getDefaultCPPValue(self)
	def getVarInstanceHeaders(self):  return  ['ncl/misc/nxs_file_path.hpp']
	def getCPPType(self): return 'NxsOutFilePath'
	def getNCLCmdParamInstanceHeaders(self): return ['ncl/command/nxs_file_cmd_param.hpp']
	def getNCLOptType(self): return  'NxsOutFileCmdOption' 
	def getMiddleNCLOptCtorArgs(self):
		if  self.base.getLabel() == 'File' or  self.base.getLabel() == 'Filename': b = 'false'
		else: b = 'true'
		return 'NxsOutFilePath(' + self.getDefaultAsString() + '), ' + b + ', ' + self.base._commandParent.getSetupCmdInstancePtrName() + '.get()'
	def getStartNCLOptCtorArgs(self):
		lab = '"' + self.base.getLabel() + '"'
		if lab == '""': lab = 'std::string()'
		return lab + ', "' + self.prefix + '", &settingsStruct->' + self.base.getCPPVarLabel()

class CmdParamIntegerTypeInfo(CmdParamNumberTypeInfo):
	def isUnsigned(self): return self.min_val.__dict__.get('constant_val') and self.min_val.constant_val >= 0
	def getCPPType(self): 
		if self.isUnsigned(): return 'unsigned'
		return 'int'
	def getNCLOptType(self):
		if self.isUnsigned(): return 'UIntCmdOption'
		return 'IntCmdOption'
	def getMinOfCPPType(self):
		if self.isUnsigned(): return '0'
		return '-INT_MAX'
	def getMaxOfCPPType(self):
		if self.isUnsigned(): return 'UINT_MAX'
		return 'INT_MAX'

class CmdParamDoubleTypeInfo(CmdParamNumberTypeInfo):
	def getCPPType(self): return 'double'
	def getNCLOptType(self): return 'DblCmdOption'
	def getMinOfCPPType(self): return '-DBL_MAX'
	def getMaxOfCPPType(self): return 'DBL_MAX'

class CmdParamSetTypeInfo(CmdParamTypeInfo):
	def getVarInstanceHeaders(self): return ['ncl/misc/nxs_index_set.hpp']
	def getCPPType(self): return 'NxsIndexSet'
	def getNCLCmdParamInstanceHeaders(self): return ['ncl/command/nxs_set_cmd_param.hpp']
	def getNCLOptType(self): return 'NxsIndexSetCmdOption'
	def getDefaultAsString(self): 
		defStr = str(self.default)
		if len(defStr) > 0: return '"' + str(self.default) + '"'
		return 'std::string()'
	def getDefaultCPPValue(self): return self.getDefaultAsString()
	def writeInstantiateNCLParamObj(self, out):
		nclOptType = self.getNCLOptType()
		setupInstanceName = self.getSetupInstanceName()
		nclOptCtorArg = self.getStartNCLOptCtorArgs() + ', ' +  self.getMiddleNCLOptCtorArgs() + ', ' + self.getFinalNCLOptCtorArgs()
		out.write('\t%s * %s = Instantiate%sSetOption(%s);\n' % (nclOptType, setupInstanceName, self.getSetPrefix(), nclOptCtorArg))
	def getMiddleNCLOptCtorArgs(self): return self.getDefaultAsString() + ', ' + self.phycas_impl.manager_ptr
	
class CmdParamTaxSetTypeInfo(CmdParamSetTypeInfo):
	def getSetPrefix(self):	return 'Tax'
	
class CmdParamTreeSetTypeInfo(CmdParamSetTypeInfo):
	def getSetPrefix(self):	return 'Tree'

class CmdParamCharSetTypeInfo(CmdParamSetTypeInfo):
	def getSetPrefix(self):	return 'Char'

class CmdParamDistributionTypeInfo(CmdParamStringTextField):
	def getVarInstanceHeaders(self): return ['phycas/rand/distribution_description.hpp']
	def getCPPType(self): return 'DistributionDescription'
	def getNCLCmdParamInstanceHeaders(self): return ['phycas/rand/distribution_command_param.hpp']
	def getNCLOptType(self): return 'DistributionCmdOption'
	def writeVarDefinition(self, out, label): out.write('%s()' % label)
	def getDistribClassCPPEnum(self):
		if self.distrib_class == 'Discrete': enVal =  'kDiscreteNonNegative'
		elif self.distrib_class == 'Continuous':  
			if self.range_constraint.constraint == 'NonNegative' : enVal = 'kNonNegative'
			elif self.range_constraint.constraint == 'SumToOne' : enVal = 'kSumToOne'
			elif self.range_constraint.constraint == 'Bounded': raise TypeError, 'Arbitrarily bounded distributions are not supported by NCL'
			elif self.range_constraint.constraint == 'Unbounded': enVal = 'kContinuous'
			elif self.range_constraint.constraint == 'Any' : enVal = 'kContinuous'
			else:
				raise TypeError , str('Unknown distrib constraint:  %s' % self.range_constraint.constraint)
		elif self.distrib_class == 'Any': enVal = 'kDiscreteOrContinuous'
		else:
			raise TypeError, 'Unknown distrib class:  %s' % self.distrib_class
		return 'DistributionCmdOption::' + enVal
	def getDistributionProviderFunctor(self):
		return 'boost::bind(&' + self.phycas_impl.getReceiverClass() + '::' + self.phycas_impl.getFuncName() + ', ' +  self.phycas_impl.getPtrName() + ', _1)'
	def getMiddleNCLOptCtorArgs(self): 
		return ', '.join([self.getDefaultAsString(),  self.getDistribClassCPPEnum(), self.num_variates, self.getDistributionProviderFunctor()])
		
class CmdParamMixedTypeInfo(CmdParamTypeInfo):
	def endSelfElement(self, name): self.setMixed()
	def setMixed(self):
		for cmdParam in self.cmd_param:
			cmdParam.appendMixed(cmdParam)
	def writeVarDeclaration(self, out, label):
		for cmdParam in self.cmd_param:
			cmdParam.writeVarDeclaration(out)
	def getNCLCmdParamInstanceHeaders(self):
		allh = ['ncl/command/nxs_mixed_cmd_param.hpp']
		allh.extend([i.getNCLCmdParamInstanceHeaders() for i in self.cmd_param])
		return allh
	def writeInstantiateNCLParamObj(self, out):
		for cmdParam in self.cmd_param:
			cmdParam.writeInstantiateNCLParamObj(out)
		CmdParamTypeInfo.writeInstantiateNCLParamObj(self, out)
	def getNCLOptType(self): return 'NxsMixedCmdOption'

class DistribRangeConstraint(SAXConstructible):
	pass	
class FileOutputElement(SAXConstructible):
	def getODDConstructorArgs(self):
		appArg = 'false'
		repArg = 'false'
		if self.append.upper() == 'TRUE': appArg = 'true'
		if self.replace.upper() == 'TRUE': repArg = 'true'
		return '"%s", %s, %s' % (self.path, appArg, repArg)
class OutputValue(SAXConstructible):
	def toCPPODDEnum(self, outstreamName):
		if outstreamName.upper() == 'OUTPUT': return 'NxsOutputDestinationDescription::kRedirectToOutputStream'
		if outstreamName.upper() == 'COMMENT': return 'NxsOutputDestinationDescription::kRedirectToCommentStream'
		if outstreamName.upper() == 'ERROR': return 'NxsOutputDestinationDescription::kRedirectToErrorStream'
		raise ValueError, 'Unknown Outstream Redirectio: %s\n' % outstreamName
	def getODDConstructorArgs(self):
		if self.suppress.upper() == "TRUE": return ''
		if len(self.redirect) > 0:
			if len(self.file) > 0:
				return self.toCPPODDEnum(self.redirect[0]) + ', ' + self.file[0].getODDConstructorArgs()
			return self.toCPPODDEnum(self.redirect[0])
		return self.file[0].getODDConstructorArgs()
class CmdParamOutputTypeInfo(CmdParamStringTextField):
	def getFinalNCLOptCtorArgs(self):
		return self.getPersistenceAndUser()
	def getDefaultCPPValue(self): return ''
	def getVarInstanceHeaders(self): return ['ncl/output/nxs_output_destination_description.hpp']
	def getCPPType(self): return 'NxsOutputDestinationDescription'
	def getNCLCmdParamInstanceHeaders(self): return ['ncl/command/nxs_file_cmd_param.hpp']
	def getNCLOptType(self): return 'NxsOutputCmdOption'
	def getMiddleNCLOptCtorArgs(self): return 'NxsOutputDestinationDescription(' + self.default.getODDConstructorArgs() + ')'
	def getStartNCLOptCtorArgs(self):
		lab = '"' + self.base.getLabel() + '"'
		return lab + ', &settingsStruct->' + self.base.getCPPVarLabel()



