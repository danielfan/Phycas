#! /usr/bin/python2.2
from sax_constructible import SAXConstructible, ChildrenToHandle, TextOnlyElement, IgnoreContentElement
import sys, cStringIO, re

def performTestInGUI(t):
	return 'allgui'.find(t.implemented_at) != -1

class CmdParam(SAXConstructible):
	_typeInfoMatchPattern = re.compile(r'.+_type_info$')
	def __init__(self, elName, parseContext):
		super(CmdParam, self).__init__(elName, parseContext)
		
	def setCommandParent(self, parent):
		self._commandParent = parent
		if self.type_info.__class__ == CmdParamMixedTypeInfo:
			for cmdParam in self.type_info.cmd_param:
				cmdParam._commandParent = parent
	
	def getComponentName(self, name = 'lbl_'):
		''' Name defaults to label, pass in other types '''
		if self.type_info.__class__ == CmdParamMixedTypeInfo:
			return name+'cmd'+self._commandParent.getCommandName()+'_'+self.type_info.cmd_param[0].getGUILabel()
		elif self._commandParent != None:
			return name+'cmd'+self._commandParent.getCommandName()+'_'+self.getGUILabel()
		else:
			return 'ERROR'
	
	def getType(self):
		return self.type_info.getType()
		
	def getNexusLabel(self):
		return str(self.label)
	
	def getGUILabel(self):
		return self.phycas_impl.manipulated_var

	def getGUIDesc(self):
		return self.gui_description

	def getValueName(self):
		return self.getComponentName(self.getType())+'_Value'
	
	def childDone(self, obj, name):
		SAXConstructible.childDone(self, obj, name)
		if CmdParam._typeInfoMatchPattern.match(name):
			self.type_info = obj
			
	def endSelfElement(self, name):
		self.type_info.setBase(self)
	
	def writeIntro(self, out):
		 self.type_info.writeIntro(out)

	def writeGetParams(self, out, type = 'reg'):
		 self.type_info.writeGetParams(out, type)

	def writeSetParams(self,out, type = 'reg'):
		 self.type_info.writeSetParams(out, type)

	def writeRequirements(self, out, par = None):
		for req in self.requirement:
			if performTestInGUI(req):
				self.writeTest(out, req.test_info, "R", par);
		if(self.type_info.__class__==CmdParamMixedTypeInfo):
			for cmdParam in self.type_info.cmd_param:
				cmdParam.writeRequirements(out, self)

	def writeAvailability(self, out, par = None):
		for avail in self.availability:
			if performTestInGUI(avail):
				self.writeTest(out, avail.test_info, "A");
		if(self.type_info.__class__ == CmdParamMixedTypeInfo):
			for cmdParam in self.type_info.cmd_param:
				cmdParam.writeAvailability(out, self)

	def writeChoice(self, out):
		self.type_info.writeChoice(out)

	def writeActionGetValues(self, out, outer = None):
		self.type_info.writeActionGetValues(out, outer)

	def writeActionValidation(self, out, outer = None):
		self.type_info.writeActionValidation(out, outer)

	def writeActionCmdString(self, out):
		self.type_info.writeActionCmdString(out)
	
	def writeDefActionCmdString(self, out):
		if len(self.getNexusLabel()) == 0:
			out.write('\t\tcommandString.append(" "+%s+" ");\n' % self.getValueName())
		else:
			out.write('\t\tcommandString.append(" %s = \'"+ %s +"\'");\n' % (self.getNexusLabel(), self.getValueName()))
	
	def writeGetErrMsgs(self, out, outer = None):
		self.type_info.writeGetErrMsgs(out,outer)

	def writeItemStateChanged(self, out):
		self.type_info.writeItemStateChanged(out)

	def writeDisallowedValues(self, out):
		self.type_info.writeDisallowedValues(out)

	def writeActionFileChooser(self, out):
		self.type_info.writeActionFileChooser(out)

	def writeActionDistribChooser(self, out):
		self.type_info.writeActionDistribChooser(out)

	def writeActionComboListSelect(self, out):
		self.type_info.writeActionComboListSelect(out)


	# VKJ 12/7/04 - The following have not been tested:  callback, mixed command, compound tests
	def writeTest(self, out, aTestInfo, type, parentCmdParam=None):
		if aTestInfo.callback != None:
			return
		aTest = ""
		predicate = 'false'
		if aTestInfo.predicate != None:
			predicate = 'true'
			aTest = aTestInfo.predicate
			value = ""
			valueGuiLabel = str(self.getRawChars())
			if aTestInfo.predicate.label == self.label or (aTestInfo.predicate.constant_val != None and aTestInfo.predicate.constant_val == 'this'):
				#if self.__class__ == TypeInfoMixed:
				#	for mixedCmdParam in self.cmd_param:
				#		if mixedCmdParam.id == aTest.left_operand.command_param.sub_cmd_opt_id:
				#			leftValueName = mixedCmdParam.getValueName()
				#			leftDataType = mixedCmdParam.getDataType()
				#else:
				if self.type_info.__class__ == CmdParamBoolTypeInfo:
					value = self.type_info.getSelectedName()
				else:
					value = self.getValueName()
			else:
				value = aTestInfo.predicate.constant_val;
		elif aTestInfo.logicalOp != None:
			aTest = aTestInfo.logicalOp
		elif aTestInfo.comparison != None:
			aTest = aTestInfo.comparison
		leftLabel = ""
		rightLabel = ""
		if predicate == 'false':
			operator = ""
			if aTest.operator == 'less_than':
				operator = '<'
			elif aTest.operator == 'less_or_equal':
				operator = '<='
			elif aTest.operator == 'greater_than':
				operator = '>'
			elif aTest.operator == 'greater_or_equal':
				operator = '>='
			elif aTest.operator == 'equals':
				operator = '=='
			elif aTest.operator == 'not_equal':
				operator = '!='
			elif aTest.operator == 'and':
				operator = '&&'
			elif aTest.operator == 'or':
				operator = '||'
			else:
				operator = aTest.operator
			if aTest.left_operand.label != None:
				leftLabel = aTest.left_operand.label
			elif aTest.left_operand.constant_val != None:
				leftLabel = aTest.left_operand.constant_val
			if aTest.right_operand.label != None:
				rightLabel = aTest.right_operand.label
			elif aTest.right_operand.constant_val != None:
				rightLabel = aTest.right_operand.constant_val
			leftValueName = ""
			rightValueName = ""
			leftDataType = ""
			rightDataType = ""
			if leftLabel == self.label or leftLabel == 'this':
				if self.type_info.__class__ == CmdParamBoolTypeInfo:
					leftValueName = self.type_info.getSelectedName()
				else:
					leftValueName = self.getValueName()
				leftDataType = self.type_info.getDataType()
			else:
				leftValueName = leftLabel;
			if rightLabel == self.label or rightLabel == 'this':
				if self.type_info.__class__ == CmdParamBoolTypeInfo:
					rightValueName = self.type_info.getSelectedName()
				else:
					rightValueName = self.getValueName()
				rightDataType = self.type_info.getDataType()
			else:
				rightValueName = rightLabel;
		valueGuiDesc = ""
		for cmdParamGroup in self._commandParent.cmd_param_group:
			for cmdParam in cmdParamGroup.cmd_param:
				if aTestInfo.predicate != None and (aTestInfo.predicate.label == cmdParam.label or (aTestInfo.predicate.constant_val != None and aTestInfo.predicate.constant_val == 'this')):
					if cmdParam.type_info.__class__ == CmdParamBoolTypeInfo:
						value = cmdParam.type_info.getSelectedName()
					else:
						value = cmdParam.getValueName()
					valueGuiDesc = cmdParam.gui_description
				elif aTestInfo.predicate != None and aTestInfo.predicate.constant_val != None:
					value = aTestInfo.predicate.constant_val;
				elif leftLabel == cmdParam.label:
					if cmdParam.type_info.__class__ == CmdParamBoolTypeInfo:
						leftValueName = cmdParam.type_info.getSelectedName()
					else:
						leftValueName = cmdParam.getValueName()
					leftDataType = cmdParam.type_info.getDataType()
				elif rightLabel == cmdParam.label:
					if cmdParam.type_info.__class__ == CmdParamBoolTypeInfo:
						rightValueName = cmdParam.type_info.getSelectedName()
					else:
						rightValueName = cmdParam.getValueName()
					rightDataType = cmdParam.type_info.getDataType()
				elif cmdParam.type_info.__class__ == CmdParamMixedTypeInfo:
					for mixedCmdParam in cmdParam.type_info.cmd_param:
						if leftLabel == mixedCmdParam.label:
							if mixedCmdParam.type_info.__class__ == CmdParamBoolTypeInfo:
								leftValueName = mixedCmdParam.type_info.getSelectedName()
							else:
								leftValueName = mixedCmdParam.getValueName()
							leftDataType = mixedCmdParam.type_info.getDataType()
						elif rightLabel == mixedCmdParam.label:
							if mixedCmdParam.type_info.__class__ == CmdParamBoolTypeInfo:
								rightValueName = mixedCmdParam.type_info.getSelectedName()
							else:
								rightValueName = mixedCmdParam.getValueName()
							rightDataType = mixedCmdParam.type_info.getDataType()
		if predicate == 'true':
			if aTestInfo.predicate.evaluate_to == 'true':
				out.write('\t\tif (%s) {\n' % value)
			elif aTestInfo.predicate.evaluate_to == 'false':
				out.write('\t\tif (!%s) {\n' % value)
			out.write('\t\t\tmain.%s.setEnabled(true);\n' % self.getComponentName(self.type_info.getType()))
			out.write('\t\t\tmain.%s.setEnabled(true);\n' % self.getComponentName())
			out.write('\t\t}\n')
			out.write('\t\telse {\n')
			out.write('\t\t\tmain.%s.setEnabled(false);\n' % self.getComponentName(self.type_info.getType()))
			out.write('\t\t\tmain.%s.setEnabled(false);\n' % self.getComponentName())
			out.write('\n\t\t\t')
			out.write(r'main.taMsgs.append("\n%s only available if %s is %s");' % (self.gui_description,valueGuiDesc,aTestInfo.predicate.evaluate_to))
			out.write('\n\t\t\tmain.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());\n')
			out.write('\t\t}\n')
		else:
			negativeTest = 1 #check for negative
			if aTestInfo.comparison != None:
				#out.write('\t\tif (')
				#leftDigit = 'false'
				#if(aTestInfo.comparison.left_operand.constant_val==None):
				#	if(leftValueName.isdigit()):
				#		leftDigit = 'true'
				#	else:
				#		out.write('!%s.equals("") ' % (leftValueName))
				#	if(rightValueName.isdigit() or aTestInfo.comparison.right_operand.constant_val!=None):
				#		pass
				#	else:
				#		out.write('&& !%s.equals("")) {\n' % (rightValueName))
				#else:	
				#	if(rightValueName.isdigit() or aTestInfo.comparison.right_operand.constant_val!=None):
				#		pass
				#	else:
				#		out.write('!%s.equals("")) {\n' % (rightValueName))
				#if((leftDigit=="false" and rightValueName.isdigit()) or aTestInfo.comparison.right_operand.constant_val!=None):
				#	out.write(') {\n')
				if(leftDataType=='String' and operator=='=='):
					out.write('\t\tif (!%s.equals(' %  leftValueName)
				elif(leftValueName.isdigit() or leftDataType=='String'):
					out.write('\t\tif (!(%s ' %  leftValueName)
				elif(aTestInfo.comparison.left_operand.constant_val!=None):
					out.write('\t\tif (!"%s".equals(' %  leftValueName)
				else:
					out.write('\t\tif (!(new %s(%s).%sValue() ' % \
						(leftDataType, leftValueName, leftDataType.lower()))
				if(leftDataType!='String' and operator!='=='):
					out.write('%s ' % operator)
				if(rightValueName.isdigit() or rightDataType=='String'):
					out.write('%s)) {\n' %  rightValueName)
				elif(aTestInfo.comparison.right_operand.constant_val!=None):
					out.write('"%s")) {\n' %  rightValueName)
				else:
					out.write('new %s(%s).%sValue())) {\n' % \
						(rightDataType, rightValueName, rightDataType.lower()))
				if(type=="R"):
					out.write('\t\t\tmain.%s.setText("");\n' % self.getComponentName(self.type_info.getType()))
			elif aTestInfo.logicalOp != None:
				negativeTest = 0; #check for positive
				if(operator=="&&"):
					out.write('\t\tif (%s && %s) {\n' % (leftValueName, rightValueName))
				elif(operator=="not_and"): 
					out.write('\t\tif (!(%s && %s)) {\n' % (leftValueName, rightValueName))
				elif(operator=="||"):
					out.write('\t\tif (%s || %s) {\n' % (leftValueName, rightValueName))
				elif(operator=="not_or"):
					out.write('\t\tif (!(%s || %s)) {\n' % (leftValueName, rightValueName))
				elif(operator=="xor"): 
					out.write('\t\tif ((%s && !%s) || (!%s && %s)) {\n' % (leftValueName, rightValueName, leftValueName, rightValueName))
				elif(operator=="not_xor"): 
					out.write('\t\tif (!((%s && !%s) || (!%s && %s))) {\n' % (leftValueName, rightValueName, leftValueName, rightValueName))
				elif(operator=="if_fir_sec"): 
					out.write('\t\tif ((%s && %s) || (!%s)) {\n' % (leftValueName, rightValueName, leftValueName))
			cnt1 = ""
			if(parentCmdParam!=None):
				cn = parentCmdParam.type_info.cmd_param[0].getComponentName() 
			else:
				cn = self.getComponentName()
			if(parentCmdParam!=None):
				cnt = parentCmdParam.type_info.cmd_param[0].getComponentName(parentCmdParam.type_info.cmd_param[0].getType()) 
			else:
				cnt = self.getComponentName(self.type_info.getType())
				if(self.type_info.__class__==CmdParamMixedTypeInfo):
					cnt1 = self.type_info.cmd_param[1].getComponentName(self.type_info.cmd_param[1].getType()) 
			if(type=="R"):
				if(not negativeTest):
					out.write('\t\t\tmain.%s.setForeground(java.awt.Color.BLACK);\n' % cn)
					out.write('\t\t}\n')
					out.write('\t\telse {\n')
					out.write('\t\t\tmain.%s.setForeground(java.awt.Color.RED);\n' % cn)
					out.write('\t\t\t')
					out.write(r'main.taMsgs.append("\n%s");' % aTestInfo.message.failure)
					out.write('\n\t\t\tmain.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());\n')
					out.write('\t\t\tsuccess = false;\n')
					out.write('\t\t}\n')
				else:
					out.write('\t\t\tmain.%s.setForeground(java.awt.Color.RED);\n' % cn)
					out.write('\t\t\t')
					out.write(r'main.taMsgs.append("\n%s");' % aTestInfo.message.failure)
					out.write('\n\t\t\tmain.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());\n')
					out.write('\t\t\tsuccess = false;\n')
					out.write('\t\t}\n')
					out.write('\t\telse {\n')
					out.write('\t\t\tmain.%s.setForeground(java.awt.Color.BLACK);\n' % cn)
					out.write('\t\t}\n')
			elif(type=="A"):
				if(not negativeTest):
					out.write('\t\t\tmain.%s.setEnabled(true);\n' % cnt)
					if(cnt1!=""):
						out.write('\t\t\tmain.%s.setEnabled(true);\n' % cnt1)
					out.write('\t\t\tmain.%s.setEnabled(true);\n' % cn)
					out.write('\t\t}\n')
					out.write('\t\telse {\n')
					out.write('\t\t\tmain.%s.setEnabled(false);\n' % cnt)
					if(cnt1!=""):
						out.write('\t\t\tmain.%s.setEnabled(false);\n' % cnt1)
					out.write('\t\t\tmain.%s.setEnabled(false);\n' % cn)
					out.write('\n\t\t\t')
					out.write(r'main.taMsgs.append("\n%s");' % (aTestInfo.message.failure))
					out.write('\n\t\t\tmain.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());\n')
					out.write('\t\t}\n')
				else:
					out.write('\t\t\tmain.%s.setEnabled(false);\n' % cnt)
					if(cnt1!=""):
						out.write('\t\t\tmain.%s.setEnabled(false);\n' % cnt1)
					out.write('\t\t\tmain.%s.setEnabled(false);\n' % cn)
					out.write('\n\t\t\t')
					out.write(r'main.taMsgs.append("\n%s");' % (aTestInfo.message.failure))
					out.write('\n\t\t\tmain.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());\n')
					out.write('\t\t}\n')
					out.write('\t\telse {\n')
					out.write('\t\t\tmain.%s.setEnabled(true);\n' % cnt)
					if(cnt1!=""):
						out.write('\t\t\tmain.%s.setEnabled(true);\n' % cnt1)
					out.write('\t\t\tmain.%s.setEnabled(true);\n' % cn)
					out.write('\t\t}\n')
			#if aTestInfo.comparison != None:
			#	out.write('\t\t}\n')
		
class CmdParamTypeInfo(SAXConstructible):
	def setBase(self, b):
		self.base = b

	def writeActionCmdString(self, out):
		self.base.writeDefActionCmdString(out)

#	def writeIntro(self, out): pass

	def writeGetErrMsgs(self, out, outer = None): pass

	def writeItemStateChanged(self, out): pass

	def writeDisallowedValues(self, out): pass

	def writeActionFileChooser(self, out): pass

	def writeActionDistribChooser(self, out): pass

	def writeActionComboListSelect(self, out): pass

	def writeGetParams(self, out, type = 'reg'): pass

	def writeSetParams(self,out, type = 'reg'): pass

class CmdParamBoolTypeInfo(CmdParamTypeInfo):
	def getType(self):
		return 'cb_Bool_'

	def getDataType(self):
		return 'Boolean'

	def getSelectedName(self):
		return self.base.getComponentName(self.getType())+'_Selected'

	def writeIntro(self, out):
		vn = self.base.getValueName()
		sn = self.getSelectedName()
		selected = 'false'
		if not self.default.isVariable and len(self.default.constant_val) > 0:
			selected = self.default
		out.write('\tprivate String %s = "";\n' % vn)
		out.write('\tprivate boolean %s = %s;\n' % (sn, selected))

	def writeActionGetValues(self, out, outer = None):
		pass

	def writeActionValidation(self, out, outer = None):
		pass

	def writeActionCmdString(self, out):
		out.write('\t\t\tcommandString.append(" %s = ");\n' % self.base.getNexusLabel())
		out.write('\t\t\tif(%s)\n' % self.getSelectedName())
		out.write('\t\t\t\tcommandString.append("true");\n')
		out.write('\t\t\telse\n')
		out.write('\t\t\t\tcommandString.append("false");\n')

	def writeItemStateChanged(self, out):
		out.write('\t\tif( (itemEvent.getStateChange() == java.awt.event.ItemEvent.SELECTED) ')
		out.write('&& (((JCheckBox)itemEvent.getItem()).getName().equals("%s")) )\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\t%s = true;\n' % self.getSelectedName())
		out.write('\t\tif( (itemEvent.getStateChange() == java.awt.event.ItemEvent.DESELECTED) ')
		out.write('&& (((JCheckBox)itemEvent.getItem()).getName().equals("%s")) )\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\t%s = false;\n' % self.getSelectedName())

	def writeGetParams(self, out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if self.base.label == "" or self.base.label == None:
			out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s") {\n' % (self.base.label,self.base.placement))
		else:
			out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\t\tBoolTypeInfo typeInfo = allParams[i].getBoolTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefBoolValueElement value = typeInfo.getDefault();\n')
                out.write('\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tmain.%s.setSelected(value.getConstantVal());\n' % cn)
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

	def writeSetParams(self,out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\t\tBoolTypeInfo typeInfo = allParams[i].getBoolTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefBoolValueElement value = typeInfo.getDefault();\n')
                out.write('\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tvalue.setConstantVal(%s);\n' % self.getSelectedName())
                out.write('\t\t\t\ttypeInfo.setDefault(value);\n') 
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

class CmdParamNumberTypeInfo(CmdParamTypeInfo):
	def writeGetErrMsgs(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		gd = self.base.getGUIDesc()
	        out.write('\t\tif(tfName.equals("%s"))\n\t\t\t' % cn)
	        if outer == None:
			out.write(r'ret =  "\n%s must be a number between "' % gd);
		else:
			out.write(r'ret =  "\n%s must be a number between "' % outer.base.getGUIDesc());
		out.write('+ %s_min +" and " + %s_max;' % (self.base.getGUILabel(), self.base.getGUILabel()));
		out.write('\n');

	def writeActionGetValues(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		vn = self.base.getValueName()
		out.write('\t\t%s = main.%s.getText();\n' % (vn, cn))

	def writeActionValidation(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		vn = self.base.getValueName()
		#out.write('\t\t%s = main.%s.getText();\n' % (vn, cn))
		out.write('\t\tif(!validate%s(%s, %s_min, ' % (self.getDataType(), vn, self.base.getGUILabel()))
		tfValueMin = ''
		tfValueMax = ''
		if self.min_val != None and self.min_val.label != None:
			for cmdParam in self.base._commandParent.cmd_param:
				if(cmdParam.label==self.min_val.label):
					tfValueMin = cmdParam.getValueName()
			out.write('"%s", %s' % (self.min_val.label,tfValueMin))
		else:
			out.write('null, null, ')
		out.write('%s_max' % self.base.getGUILabel())
		if self.max_val != None and self.max_val.label != None:
			tfValue = ''
			for cmdParamGroup in self.base._commandParent.cmd_param_group:
				for cmdParam in cmdParamGroup.cmd_param:
					if(cmdParam.label==self.max_val.label):
						tfValueMax = cmdParam.getValueName()
			out.write(', "%s", %s' % (self.max_val.label,tfValueMax))
		else:
			out.write(', null, null')
		out.write(')) {\n ')
		#out.write('\t\tif(!validate%s(%s, %s_min, %s_max)) {\n' % (self.getDataType(), vn, self.base.getGUILabel(), self.base.getGUILabel()))
		if outer == None:
			if self.min_val != None and self.min_val.label != None:
				out.write('\t\t\t%s_min = new Integer(%s).intValue();\n' % (self.base.getGUILabel(), tfValueMin))
			if self.max_val != None and self.max_val.label != None:
				out.write('\t\t\t%s_max = new Integer(%s).intValue();\n' % (self.base.getGUILabel(), tfValueMax))
			out.write('\t\t\tsetValidationError(main.%s, main.%s);\n' % (cn, self.base.getComponentName()))
			out.write('\t\t\tsuccess = false;\n')
		else:
			outer.cmd_param[1].writeActionGetValues(out, outer)
			outer.cmd_param[1].writeActionValidation(out, outer)
		out.write('\t\t}\n')
		out.write('\t\telse\n')
		out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())
		
	def getType(self):
		return 'tf_Number_'


class CmdParamStringTextField(CmdParamTypeInfo):
	def getType(self):
		return 'tf_String_'

	def getDataType(self):
		return 'String'
		
	def writeIntro(self, out):
		vn = self.base.getValueName()
		out.write('\tprivate String %s = "";\n' % vn)
		
	def writeGetErrMsgs(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		gd = self.base.getGUIDesc()
	        out.write('\t\tif(tfName.equals("%s"))\n\t\t\t' % cn)
	        if outer == None:
			out.write(r'ret =  "\n%s must be text";' % gd);
		else:
			out.write(r'ret =  "\n%s must be text";' % outer.base.getGUIDesc());
		out.write('\n');

	def writeActionGetValues(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		vn = self.base.getValueName()
		out.write('\t\t%s = main.%s.getText();\n' % (vn, cn))

	def writeActionValidation(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		vn = self.base.getValueName()
		#out.write('\t\t%s = main.%s.getText();\n' % (vn, cn))
		out.write('\t\tif(!validate%s(%s)) {\n' % (self.getDataType(), vn))
		if outer == None:
			out.write('\t\t\tsetValidationError(main.%s, main.%s);\n' % (cn, self.base.getComponentName()))
			out.write('\t\t\tsuccess = false;\n')
		else:
			outer.cmd_param[1].writeActionGetValues(out, outer)
			outer.cmd_param[1].writeActionValidation(out, outer)
		out.write('\t\t}\n')
		out.write('\t\telse\n')
		out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())
		
class CmdParamDistributionTypeInfo(CmdParamStringTextField):
	def writeActionCmdString(self, out):
		if len(self.base.getNexusLabel()) == 0:
			out.write('\t\tcommandString.append(" "+%s+" ");\n' % self.base.getValueName())
		else:
			out.write('\t\tcommandString.append(" %s = \"+ %s +"\");\n' % (self.base.getNexusLabel(), self.base.getValueName()))

	def writeActionDistribChooser(self, out):
		out.write('\tif(e.getActionCommand().equals("chooseAction_%s")) {\n' % self.base.getComponentName(''))
                out.write('\t\t\tphycasGUI.swixml.PhycasDistributionChooser dialog = new phycasGUI.swixml.PhycasDistributionChooser(main.frame, allParams);\n')
                out.write('\t\t\tdialog.pack();\n')
                out.write('\t\t\tdialog.show();\n')
                out.write('\t\t\tmain.%s.setText(dialog.getChoice());\n' % self.base.getComponentName(self.base.getType()))
		out.write('\t\t}\n')

	def writeGetErrMsgs(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		gd = self.base.getGUIDesc()
	        out.write('\t\tif(tfName.equals("%s"))\n\t\t\t' % cn)
	        if outer == None:
			out.write(r'ret =  "\n%s must be text in the form [distribution name]([parameters]), i.e. Uniform(0,1).";' % gd);
		else:
			out.write(r'ret =  "\n%s must be text in the form [distribution name]([parameters]), i.e. Uniform(0,1).";' % outer.base.getGUIDesc());
		out.write('\n');

	def writeActionValidation(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		vn = self.base.getValueName()
		out.write('\t\tif(!validateString(%s)) {\n' % vn)
		if outer == None:
			out.write('\t\t\tsetValidationError(main.%s, main.%s);\n' % (cn, self.base.getComponentName()))
			out.write('\t\t\tsuccess = false;\n')
		else:
			outer.cmd_param[1].writeActionGetValues(out, outer)
			outer.cmd_param[1].writeActionValidation(out, outer)
		out.write('\t\t}\n')
		out.write('\t\telse\n')
		out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())

class DistribRangeConstraint(SAXConstructible):
	pass

class CmdParamNxsStringTextField(CmdParamStringTextField):
	
	def writeGetParams(self, out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if type == 'reg':
			if self.base.label == "" or self.base.label == None:
				out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s") {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tStringTypeInfo typeInfo = allParams[i].getStringTypeInfo();\n')
		else:
			if self.base.label == "" or self.base.label == None:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s") && allSubParams[j].getPlacement().toString().equals("%s") {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tStringTypeInfo typeInfo = allSubParams[j].getStringTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefStringValueElement value = typeInfo.getDefault();\n')
                out.write('\t\t\tif(value!=null && value.getConstantVal()!=null) {\n')
                out.write('\t\t\t\tmain.%s.setText(value.getConstantVal());\n' % cn)
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

	def writeSetParams(self,out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if type == 'reg':
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tStringTypeInfo typeInfo = allParams[i].getStringTypeInfo();\n')
		else:
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s") && allSubParams[j].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tStringTypeInfo typeInfo = allSubParams[j].getStringTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefStringValueElement value = typeInfo.getDefault();\n')
                out.write('\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tvalue.setConstantVal(%s);\n' % self.base.getValueName())
                out.write('\t\t\t\ttypeInfo.setDefault(value);\n')
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

class CmdParamRestrictedStringTypeInfo(CmdParamStringTextField):

	def getType(self):
		return 'tf_RestrictedString_'
		
	def writeIntro(self, out):
		vn0 = self.base.getValueName()
		out.write('\tprivate String %s = "";\n' % vn0)
		out.write('\tprivate String[] %s_Disallowed = null;\n' % self.base.getComponentName(self.getType()))

	def writeDisallowedValues(self, out):
		cn = self.base.getComponentName(self.getType())
		out.write('\t/**\n')
		out.write('\t* Initializes any disallowed values for a restricted string type parameter\n')
		out.write('\t*/\n')  
		out.write('\tpublic void setInitialDisallowedValues() {\n')
		listOfValues = []
		for i in range(0, len(self.disallowed_values.constant_val)):
			listOfValues.append(self.disallowed_values.constant_val[i])
		out.write('\t\t%s_Disallowed = new String[%s];\n' % (cn,len(listOfValues)))
		for i in range(0, len(listOfValues)):
			out.write('\t\t%s_Disallowed[%s] = "%s";\n' % (cn,i,listOfValues[i]))
		out.write('\t}\n')
		out.write('\n')

		out.write('\t/**\n')
		out.write('\t* Sets the disallowed values for a restricted string type parameter\n')
		out.write('\t* @param disallowed disallowed values\n')
		out.write('\t*/\n')  
		out.write('\tpublic void setDisallowedValues(String[] disallowed) {\n')
		out.write('\t\tif(%s_Disallowed == null) {\n' % cn)
		out.write('\t\t\t%s_Disallowed = new String[disallowed.length];\n' % cn)
		out.write('\t\t\tfor(int i=0;i<disallowed.length;i++)\n')
		out.write('\t\t\t\t%s_Disallowed[i] = disallowed[i];\n' % cn)
		out.write('\t\t}\n')
		out.write('\t\telse {\n')
		out.write('\t\t\t%s_Disallowed = mergeArrays(%s_Disallowed,disallowed);\n' % (cn,cn))
		out.write('\t\t}\n')
		out.write('\t}\n')
		out.write('\n')
		
		out.write('\tprivate String[] mergeArrays(String[] pa, String[] pb) {\n')
		out.write('\t\tString[] arr = new String[pa.length+pb.length];\n')
		out.write('\t\tfor (int x=0; x < pa.length; x++) {\n')
		out.write('\t\t\tarr[x] = pa[x];\n')
		out.write('\t\t}\n')
		out.write('\t\tfor (int x=0; x < pb.length; x++) {\n')
		out.write('\t\t\tarr[x+pa.length] = pb[x];\n')
		out.write('\t\t}\n')
		out.write('\t\treturn arr;\n')
		out.write('\t}\n')

	def writeActionGetValues(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		vn = self.base.getValueName()
		out.write('\t\t%s = main.%s.getText();\n' % (vn, cn))

	def writeActionValidation(self, out, outer = None):
		cn0 = self.base.getComponentName(self.getType())
		vn0 = self.base.getValueName()
		#out.write('\t\t%s = main.%s.getText();\n' % (vn0, cn0))
		out.write('\t\tif(!validate%s(%s, %s_Disallowed)) {\n' % (self.getDataType(), vn0, cn0))
		out.write('\t\t\tsetValidationError(main.%s, main.%s);\n' % (cn0, self.base.getComponentName()))
		out.write('\t\t\tsuccess = false;\n')
		out.write('\t\t}\n')
		out.write('\t\telse\n')
		out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())

	def writeActionCmdString(self, out):
		out.write('\t\t\tcommandString.append(" "+%s+" = ");\n' % self.base.getValueName())

	def writeGetParams(self, out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\t\tRestrictedStringTypeInfo typeInfo = allParams[i].getRestrictedStringTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tStringValueList values = typeInfo.getDisallowedValues();\n')
                out.write('\t\t\tsetDisallowedValues(values.getConstantValArray());\n')
                out.write('\t\t}\n')

	#def writeSetParams(self,out, type = 'reg'):
	#	pass
	
	def writeGetErrMsgs(self, out, outer = None):
		cn0 = self.base.getComponentName(self.getType())
		gd = self.base.getGUIDesc()
	        out.write('\t\tif(tfName.equals("%s"))\n\t\t\t' % cn0)
		out.write(r'ret =  "\nThe value you have chosen for %s is invalid.";' % gd);
		out.write('\n');

class CmdParamNameTypeInfo(CmdParamRestrictedStringTypeInfo):
	def writeGetParams(self, out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\t\tRestrictedStringTypeInfo typeInfo = allParams[i].getNameTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tStringValueList values = typeInfo.getDisallowedValues();\n')
                out.write('\t\t\tsetDisallowedValues(values.getConstantValArray());\n')
                out.write('\t\t}\n')

	
class CmdParamChoiceTypeInfo(CmdParamTypeInfo):
	def endSelfElement(self, name):
		if len(self.choices.constant_val) > 4 or self.choices.labile != None: 
			self.__class__ = CmdParamComboBox
		else: 
			self.__class__ = CmdParamRadioButton

class CmdParamFileTypeInfo(CmdParamStringTextField):
	def getType(self):
		return 'tf_File_'
	
class CmdParamInfileTypeInfo(CmdParamFileTypeInfo):
	def writeIntro(self, out):
		vn = self.base.getValueName()
		out.write('\tprivate String %s = "";\n' % vn)

	def writeActionFileChooser(self, out):
		out.write('\tif(e.getActionCommand().equals("chooseAction_%s")) {\n' % self.base.getComponentName(''))
		out.write('\t\t\tint retVal = main.fc.showOpenDialog(main.%s);\n' % self.base.getComponentName(self.getType()))
                out.write('\t\t\tif (retVal == JFileChooser.APPROVE_OPTION) {\n')
		out.write('\t\t\t\tFile file = main.fc.getSelectedFile();\n')
		out.write('\t\t\t\tmain.%s.setText(file.getAbsolutePath());\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\t\tmain.fc.setSelectedFile(new File(""));\n')
		out.write('\t\t\t}\n')
		out.write('\t\t}\n')

	def writeActionCmdString(self, out):
		out.write('\t\t\tcommandString.append(" %s = \'"+ %s +"\'");\n' % (self.base.getNexusLabel(), self.base.getValueName()))

	def writeGetParams(self, out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\t\tStringTypeInfo typeInfo = allParams[i].getInfileTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefStringValueElement value = typeInfo.getDefault();\n')
                out.write('\t\t\tif(value!=null && value.getConstantVal()!=null) {\n')
                out.write('\t\t\t\tmain.%s.setText(value.getConstantVal());\n' % cn)
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

	def writeSetParams(self,out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\t\tStringTypeInfo typeInfo = allParams[i].getInfileTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefStringValueElement value = typeInfo.getDefault();\n')
                out.write('\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tvalue.setConstantVal(%s);\n' % self.base.getValueName())
                out.write('\t\t\t\ttypeInfo.setDefault(value);\n')
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

class CmdParamOutfileTypeInfo(CmdParamFileTypeInfo):
	def writeIntro(self, out):
		vn = self.base.getValueName()
		out.write('\tprivate String %s = "";\n' % vn)

	def writeActionFileChooser(self, out):
		out.write('\tif(e.getActionCommand().equals("chooseAction_%s")) {\n' % self.base.getComponentName(''))
		out.write('\t\t\tint retVal = main.fc.showSaveDialog(main.%s);\n' % self.base.getComponentName(self.getType()))
                out.write('\t\t\tif (retVal == JFileChooser.APPROVE_OPTION) {\n')
		out.write('\t\t\t\tFile file = main.fc.getSelectedFile();\n')
                out.write('\t\t\t\tif(file.exists()) {\n')
                out.write('\t\t\t\t\tObject[] options = {"Replace",\n')
                out.write('\t\t\t\t\t"Append",\n')
                out.write('\t\t\t\t\t"Cancel"};\n')
                out.write('\t\t\t\t\tint n = JOptionPane.showOptionDialog(main.frame,\n')
                out.write('\t\t\t\t\t"That file already exists.  Do you want to replace or append to "\n')
                out.write('\t\t\t\t\t+ "that file or select another file by pressing cancel?",\n')
                out.write('\t\t\t\t\t"File Already Exists",\n')
                out.write('\t\t\t\t\tJOptionPane.YES_NO_CANCEL_OPTION,\n')
                out.write('\t\t\t\t\tJOptionPane.QUESTION_MESSAGE,\n')
                out.write('\t\t\t\t\tnull,\n')
                out.write('\t\t\t\t\toptions,\n')
                out.write('\t\t\t\t\toptions[0]);\n')
		out.write('\t\t\t\t\tif(n==0) {\n')
		out.write('\t\t\t\t\t\treplace = true;\n')
		out.write('\t\t\t\t\t\tappend = false;\n')
		out.write('\t\t\t\t\t\tmain.%s.setText(file.getAbsolutePath());\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t\telse if(n==1) {\n')
		out.write('\t\t\t\t\t\tappend = true;\n')
		out.write('\t\t\t\t\t\treplace = false;\n')
		out.write('\t\t\t\t\t\tmain.%s.setText(file.getAbsolutePath());\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t\telse {\n')
		out.write('\t\t\t\t\t\treplace = false;\n')
		out.write('\t\t\t\t\t\tappend = false;\n')
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t}\n')
		out.write('\t\t\t\telse\n')
		out.write('\t\t\t\t\tmain.%s.setText(file.getAbsolutePath());\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\t\tmain.fc.setSelectedFile(new File(""));\n')
		out.write('\t\t\t}\n')
		out.write('\t\t}\n')

	def writeActionCmdString(self, out):
		out.write('\t\t\tcommandString.append(" %s = \'"+ %s +"\'");\n' % (self.base.getNexusLabel(), self.base.getValueName()))
		out.write('\t\t\tif(replace)\n')
		out.write('\t\t\t\tcommandString.append(" %sReplace = \'yes\'");\n' % self.prefix)
		out.write('\t\t\telse if(append)\n')
		out.write('\t\t\t\tcommandString.append(" %sAppend = \'yes\'");\n' % self.prefix)

	def writeGetParams(self, out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\t\tOutfileTypeInfo typeInfo = allParams[i].getOutfileTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefStringValueElement value = typeInfo.getDefault();\n')
                out.write('\t\t\tif(value!=null && value.getConstantVal()!=null) {\n')
                out.write('\t\t\t\tmain.%s.setText(value.getConstantVal());\n' % cn)
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

	def writeSetParams(self,out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\t\tOutfileTypeInfo typeInfo = allParams[i].getOutfileTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefStringValueElement value = typeInfo.getDefault();\n')
                out.write('\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tvalue.setConstantVal(%s);\n' % self.base.getValueName())
                out.write('\t\t\t\ttypeInfo.setDefault(value);\n')
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

class CmdParamOutputTypeInfo(CmdParamTypeInfo): 
	def getType(self, which=1):
		if(which==1):
			return 'cb_Out_'		
		elif(which==2):
			return 'tf_Out_'		
		elif(which==3):
			return 'rb_Out_'		

	def getDataType(self, which=1):
		if(which==1):
			return 'Boolean'		
		elif(which==2):
			return 'String'		
		elif(which==3):
			return ''		

	def getSelectedName(self, which='a'):
		type=''
		if(which=='a'):
			type='_File'		
		elif(which=='b'):
			type='_Out'		
		elif(which=='c'):
			type='_Plot'		
		return self.base.getComponentName(self.getType())+type+'_Selected'

	def getValueName(self, which=1):
		type=''
		num=1
		if(which=='a'):
			type='_File'		
		elif(which=='b'):
			type='_Out'		
		elif(which=='c'):
			type='_Plot'		
		elif(which==2):
			num=2;
		elif(which==3):
			num=3;
		return self.base.getComponentName(self.getType(num))+type+'_Value'

	def writeIntro(self, out):
		# File checkbox
		vn = self.getValueName('a')
		sn = self.getSelectedName()
		selected = 'false'
		if len(self.default.file)>0 and self.default.suppress=="false":
			selected = 'true'
		out.write('\tprivate String %s = "";\n' % vn)
		out.write('\tprivate boolean %s = %s;\n' % (sn, selected))
		# Out checkbox
		vn = self.getValueName('b')
		sn = self.getSelectedName('b')
		selected = 'false'
		outputFound = 0
		for redirect in self.default.redirect:
			if redirect=="Output":
				outputFound = 1
		if outputFound and self.default.suppress=="false":
			selected = 'true'
		out.write('\tprivate String %s = "";\n' % vn)
		out.write('\tprivate boolean %s = %s;\n' % (sn, selected))
		# Plot checkbox
		if self.default.plottable=="true":
			vn = self.getValueName('c')
			sn = self.getSelectedName('c')
			selected = 'false'
			plotFound = 0
			for redirect in self.default.redirect:
				if redirect=="Plot":
					plotFound = 1
			if plotFound and self.default.suppress=="false":
				selected = 'true'
			out.write('\tprivate String %s = "";\n' % vn)
			out.write('\tprivate boolean %s = %s;\n' % (sn, selected))
		# File textfield
		vn = self.getValueName(2)
		out.write('\tprivate String %s = "";\n' % vn)
		# Append & Replace radiobuttons
		vn = self.getValueName(3)
		out.write('\tprivate String %s = "";\n' % vn)

	def writeItemStateChanged(self, out):
		out.write('\t\tif( (itemEvent.getStateChange() == java.awt.event.ItemEvent.SELECTED) ')
		out.write('&& (((JCheckBox)itemEvent.getItem()).getName().equals("%s_File")) ) {\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\t%s = true;\n' % self.getSelectedName('a'))
		out.write('\t\t\tmain.choose_%s.setEnabled(true);\n' % self.base.getComponentName(""))
		out.write('\t\t\tmain.%s.setEnabled(true);\n' % self.base.getComponentName(self.getType(2)))
		out.write('\t\t\tmain.%s_Append.setEnabled(true);\n' % self.base.getComponentName(self.getType(3)))
		out.write('\t\t\tmain.%s_Replace.setEnabled(true);\n' % self.base.getComponentName(self.getType(3)))
		out.write('\t\t}\n')
		out.write('\t\tif( (itemEvent.getStateChange() == java.awt.event.ItemEvent.DESELECTED) ')
		out.write('&& (((JCheckBox)itemEvent.getItem()).getName().equals("%s_File")) ) {\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\t%s = false;\n' % self.getSelectedName('a'))
		out.write('\t\t\tmain.choose_%s.setEnabled(false);\n' % self.base.getComponentName(""))
		out.write('\t\t\tmain.%s.setEnabled(false);\n' % self.base.getComponentName(self.getType(2)))
		out.write('\t\t\tmain.%s_Append.setEnabled(false);\n' % self.base.getComponentName(self.getType(3)))
		out.write('\t\t\tmain.%s_Replace.setEnabled(false);\n' % self.base.getComponentName(self.getType(3)))
		out.write('\t\t}\n')
		out.write('\t\tif( (itemEvent.getStateChange() == java.awt.event.ItemEvent.SELECTED) ')
		out.write('&& (((JCheckBox)itemEvent.getItem()).getName().equals("%s_Out")) )\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\t%s = true;\n' % self.getSelectedName('b'))
		out.write('\t\tif( (itemEvent.getStateChange() == java.awt.event.ItemEvent.DESELECTED) ')
		out.write('&& (((JCheckBox)itemEvent.getItem()).getName().equals("%s_Out")) )\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\t%s = false;\n' % self.getSelectedName('b'))
		if self.default.plottable=="true":
			out.write('\t\tif( (itemEvent.getStateChange() == java.awt.event.ItemEvent.SELECTED) ')
			out.write('&& (((JCheckBox)itemEvent.getItem()).getName().equals("%s_Plot")) )\n' % self.base.getComponentName(self.getType()))
			out.write('\t\t\t%s = true;\n' % self.getSelectedName('c'))
			out.write('\t\tif( (itemEvent.getStateChange() == java.awt.event.ItemEvent.DESELECTED) ')
			out.write('&& (((JCheckBox)itemEvent.getItem()).getName().equals("%s_Plot")) )\n' % self.base.getComponentName(self.getType()))
			out.write('\t\t\t%s = false;\n' % self.getSelectedName('c'))

	def writeGetErrMsgs(self, out, outer = None):
		cn = self.base.getComponentName(self.getType(2))
		gd = self.base.getGUIDesc()
	        out.write('\t\tif(tfName.equals("%s"))\n\t\t\t' % cn)
	        if outer == None:
			out.write(r'ret =  "\n%s filename must be text";' % gd);
		else:
			out.write(r'ret =  "\n%s filename must be text";' % outer.base.getGUIDesc());
		out.write('\n');
		cn = self.base.getComponentName(self.getType(3))
	        out.write('\t\tif(tfName.equals("%s_Append") || tfName.equals("%s_Replace"))\n\t\t\t' % (cn,cn))
	        if outer == None:
			out.write(r'ret =  "\nPlease select either Append or Replace for %s";' % gd);
		else:
			out.write(r'ret =  "\nPlease select either Append or Replace for %s";' % outer.base.getGUIDesc());
		out.write('\n');

	def writeActionCmdString(self, out):
		out.write('\t\t\tcommandString.append(" %s = (");\n' % self.base.getNexusLabel())
		out.write('\t\t\tif(%s) {\n' % self.getSelectedName('a'))
		out.write('\t\t\tcommandString.append("file = \'"+ %s +"\'");\n' % self.getValueName(2))
		out.write('\t\t\tif(%s.equals("Append")) \n' % self.getValueName(3))
		out.write('\t\t\tcommandString.append(" append");\n')
		out.write('\t\t\telse if(%s.equals("Replace")) \n' % self.getValueName(3))
		out.write('\t\t\tcommandString.append(" replace");\n')
		out.write('\t\t\t}\n')
		out.write('\t\t\tif(%s) {\n' % self.getSelectedName('b'))
		out.write('\t\t\tcommandString.append(" redirect = Output ");\n')
		out.write('\t\t\t}\n')
		if self.default.plottable=="true":
			out.write('\t\t\tif(%s) {\n' % self.getSelectedName('c'))
			out.write('\t\t\tcommandString.append(" redirect = Plot ");\n')
			out.write('\t\t\t}\n')
			out.write('\t\t\tif(!%s && !%s && !%s) {\n' % (self.getSelectedName('a'),self.getSelectedName('b'),self.getSelectedName('c')))
		else:
			out.write('\t\t\tif(!%s && !%s) {\n' % (self.getSelectedName('a'),self.getSelectedName('b')))
		out.write('\t\t\tcommandString.append("Suppress");\n')
		out.write('\t\t\t}\n')
		out.write('\t\t\tcommandString.append(")");\n')

	def writeActionGetValues(self, out, outer = None):
		cn = self.base.getComponentName(self.getType(2))
		vn = self.getValueName(2)
		out.write('\t\t%s = main.%s.getText();\n' % (vn, cn))


	def writeActionValidation(self, out, outer = None):
		cn = self.base.getComponentName(self.getType(2))
		vn = self.getValueName(2)
		out.write('\t\tif(%s) {\n' % self.getSelectedName('a'))
		out.write('\t\tif(!validate%s(%s)) {\n' % (self.getDataType(2), vn))
		if outer == None:
			out.write('\t\t\tsetValidationError(main.%s, main.%s);\n' % (cn, self.base.getComponentName()))
			out.write('\t\t\tsuccess = false;\n')
		else:
			outer.cmd_param[1].writeActionGetValues(out, outer)
			outer.cmd_param[1].writeActionValidation(out, outer)
		out.write('\t\t}\n')
		out.write('\t\telse\n')
		out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())
		cn = self.base.getComponentName(self.getType(3))
		vn = self.getValueName(3)
		out.write('\t\tif(main.%s_Append.isSelected()) {\n' % (cn))
		out.write('\t\t\t%s = "Append";\n' % (vn))
		if outer == None:
			out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())
		else:
			out.write('\t\t\tresetLabel(main.%s);\n' % outer.cmd_param[0].getComponentName())
		out.write('\t\t}\n')
		out.write('\t\telse if(main.%s_Replace.isSelected()) {\n' % (cn))
		out.write('\t\t\t%s = "Replace";\n' % (vn))
		if outer == None:
			out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())
		else:
			out.write('\t\t\tresetLabel(main.%s);\n' % outer.cmd_param[0].getComponentName())
		out.write('\t\t}\n')
		out.write('\t\telse {\n')
		if outer == None:
			out.write('\t\tsetValidationError(main.%s_Append, main.%s);\n' % (cn,self.base.getComponentName()))
			out.write('\t\tsetValidationError(main.%s_Replace, main.%s);\n' % (cn,self.base.getComponentName()))
		else:
			out.write('\t\tsetValidationError(main.%s_Append, main.%s, main.%s);\n' % (cn,outer.cmd_param[0].getComponentName(outer.cmd_param[0].getType()),outer.cmd_param[0].getComponentName()))
			out.write('\t\tsetValidationError(main.%s_Replace, main.%s, main.%s);\n' % (cn,outer.cmd_param[0].getComponentName(outer.cmd_param[0].getType()),outer.cmd_param[0].getComponentName()))
		out.write('\t\t\tsuccess = false;\n')
		out.write('\t\t}\n')
		out.write('\t\t}\n')
		out.write('\t\telse\n')
		out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())

	def writeGetParams(self, out, type = 'reg'):
		cn1 = self.base.getComponentName(self.getType())+'_File'
		cn1a = self.base.getComponentName(self.getType(3))+'_Append'
		cn1r = self.base.getComponentName(self.getType(3))+'_Replace'
		cn2 = self.base.getComponentName(self.getType())+'_Out'
		cn3 = self.base.getComponentName(self.getType())+'_Plot'
		if self.base.label == "" or self.base.label == None:
			out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s") {\n' % (self.base.label,self.base.placement))
		else:
			out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\t\tOutputTypeInfo typeInfo = allParams[i].getOutputTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tOutputTypeInfo.Default value = typeInfo.getDefault();\n')
                out.write('\t\t\tif(value!=null) {\n')
		out.write('\t\t\t\tif(!value.getSuppress()) {\n')
		out.write('\t\t\t\t\tOutputTypeInfo.Default.File f = value.getFile();\n')
		out.write('\t\t\t\t\tif(f!=null) {\n')
		out.write('\t\t\t\t\t\tmain.%s.setSelected(true);\n' % cn1)
		out.write('\t\t\t\t\t\tmain.%s.setSelected(f.getAppend());\n' % cn1a)
		out.write('\t\t\t\t\t\tmain.%s.setSelected(f.getReplace());\n' % cn1r)
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t\telse\n')
		out.write('\t\t\t\t\t\tmain.%s.setSelected(false);\n' % cn1)
		out.write('\t\t\t\t\tOutputRedirectionEnum.Enum[] redirect = value.getRedirectArray();\n')
		out.write('\t\t\t\t\tboolean foundOut = false;\n')
		out.write('\t\t\t\t\tboolean foundPlot = false;\n')
		out.write('\t\t\t\t\tfor(int r=0;r<redirect.length;r++) {\n')
		out.write('\t\t\t\t\tif(redirect[r]!=null && redirect[r].toString().equals("Output"))\n')
		out.write('\t\t\t\t\tfoundOut = true;\n')
		if self.default.plottable=="true":
			out.write('\t\t\t\t\tif(redirect[r]!=null && redirect[r].toString().equals("Plot"))\n')
			out.write('\t\t\t\t\tfoundPlot = true;\n')
		out.write('\t\t\t\t}\n')
		out.write('\t\t\t\t\tif(foundOut)\n')
		out.write('\t\t\t\t\t\tmain.%s.setSelected(true);\n' % cn2)
		out.write('\t\t\t\t\telse\n')
		out.write('\t\t\t\t\t\tmain.%s.setSelected(false);\n' % cn2)
		if self.default.plottable=="true":
			out.write('\t\t\t\t\tif(foundPlot)\n')
			out.write('\t\t\t\t\t\tmain.%s.setSelected(true);\n' % cn3)
			out.write('\t\t\t\t\telse\n')
			out.write('\t\t\t\t\t\tmain.%s.setSelected(false);\n' % cn3)
		out.write('\t\t\t\t}\n')
                out.write('\t\t\t\telse {\n')
                out.write('\t\t\t\t\tmain.%s.setSelected(false);\n' % cn1)
                out.write('\t\t\t\t\tmain.%s.setSelected(false);\n' % cn2)
		if self.default.plottable=="true":
	                out.write('\t\t\t\t\tmain.%s.setSelected(false);\n' % cn3)
                out.write('\t\t\t\t}\n')
		
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

	def writeSetParams(self,out, type = 'reg'):
		cn1 = self.base.getComponentName(self.getType())+'_File'
		cn1a = self.base.getComponentName(self.getType(3))+'_Append'
		cn1r = self.base.getComponentName(self.getType(3))+'_Replace'
		cn2 = self.base.getComponentName(self.getType())+'_Out'
		cn3 = self.base.getComponentName(self.getType())+'_Plot'
		if self.base.label == "" or self.base.label == None:
			out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s") {\n' % (self.base.label,self.base.placement))
		else:
			out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\t\tOutputTypeInfo typeInfo = allParams[i].getOutputTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tOutputTypeInfo.Default value = typeInfo.getDefault();\n')
                out.write('\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tif(main.%s.isSelected()) {\n' % cn1)
                out.write('\t\t\t\t\tOutputTypeInfo.Default.File f = value.getFile();\n')
                out.write('\t\t\t\t\tif(f==null)\n')
                out.write('\t\t\t\t\t\tf = value.addNewFile();\n')
                out.write('\t\t\t\t\tf.setAppend(main.%s.isSelected());\n' % cn1a)
                out.write('\t\t\t\t\tf.setReplace(main.%s.isSelected());\n' % cn1r)
                out.write('\t\t\t\t\tvalue.setSuppress(false);\n')
                out.write('\t\t\t\t}\n')
                out.write('\t\t\t\telse {\n')
                out.write('\t\t\t\t\tOutputTypeInfo.Default.File f = value.getFile();\n')
                out.write('\t\t\t\t\tif(f!=null)\n')
                out.write('\t\t\t\t\t\tvalue.unsetFile();\n')
                out.write('\t\t\t\t}\n')
                out.write('\t\t\t\tif(main.%s.isSelected()) {\n' % cn2)
		out.write('\t\t\t\t\tOutputRedirectionEnum.Enum[] redirect = value.getRedirectArray();\n')
		out.write('\t\t\t\t\tboolean found = false;\n')
		out.write('\t\t\t\t\tfor(int r=0;r<redirect.length;r++) {\n')
		out.write('\t\t\t\t\t\tif(redirect[r]!=null && redirect[r].toString().equals("Output")) {\n')
                out.write('\t\t\t\t\t\t\tfound = true;\n')
                out.write('\t\t\t\t\t\t\tbreak;\n')
		out.write('\t\t\t\t\t\t}\n')
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t\tif(!found) {\n')
                out.write('\t\t\t\t\t\tvalue.addRedirect(OutputRedirectionEnum.Enum.forString("Output"));\n')
                out.write('\t\t\t\t\t\tvalue.setSuppress(false);\n')
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t}\n')
		out.write('\t\t\t\telse {\n')
		out.write('\t\t\t\tOutputRedirectionEnum.Enum[] redirect = value.getRedirectArray();\n')
		out.write('\t\t\t\tboolean found = false;\n')
		out.write('\t\t\t\tint index = -1;\n')
		out.write('\t\t\t\tfor(int r=0;r<redirect.length;r++) {\n')
		out.write('\t\t\t\t\tif(redirect[r]!=null && redirect[r].toString().equals("Output")) {\n')
                out.write('\t\t\t\t\t\tfound = true;\n')
                out.write('\t\t\t\t\t\tindex = r;\n')
                out.write('\t\t\t\t\t\tbreak;\n')
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t\tif(found) \n')
		out.write('\t\t\t\t\t\tvalue.removeRedirect(index);\n')
		out.write('\t\t\t\t}\n')
		if self.default.plottable=="true":
			out.write('\t\t\t\tif(main.%s.isSelected()) {\n' % cn3)
			out.write('\t\t\t\t\tOutputRedirectionEnum.Enum[] redirect = value.getRedirectArray();\n')
			out.write('\t\t\t\t\tboolean found = false;\n')
			out.write('\t\t\t\t\tfor(int r=0;r<redirect.length;r++) {\n')
			out.write('\t\t\t\t\t\tif(redirect[r]!=null && redirect[r].toString().equals("Plot")) {\n')
			out.write('\t\t\t\t\t\t\tfound = true;\n')
			out.write('\t\t\t\t\t\t\tbreak;\n')
			out.write('\t\t\t\t\t\t}\n')
			out.write('\t\t\t\t\t}\n')
			out.write('\t\t\t\t\tif(!found) {\n')
			out.write('\t\t\t\t\t\tvalue.addRedirect(OutputRedirectionEnum.Enum.forString("Plot"));\n')
			out.write('\t\t\t\t\t\tvalue.setSuppress(false);\n')
			out.write('\t\t\t\t\t}\n')
			out.write('\t\t\t\t}\n')
			out.write('\t\t\t\telse {\n')
			out.write('\t\t\t\tOutputRedirectionEnum.Enum[] redirect = value.getRedirectArray();\n')
			out.write('\t\t\t\tboolean found = false;\n')
			out.write('\t\t\t\tint index = -1;\n')
			out.write('\t\t\t\tfor(int r=0;r<redirect.length;r++) {\n')
			out.write('\t\t\t\t\tif(redirect[r]!=null && redirect[r].toString().equals("Plot")) {\n')
			out.write('\t\t\t\t\t\tfound = true;\n')
			out.write('\t\t\t\t\t\tindex = r;\n')
			out.write('\t\t\t\t\t\tbreak;\n')
			out.write('\t\t\t\t\t}\n')
			out.write('\t\t\t\t\t}\n')
			out.write('\t\t\t\t\tif(found) \n')
			out.write('\t\t\t\t\t\tvalue.removeRedirect(index);\n')
			out.write('\t\t\t\t}\n')
			out.write('\t\t\t\tif(!main.%s.isSelected() && !main.%s.isSelected() && !main.%s.isSelected()) {\n' %(cn1,cn2,cn3))
		else:
			out.write('\t\t\t\tif(!main.%s.isSelected() && !main.%s.isSelected()) {\n' %(cn1,cn2))
                out.write('\t\t\t\t\tvalue.setSuppress(true);\n')
                out.write('\t\t\t\t}\n')
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')
                

	def writeActionFileChooser(self, out):
		out.write('\tif(e.getActionCommand().equals("chooseAction_%s")) {\n' % self.base.getComponentName(''))
		out.write('\t\t\tint retVal = main.fc.showSaveDialog(main.%s);\n' % self.base.getComponentName(self.getType(2)))
                out.write('\t\t\tif (retVal == JFileChooser.APPROVE_OPTION) {\n')
		out.write('\t\t\t\tFile file = main.fc.getSelectedFile();\n')
                out.write('\t\t\t\tif(file.exists()) {\n')
                out.write('\t\t\t\t\tObject[] options = {"Replace",\n')
                out.write('\t\t\t\t\t"Append",\n')
                out.write('\t\t\t\t\t"Cancel"};\n')
                out.write('\t\t\t\t\tint n = JOptionPane.showOptionDialog(main.frame,\n')
                out.write('\t\t\t\t\t"That file already exists.  Do you want to replace or append to "\n')
                out.write('\t\t\t\t\t+ "that file or select another file by pressing cancel?",\n')
                out.write('\t\t\t\t\t"File Already Exists",\n')
                out.write('\t\t\t\t\tJOptionPane.YES_NO_CANCEL_OPTION,\n')
                out.write('\t\t\t\t\tJOptionPane.QUESTION_MESSAGE,\n')
                out.write('\t\t\t\t\tnull,\n')
                out.write('\t\t\t\t\toptions,\n')
                out.write('\t\t\t\t\toptions[0]);\n')
		out.write('\t\t\t\t\tif(n==0) {\n')
		out.write('\t\t\t\t\t\tmain.%s_Append.setSelected(false);\n' % self.base.getComponentName(self.getType(3)))
		out.write('\t\t\t\t\t\tmain.%s_Replace.setSelected(true);\n' % self.base.getComponentName(self.getType(3)))
		out.write('\t\t\t\t\t\tmain.%s.setText(file.getAbsolutePath());\n' % self.base.getComponentName(self.getType(2)))
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t\telse if(n==1) {\n')
		out.write('\t\t\t\t\t\tmain.%s_Append.setSelected(true);\n' % self.base.getComponentName(self.getType(3)))
		out.write('\t\t\t\t\t\tmain.%s_Replace.setSelected(false);\n' % self.base.getComponentName(self.getType(3)))
		out.write('\t\t\t\t\t\tmain.%s.setText(file.getAbsolutePath());\n' % self.base.getComponentName(self.getType(2)))
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t}\n')
		out.write('\t\t\t\telse\n')
		out.write('\t\t\t\t\tmain.%s.setText(file.getAbsolutePath());\n' % self.base.getComponentName(self.getType(2)))
		out.write('\t\t\t\tmain.fc.setSelectedFile(new File(""));\n')
		out.write('\t\t\t}\n')
		out.write('\t\t}\n')

class CmdParamIntegerTypeInfo(CmdParamNumberTypeInfo):
	def getDataType(self):
		return 'Long'		

	def writeIntro(self, out):
		vn = self.base.getValueName()
		out.write('\tprivate String %s = "";\n' % vn)
		out.write('\tprivate long %s_min = -1;\n' % self.base.getGUILabel())
		out.write('\tprivate long %s_max = -1;\n' % self.base.getGUILabel())

	def writeGetParams(self, out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if type == 'reg':
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tIntegerTypeInfo typeInfo = allParams[i].getIntegerTypeInfo();\n')
		else:
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s") && allSubParams[j].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tIntegerTypeInfo typeInfo = allSubParams[j].getIntegerTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefIntegerValueElement valueD = typeInfo.getDefault();\n')
                out.write('\t\t\tif(valueD!=null && valueD.getConstantVal()!=null) {\n')
                out.write('\t\t\t\tmain.%s.setText(valueD.getConstantVal().toString());\n' % cn)
                out.write('\t\t\t}\n')
                out.write('\t\t\tIntegerValueElement value = typeInfo.getMinVal();\n')
                out.write('\t\t\ttry {\n')
                out.write('\t\t\t\t%s_min = Long.MIN_VALUE;\n' % self.base.getGUILabel())
                out.write('\t\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tif(value.getConstantVal()!=null) \n')
                out.write('\t\t\t\t\t%s_min = Long.parseLong(value.getConstantVal().toString());\n' % self.base.getGUILabel())
                out.write('\t\t\t\t}\n')
                out.write('\t\t\t}\n')
                out.write('\t\t\tcatch (Exception e) {\n')
                out.write('\t\t\t\t%s_min = Long.MIN_VALUE;\n' % self.base.getGUILabel())
                out.write('\t\t\t}\n')
                out.write('\t\t\tvalue = typeInfo.getMaxVal();\n')
                out.write('\t\t\ttry {\n')
                out.write('\t\t\t\t%s_max = Long.MAX_VALUE;\n' % self.base.getGUILabel())
                out.write('\t\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tif(value.getConstantVal()!=null) \n')
                out.write('\t\t\t\t\t%s_max = Long.parseLong(value.getConstantVal().toString());\n' % self.base.getGUILabel())
                out.write('\t\t\t\t}\n')
                out.write('\t\t\t}\n')
                out.write('\t\t\tcatch (Exception e) {\n')
                out.write('\t\t\t\t%s_min = Long.MAX_VALUE;\n' % self.base.getGUILabel())
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

	def writeSetParams(self,out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if type == 'reg':
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tIntegerTypeInfo typeInfo = allParams[i].getIntegerTypeInfo();\n')
		else:
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s") && allSubParams[j].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tIntegerTypeInfo typeInfo = allSubParams[j].getIntegerTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefIntegerValueElement valueD = typeInfo.getDefault();\n')
                out.write('\t\t\tif(valueD!=null && !%s.trim().equals("")) {\n' % self.base.getValueName())
                out.write('\t\t\t\tvalueD.setConstantVal(new java.math.BigInteger(%s));\n' % self.base.getValueName())
                out.write('\t\t\t\ttypeInfo.setDefault(valueD);\n')
                out.write('\t\t\t}\n')
                out.write('\t\t\tIntegerValueElement value = typeInfo.getMinVal();\n')
                out.write('\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tvalue.setConstantVal(new java.math.BigInteger(new Long(%s_min).toString()));\n' % self.base.getGUILabel())
                out.write('\t\t\t\ttypeInfo.setMinVal(value);\n')
                out.write('\t\t\t}\n')
                out.write('\t\t\tvalue = typeInfo.getMaxVal();\n')
                out.write('\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tvalue.setConstantVal(new java.math.BigInteger(new Long(%s_max).toString()));\n' % self.base.getGUILabel())
                out.write('\t\t\t\ttypeInfo.setMaxVal(value);\n')
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

class CmdParamDoubleTypeInfo(CmdParamNumberTypeInfo):
	def getDataType(self):
		return 'Double'

	def writeIntro(self, out):
		vn = self.base.getValueName()
		out.write('\tprivate String %s = "";\n' % vn)
		out.write('\tprivate double %s_min = -1;\n' % self.base.getGUILabel())
		out.write('\tprivate double %s_max = -1;\n' % self.base.getGUILabel())

	def writeGetParams(self, out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if type == 'reg':
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tDoubleTypeInfo typeInfo = allParams[i].getDoubleTypeInfo();\n')
		else:
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s") && allSubParams[j].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tDoubleTypeInfo typeInfo = allSubParams[j].getDoubleTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefDoubleValueElement valueD = typeInfo.getDefault();\n')
                out.write('\t\t\tif(valueD!=null && valueD.getConstantVal()!=null) {\n')
                out.write('\t\t\t\tmain.%s.setText(valueD.getConstantVal());\n' % cn)
                out.write('\t\t\t}\n')
                out.write('\t\t\tDoubleValueElement value = typeInfo.getMinVal();\n')
                out.write('\t\t\ttry {\n')
                out.write('\t\t\t\t%s_min = Long.MIN_VALUE;\n' % self.base.getGUILabel())
                out.write('\t\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tif(value.getConstantVal()!=null) \n')
                out.write('\t\t\t\t\t%s_min = Double.parseDouble(value.getConstantVal());\n' % self.base.getGUILabel())
                out.write('\t\t\t\t}\n')
                out.write('\t\t\t}\n')
                out.write('\t\t\tcatch (Exception e) {\n')
                out.write('\t\t\t\t%s_min = Long.MIN_VALUE;\n' % self.base.getGUILabel())
                out.write('\t\t\t}\n')
                out.write('\t\t\tvalue = typeInfo.getMaxVal();\n')
                out.write('\t\t\ttry {\n')
                out.write('\t\t\t\t%s_max = Double.MAX_VALUE;\n' % self.base.getGUILabel())
                out.write('\t\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tif(value.getConstantVal()!=null) \n')
                out.write('\t\t\t\t\t%s_max = Double.parseDouble(value.getConstantVal());\n' % self.base.getGUILabel())
                out.write('\t\t\t\t}\n')
                out.write('\t\t\t}\n')
                out.write('\t\t\tcatch (Exception e) {\n')
                out.write('\t\t\t\t%s_min = Double.MAX_VALUE;\n' % self.base.getGUILabel())
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

	def writeSetParams(self,out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if type == 'reg':
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tDoubleTypeInfo typeInfo = allParams[i].getDoubleTypeInfo();\n')
		else:
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s") && allSubParams[j].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tDoubleTypeInfo typeInfo = allSubParams[j].getDoubleTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tDefDoubleValueElement valueD = typeInfo.getDefault();\n')
                out.write('\t\t\tif(valueD!=null) {\n')
                out.write('\t\t\t\tvalueD.setConstantVal(%s);\n' % self.base.getValueName())
                out.write('\t\t\t\ttypeInfo.setDefault(valueD);\n')
                out.write('\t\t\t}\n')
                out.write('\t\t\tDoubleValueElement value = typeInfo.getMinVal();\n')
                out.write('\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tvalue.setConstantVal(new Double(%s_min).toString());\n' % self.base.getGUILabel())
                out.write('\t\t\t\ttypeInfo.setMinVal(value);\n')
                out.write('\t\t\t}\n')
                out.write('\t\t\tvalue = typeInfo.getMaxVal();\n')
                out.write('\t\t\tif(value!=null) {\n')
                out.write('\t\t\t\tvalue.setConstantVal(new Double(%s_max).toString());\n' % self.base.getGUILabel())
                out.write('\t\t\t\ttypeInfo.setMaxVal(value);\n')
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

class CmdParamSetTypeInfo(CmdParamTypeInfo):
	def getType(self, num=0):
		if num == 0:
			return 'list_Set_'
		elif num == 1:
			return 'cmb_Set_'
		elif num == 2:
			return 'tf_Set_'

	def getValueName(self):
		return self.base.getComponentName(self.getType(2))+'_Value'

	def writeIntro(self, out):
		vn = self.getValueName()
		out.write('\tprivate String %s = "";\n' % vn)

	def writeActionGetValues(self, out, outer = None):
		cn = self.base.getComponentName(self.getType(2))
		vn = self.getValueName()
		out.write('\t\t%s = main.%s.getText();\n' % (vn, cn))

	def writeActionValidation(self, out, outer = None):
		out.write('\t\tif(main.%s.isSelectionEmpty() && %s.equals("")) {\n' % (self.base.getComponentName(self.getType()),self.getValueName()))
		out.write('\t\t\tsetValidationError(main.%s,main.%s);\n' % (self.base.getComponentName(self.getType()),self.base.getComponentName()))
		out.write('\t\t\tsuccess = false;\n')
		out.write('\t\t}\n')
		out.write('\t\telse if(!main.%s.getText().equals("") && !validateString(%s)) {\n' % (self.base.getComponentName(self.getType(2)),self.getValueName()))
		out.write('\t\t\tsetValidationError(main.%s,main.%s);\n' % (self.base.getComponentName(self.getType(2)),self.base.getComponentName()))
		out.write('\t\t\tsuccess = false;\n')
		out.write('\t\t}\n')
		out.write('\t\telse\n')
		out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())

	def writeActionCmdString(self, out):
		out.write('\t\t\tif(%s.equals("")) {\n' % self.getValueName())
		out.write('\t\t\tint[] selectedIndices = main.%s.getSelectedIndices();\n' % self.base.getComponentName(self.getType()))
		out.write('\t\t\tStringBuffer choices = new StringBuffer();\n')
		out.write('\t\t\tfor(int i=0;i<selectedIndices.length;i++) {\n')
		out.write('\t\t\t\tchoices.append(" "+(selectedIndices[i]+1));\n')
		out.write('\t\t\t}\n')
		out.write('\t\t\tcommandString.append(choices.toString());\n')
		out.write('\t\t\t}\n')
		out.write('\t\t\telse {\n')
		out.write('\t\t\tcommandString.append(" "+%s);\n' % self.getValueName())
		out.write('\t\t\t}\n')

	def writeGetErrMsgs(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
	        out.write('\t\tif(tfName.equals("%s"))\n\t\t\t' % cn)
		out.write(r'ret =  "\nA value for %s must be selected.";' % self.base.getGUIDesc());
		out.write('\n');


class CmdParamTaxSetTypeInfo(CmdParamSetTypeInfo):
	def __str__(self): return 'CmdParamTaxSetTypeInfo'

	def getName(self): return 'TaxSet'

	def writeGetParams(self, out, type = 'reg'):
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
		out.write('\t\t\tSetTypeInfo typeInfo = allParams[i].getTaxSetTypeInfo();\n')
		out.write('\t\t\tif(typeInfo==null)\n')
		out.write('\t\t\t\tcontinue;\n')
		out.write('\t\t\tDefStringValueElement value = typeInfo.getDefault();\n')
		out.write('\t\t\tif(value!=null && value.getConstantVal()!=null) {\n')
		out.write('\t\t\t\tString cv = value.getConstantVal();\n')
		out.write('\t\t\t\tPhycasSet[] allSets = main.taxSetMgr.getAllSets();\n')
		out.write('\t\t\t\tboolean found = false;\n')
		out.write('\t\t\t\tfor(int a=0;a<allSets.length;a++) {\n')
		out.write('\t\t\t\t\tif(allSets[a].getLabel().equals("default")) {\n')
		out.write('\t\t\t\t\t\tfound = true;\n')
		out.write('\t\t\t\t\t\tbreak;\n')
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t}\n')
		out.write('\t\t\t\tif(!found && !value.getConstantVal().equals("")) {\n')
		out.write('\t\t\t\t\tPhycasSet s = new PhycasSet("default");\n')
		out.write('\t\t\t\t\ts.setMembers(value.getConstantVal());\n')
		out.write('\t\t\t\t\tmain.taxSetMgr.addSet(s);\n')
		out.write('\t\t\t\t}\n')
		out.write('\t\t\t}\n')
	        if self.check_indeces == 'false':
	                out.write('\t\t\tmain.%s.addItem("TaxSets");\n' % self.base.getComponentName(self.getType(1)))
	        else:
			out.write('\t\t\tmain.%s.setListData(main.taxSetMgr.getListData());\n' % self.base.getComponentName(self.getType()))
			out.write('\t\t\tphycasGUI.swixml.PhycasSet[] allSets = main.taxSetMgr.getAllSets();\n')
			out.write('\t\t\tmain.%s.removeAllItems();\n' % self.base.getComponentName(self.getType(1)))
			out.write('\t\t\tmain.%s.addItem("TaxSets");\n' % self.base.getComponentName(self.getType(1)))
			out.write('\t\t\tfor(int j=0;j<allSets.length;j++) {\n')
			out.write('\t\t\t\tmain.%s.addItem(allSets[j].getLabel());\n' % self.base.getComponentName(self.getType(1)))
			out.write('\t\t\t}\n')
		out.write('\t\t\tif (!value.getConstantVal().equals(""))\n')
		out.write('\t\t\t\tmain.%s.setSelectedItem("default");\n' % self.base.getComponentName(self.getType(1)))
                out.write('\t\t}\n')

	def writeActionComboListSelect(self, out):
		out.write('\tif(e.getActionCommand().equals("%s")) {\n' % self.base.getComponentName(self.getType(1)))
                out.write('\t\t\tString selected = (String)main.%s.getSelectedItem();\n' % self.base.getComponentName(self.getType(1)))
                out.write('\t\t\tif(selected==null) {return;}\n')
                out.write('\t\t\tphycasGUI.swixml.PhycasSet s = main.taxSetMgr.getSet(selected);\n')
                out.write('\t\t\tString[] listData = main.taxSetMgr.getListData();\n')
                out.write('\t\t\tif(s!=null) {\n')
                out.write('\t\t\t\tString[] members = s.getMembers();\n')
                out.write('\t\t\t\tint[] indices = new int[members.length];\n')
                out.write('\t\t\t\tint count = 0;\n')
                out.write('\t\t\t\tfor(int i=0;i<listData.length;i++) {\n')
                out.write('\t\t\t\t\tfor(int j=0;j<members.length;j++) {\n')
                out.write('\t\t\t\t\t\tint index = main.taxSetMgr.getIndexForLabel(listData[i]);\n')
                out.write('\t\t\t\t\t\tif(members[j].equals(listData[i]) || (index!=-1 && new Integer(members[j]).intValue()==index)) {\n')
                out.write('\t\t\t\t\t\t\tindices[count] = i;\n')
                out.write('\t\t\t\t\t\t\tcount++;\n')
                out.write('\t\t\t\t\t\t}\n')
                out.write('\t\t\t\t\t}\n')
                out.write('\t\t\t\t}\n')
                out.write('\t\t\t\tmain.%s.setSelectedIndices(indices);\n'  % self.base.getComponentName(self.getType()))
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')


class CmdParamTreeSetTypeInfo(CmdParamSetTypeInfo):
	def __str__(self): return 'CmdParamTreeSetTypeInfo'

	def getName(self): return 'TreeSet'

	def writeGetParams(self, out, type = 'reg'):
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
		out.write('\t\t\tSetTypeInfo typeInfo = allParams[i].getTreeSetTypeInfo();\n')
		out.write('\t\t\tif(typeInfo==null)\n')
		out.write('\t\t\t\tcontinue;\n')
		out.write('\t\t\tDefStringValueElement value = typeInfo.getDefault();\n')
		out.write('\t\t\tif(value!=null && value.getConstantVal()!=null) {\n')
		out.write('\t\t\t\tString cv = value.getConstantVal();\n')
		out.write('\t\t\t\tPhycasSet[] allSets = main.treeSetMgr.getAllSets();\n')
		out.write('\t\t\t\tboolean found = false;\n')
		out.write('\t\t\t\tfor(int a=0;a<allSets.length;a++) {\n')
		out.write('\t\t\t\t\tif(allSets[a].getLabel().equals("default")) {\n')
		out.write('\t\t\t\t\t\tfound = true;\n')
		out.write('\t\t\t\t\t\tbreak;\n')
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t}\n')
		out.write('\t\t\t\tif(!found && !value.getConstantVal().equals("")) {\n')
		out.write('\t\t\t\t\tPhycasSet s = new PhycasSet("default");\n')
		out.write('\t\t\t\t\ts.setMembers(value.getConstantVal());\n')
		out.write('\t\t\t\t\tmain.treeSetMgr.addSet(s);\n')
		out.write('\t\t\t\t}\n')
		out.write('\t\t\t}\n')
	        if self.check_indeces == 'false':
	                out.write('\t\t\tmain.%s.addItem("TreeSets");\n' % self.base.getComponentName(self.getType(1)))
	        else:
			out.write('\t\t\tmain.%s.setListData(main.treeSetMgr.getListData());\n' % self.base.getComponentName(self.getType()))
			out.write('\t\t\tphycasGUI.swixml.PhycasSet[] allSets = main.treeSetMgr.getAllSets();\n')
			out.write('\t\t\tmain.%s.removeAllItems();\n' % self.base.getComponentName(self.getType(1)))
			out.write('\t\t\tmain.%s.addItem("TreeSets");\n' % self.base.getComponentName(self.getType(1)))
			out.write('\t\t\tfor(int j=0;j<allSets.length;j++) {\n')
			out.write('\t\t\t\tmain.%s.addItem(allSets[j].getLabel());\n' % self.base.getComponentName(self.getType(1)))
			out.write('\t\t\t}\n')
		out.write('\t\t\tif (!value.getConstantVal().equals(""))\n')
		out.write('\t\t\t\tmain.%s.setSelectedItem("default");\n' % self.base.getComponentName(self.getType(1)))
                out.write('\t\t}\n')


	def writeActionComboListSelect(self, out):
		out.write('\tif(e.getActionCommand().equals("%s")) {\n' % self.base.getComponentName(self.getType(1)))
                out.write('\t\t\tString selected = (String)main.%s.getSelectedItem();\n' % self.base.getComponentName(self.getType(1)))
                out.write('\t\t\tif(selected==null) {return;}\n')
                out.write('\t\t\tphycasGUI.swixml.PhycasSet s = main.treeSetMgr.getSet(selected);\n')
                out.write('\t\t\tString[] listData = main.treeSetMgr.getListData();\n')
                out.write('\t\t\tif(s!=null) {\n')
                out.write('\t\t\t\tString[] members = s.getMembers();\n')
                out.write('\t\t\t\tint[] indices = new int[members.length];\n')
                out.write('\t\t\t\tint count = 0;\n')
                out.write('\t\t\t\tfor(int i=0;i<listData.length;i++) {\n')
                out.write('\t\t\t\t\tfor(int j=0;j<members.length;j++) {\n')
                out.write('\t\t\t\t\t\tint index = main.treeSetMgr.getIndexForLabel(listData[i]);\n')
                out.write('\t\t\t\t\t\tif(members[j].equals(listData[i]) || (index!=-1 && new Integer(members[j]).intValue()==index)) {\n')
                out.write('\t\t\t\t\t\t\tindices[count] = i;\n')
                out.write('\t\t\t\t\t\t\tcount++;\n')
                out.write('\t\t\t\t\t\t}\n')
                out.write('\t\t\t\t\t}\n')
                out.write('\t\t\t\t}\n')
                out.write('\t\t\t\tmain.%s.setSelectedIndices(indices);\n'  % self.base.getComponentName(self.getType()))
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

class CmdParamCharSetTypeInfo(CmdParamSetTypeInfo):
	def __str__(self): return 'CmdParamCharSetTypeInfo'

	def getName(self): return 'CharSet'

	def writeGetParams(self, out, type = 'reg'):
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
		out.write('\t\t\tSetTypeInfo typeInfo = allParams[i].getCharSetTypeInfo();\n')
		out.write('\t\t\tif(typeInfo==null)\n')
		out.write('\t\t\t\tcontinue;\n')
		out.write('\t\t\tDefStringValueElement value = typeInfo.getDefault();\n')
		out.write('\t\t\tif(value!=null && value.getConstantVal()!=null) {\n')
		out.write('\t\t\t\tString cv = value.getConstantVal();\n')
		out.write('\t\t\t\tPhycasSet[] allSets = main.charSetMgr.getAllSets();\n')
		out.write('\t\t\t\tboolean found = false;\n')
		out.write('\t\t\t\tfor(int a=0;a<allSets.length;a++) {\n')
		out.write('\t\t\t\t\tif(allSets[a].getLabel().equals("default")) {\n')
		out.write('\t\t\t\t\t\tfound = true;\n')
		out.write('\t\t\t\t\t\tbreak;\n')
		out.write('\t\t\t\t\t}\n')
		out.write('\t\t\t\t}\n')
		out.write('\t\t\t\tif(!found && !value.getConstantVal().equals("")) {\n')
		out.write('\t\t\t\t\tPhycasSet s = new PhycasSet("default");\n')
		out.write('\t\t\t\t\ts.setMembers(value.getConstantVal());\n')
		out.write('\t\t\t\t\tmain.charSetMgr.addSet(s);\n')
		out.write('\t\t\t\t}\n')
		out.write('\t\t\t}\n')
	        if self.check_indeces == 'false':
	                out.write('\t\t\tmain.%s.addItem("CharSets");\n' % self.base.getComponentName(self.getType(1)))
	        else:
			out.write('\t\t\tmain.%s.setListData(main.charSetMgr.getListData());\n' % self.base.getComponentName(self.getType()))
			out.write('\t\t\tphycasGUI.swixml.PhycasSet[] allSets = main.charSetMgr.getAllSets();\n')
			out.write('\t\t\tmain.%s.removeAllItems();\n' % self.base.getComponentName(self.getType(1)))
			out.write('\t\t\tmain.%s.addItem("CharSets");\n' % self.base.getComponentName(self.getType(1)))
			out.write('\t\t\tfor(int j=0;j<allSets.length;j++) {\n')
			out.write('\t\t\t\tmain.%s.addItem(allSets[j].getLabel());\n' % self.base.getComponentName(self.getType(1)))
			out.write('\t\t\t}\n')
		out.write('\t\t\tif (!value.getConstantVal().equals(""))\n')
		out.write('\t\t\t\tmain.%s.setSelectedItem("default");\n' % self.base.getComponentName(self.getType(1)))
                out.write('\t\t}\n')

	def writeActionComboListSelect(self, out):
		out.write('\tif(e.getActionCommand().equals("%s")) {\n' % self.base.getComponentName(self.getType(1)))
                out.write('\t\t\tString selected = (String)main.%s.getSelectedItem();\n' % self.base.getComponentName(self.getType(1)))
                out.write('\t\t\tif(selected==null) {return;}\n')
                out.write('\t\t\tphycasGUI.swixml.PhycasSet s = main.charSetMgr.getSet(selected);\n')
                out.write('\t\t\tString[] listData = main.charSetMgr.getListData();\n')
                out.write('\t\t\tif(s!=null) {\n')
                out.write('\t\t\t\tString[] members = s.getMembers();\n')
                out.write('\t\t\t\tint[] indices = new int[members.length];\n')
                out.write('\t\t\t\tint count = 0;\n')
                out.write('\t\t\t\tfor(int i=0;i<listData.length;i++) {\n')
                out.write('\t\t\t\t\tfor(int j=0;j<members.length;j++) {\n')
                out.write('\t\t\t\t\t\tint index = main.charSetMgr.getIndexForLabel(listData[i]);\n')
                out.write('\t\t\t\t\t\tif(members[j].equals(listData[i]) || (index!=-1 && new Integer(members[j]).intValue()==index)) {\n')
                out.write('\t\t\t\t\t\t\tindices[count] = i;\n')
                out.write('\t\t\t\t\t\t\tcount++;\n')
                out.write('\t\t\t\t\t\t}\n')
                out.write('\t\t\t\t\t}\n')
                out.write('\t\t\t\t}\n')
                out.write('\t\t\t\tmain.%s.setSelectedIndices(indices);\n'  % self.base.getComponentName(self.getType()))
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

class CmdParamComboBox(CmdParamChoiceTypeInfo):
	def __str__(self):
		return 'CmdParamComboBox'

	def getType(self):
		return 'cmb_Choice_'

	def getDataType(self):
		return 'String'

	def getSelectedName(self):
		return self.base.getComponentName(self.getType())+'_Selected'

	def writeIntro(self, out):
		vn = self.base.getValueName()
		sn = self.getSelectedName()
		out.write('\tprivate String %s = "";\n' % vn)

	def writeGetErrMsgs(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
	        out.write('\t\tif(tfName.equals("%s"))\n\t\t\t' % cn)
	        if outer == None:
			out.write(r'ret =  "\nA value for %s must be selected.";' % self.base.getGUIDesc());
		else:
			out.write(r'ret =  "\nA value for %s must be selected.";' % outer.base.getGUIDesc());
		out.write('\n');


	def writeChoice(self,out): 
		cn = self.base.getComponentName(self.getType())
		vn = self.base.getValueName()
		sn = self.getSelectedName()
		out.write('\t\tif(main.%s.getSelectedItem()!=null)\n' % cn)
		out.write('\t\t\t%s = (String)main.%s.getSelectedItem();\n' % (vn, cn))
	
	def writeActionGetValues(self, out, outer = None):
		pass

	def writeActionValidation(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		vn = self.base.getValueName()
		sn = self.getSelectedName()
		self.writeChoice(out)
		out.write('\t\tif(%s.equals("")) {\n' % vn)
		if outer == None:
			out.write('\t\tsetValidationError(main.%s, main.%s);\n' % (cn,self.base.getComponentName()))
		else:
			out.write('\t\tsetValidationError(main.%s, main.%s, main.%s);\n' % (cn,outer.cmd_param[0].getComponentName(outer.cmd_param[0].getType()),outer.cmd_param[0].getComponentName()))
		out.write('\t\t\tsuccess = false;\n')
		out.write('\t\t}\n')
		out.write('\t\telse\n')
		if outer == None:
			out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())
		else:
			out.write('\t\t\tresetLabel(main.%s);\n' % outer.cmd_param[0].getComponentName())

	def writeGetParams(self, out, type = 'reg'):
		cn = self.base.getComponentName(self.getType())
		if type == 'reg':
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tChoiceTypeInfo typeInfo = allParams[i].getChoiceTypeInfo();\n')
		else:
			if self.base.label == '' or self.base.label == None:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s") && allSubParams[j].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
			else:
				out.write('\t\tif(allSubParams[j].getLabel().equals("%s")) {\n' % self.base.label)
	                out.write('\t\t\tChoiceTypeInfo typeInfo = allSubParams[j].getChoiceTypeInfo();\n')
                out.write('\t\t\tif(typeInfo==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tStringValueList value = typeInfo.getChoices();\n')
                out.write('\t\t\tif(value==null)\n')
                out.write('\t\t\t\tcontinue;\n')
                out.write('\t\t\tString[] values = value.getConstantValArray();\n')
                out.write('\t\t\tif(values!=null) {\n')
                out.write('\t\t\t\tmain.%s.removeAllItems();\n' % cn)
		if type != 'reg':
	                out.write('\t\t\t\tmain.%s.addItem("");\n' % cn)
                out.write('\t\t\t\tfor(int k=0;k<values.length;k++) {\n')
                out.write('\t\t\t\t\tSystem.out.println("Values["+k+"] = "+values[k]);\n')
                out.write('\t\t\t\t\tmain.%s.addItem(values[k]);\n' % cn)
		out.write('\t\t\t\t}\n')
                out.write('\t\t\t\tmain.%s.setSelectedItem("%s".toUpperCase());\n' % (cn, self.default))
                out.write('\t\t\t}\n')
                out.write('\t\t}\n')

class CmdParamRadioButton(CmdParamChoiceTypeInfo):
	def __str__(self):
		return 'CmdParamRadioButton'

	def getType(self):
		return 'rb_Choice_'

	def getSelectedName(self):
		return self.base.getComponentName(self.getType())+'_Selected'

	def writeIntro(self, out):
		vn = self.base.getValueName()
		out.write('\tprivate String %s = "";\n' % vn)

	def writeGetErrMsgs(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
	        out.write('\t\tif(tfName.equals("%s_0"))\n\t\t\t' % cn)
	        if outer == None:
			out.write(r'ret =  "\nA value for %s must be selected.";' % self.base.getGUIDesc());
		else:
			out.write(r'ret =  "\nA value for %s must be selected.";' % outer.base.getGUIDesc());
		out.write('\n');

	def writeChoice(self,out,outer = None):
		cn = self.base.getComponentName(self.getType())
		vn = self.base.getValueName()
		sn = self.getSelectedName()
		listOfChoices = []
		for i in range(0, len(self.choices.constant_val)):
			listOfChoices.append(self.choices.constant_val[i]) 
		for i in range(0, len(listOfChoices)):
			if i == 0:
				out.write('\t\tif(main.%s_%s.isSelected()) {\n' % (cn, i))
			else:
				out.write('\t\telse if(main.%s_%s.isSelected()) {\n' % (cn, i))
			out.write('\t\t\t%s = "%s";\n' % (vn, listOfChoices[i]))
			if outer == None:
				out.write('\t\t\tresetLabel(main.%s);\n' % self.base.getComponentName())
			else:
				out.write('\t\t\tresetLabel(main.%s);\n' % outer.cmd_param[0].getComponentName())
			out.write('\t\t}\n')
		
	def writeActionGetValues(self, out, outer = None):
		pass

	def writeActionValidation(self, out, outer = None):
		cn = self.base.getComponentName(self.getType())
		vn = self.base.getValueName()
		sn = self.getSelectedName()
		self.writeChoice(out,outer)
		out.write('\t\telse {\n')
		if outer == None:
			out.write('\t\tsetValidationError(main.%s_0, main.%s);\n' % (cn,self.base.getComponentName()))
		else:
			out.write('\t\tsetValidationError(main.%s_0, main.%s, main.%s);\n' % (cn,outer.cmd_param[0].getComponentName(outer.cmd_param[0].getType()),outer.cmd_param[0].getComponentName()))
		out.write('\t\t\tsuccess = false;\n')
		out.write('\t\t}\n')

class CmdParamMixedTypeInfo(CmdParamTypeInfo):
	def __str__(self):
		return 'CmdParamMixedTypeInfo' + ''.join([str(x) for x in self.cmd_param])

	def getType(self):
		return self.cmd_param[0].type_info.getType()

	def writeActionCmdString(self, out):
		out.write('\t\t\tcommandString.append(" %s = ");\n' % self.base.getNexusLabel())
		self.writeCommandStringOptions(out)

	def writeCommandStringOptions(self, out):
		if self.default.constant_val == '0':
			out.write('\t\t\tif(!main.%s.getText().trim().equals("")) \n' % self.cmd_param[0].getComponentName(self.cmd_param[0].getType()))
			out.write('\t\t\t\tcommandString.append(%s);\n' % self.cmd_param[0].getValueName())
			out.write('\t\t\telse\n')
			out.write('\t\t\t\tcommandString.append(%s);\n' % self.cmd_param[1].getValueName())
		else:
			if self.cmd_param[1].type_info.__class__ == CmdParamRadioButton:
				listOfChoices = []
				for i in range(0, len(self.cmd_param[1].type_info.choices.constant_val)):
					listOfChoices.append(self.cmd_param[1].type_info.choices.constant_val[i]) 
				for i in range(0, len(listOfChoices)):
					if i == 0:
						out.write('\t\t\tif(main.%s_%s.isSelected()' % (self.cmd_param[1].getComponentName(self.cmd_param[1].getType()), i))
					else: 
						out.write(' || main.%s_%s.isSelected()' % (self.cmd_param[1].getComponentName(self.cmd_param[1].getType()), i))
				out.write(')\n')
			else: #ComboBox
				out.write('\t\t\tif(!%s.equals(""))\n' % self.cmd_param[1].getValueName())
			out.write('\t\t\t\tcommandString.append(%s);\n' % self.cmd_param[1].getValueName())
			out.write('\t\t\telse\n')
			out.write('\t\t\t\tcommandString.append(%s);\n' % self.cmd_param[0].getValueName())
			
			

	def writeActionGetValues(self, out, outer = None):
		self.cmd_param[0].writeActionGetValues(out, self)

	def writeActionValidation(self, out, outer = None):
		self.cmd_param[0].writeActionValidation(out, self)
		self.cmd_param[1].writeChoice(out)

	def writeIntro(self, out):
		for cmdParam in self.cmd_param:
			cmdParam.writeIntro(out)

	def writeGetErrMsgs(self, out, outer = None):
		for cmdParam in self.cmd_param:
			cmdParam.writeGetErrMsgs(out,self)

	def writeGetParams(self, out, type = 'reg'):
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\tCmdParam[] allSubParams = allParams[i].getMixedTypeInfo().getCmdParamArray();\n')
                out.write('\t\tfor(int j=0;j<allSubParams.length;j++) {\n')
		for cmdParam in self.cmd_param:
			cmdParam.writeGetParams(out, 'sub')
                out.write('\t\t}\n')
                out.write('\t\t}\n')

	def writeSetParams(self,out, type = 'reg'):
		if self.base.label == '' or self.base.label == None:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s") && allParams[i].getPlacement().toString().equals("%s")) {\n' % (self.base.label,self.base.placement))
	        else:
	                out.write('\t\tif(allParams[i].getLabel().equals("%s")) {\n' % self.base.label)
                out.write('\t\tCmdParam[] allSubParams = allParams[i].getMixedTypeInfo().getCmdParamArray();\n')
                out.write('\t\tfor(int j=0;j<allSubParams.length;j++) {\n')
		for cmdParam in self.cmd_param:
			cmdParam.writeSetParams(out, 'sub')
                out.write('\t\t}\n')
                out.write('\t\t}\n')


class FileOutputElement(SAXConstructible): pass
class OutputValue(SAXConstructible): pass

