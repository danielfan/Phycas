#! /usr/bin/python2.2
from sax_constructible import *
import sys, re

if not writingNCL:
	class CommandTestInfo(SAXConstructible): pass
	class CommandRequireOrAvail(SAXConstructible):pass
	class CommandTest(SAXConstructible):pass
	class CommandTestOperand(SAXConstructible):pass
	class CommandTestMessage(SAXConstructible):pass
class CommandTestPredicate(SAXConstructible):pass
class CommandCallback(SAXConstructible):pass

CommandRequireOrAvail._CTH = ChildrenToHandle(
	attr = [	'implemented_at'],
	singleEl = {'test_info' : CommandTestInfo})
CommandTestInfo._CTH = ChildrenToHandle(
	singleEl = {'message' : CommandTestMessage,
				'label' : TextOnlyElement,
				'callback': CommandCallback,
				'logicalOp': CommandTest,
				'comparison': CommandTest,
				'predicate': CommandTestPredicate}) 
CommandTestOperand._CTH = ChildrenToHandle(
	singleEl = {'label' : TextOnlyElement,
				'constant_val' : TextOnlyElement,
				'labile':	PhycasLabileListImpl}) # note PhycasLabileListImpl should be renamed

CommandTestPredicate._CTH = ChildrenToHandle(
	attr = [	'evaluate_to'],
	singleEl = {'label' : TextOnlyElement,
				'constant_val' : TextOnlyElement})
CommandTestMessage._CTH = ChildrenToHandle(
	singleEl = {'failure' : TextOnlyElement,
				'success' : TextOnlyElement})
CommandTest._CTH = ChildrenToHandle(
	attr = [	'name'], 
	singleEl = {'operator' : TextOnlyElement,
				'left_operand' : CommandTestOperand,
				'right_operand' : CommandTestOperand})
CommandCallback_CTH = ChildrenToHandle(
	singleEl = {'function' : TextOnlyElement})

CmdParam._CTH =ChildrenToHandle( 	
	attr = [	'gui_order', 
				'label',
				'persistent',
				'placement'],
	singleEl = {'user_level'	 : TextOnlyElement,
				'description'	 : TextOnlyElement,
				'gui_description'	 : TextOnlyElement,
				'bool_type_info' :	CmdParamBoolTypeInfo,
				'char_set_type_info' :	CmdParamCharSetTypeInfo,
				'choice_type_info' :	CmdParamChoiceTypeInfo,
				'double_type_info' :	CmdParamDoubleTypeInfo,
				'infile_type_info' :	CmdParamInfileTypeInfo,
				'integer_type_info' :	CmdParamIntegerTypeInfo,
				'outfile_type_info' :	CmdParamOutfileTypeInfo,
				'output_type_info' :	CmdParamOutputTypeInfo,
				'restricted_string_type_info' :	CmdParamRestrictedStringTypeInfo,
				'string_type_info' :	CmdParamNxsStringTextField,
				'name_type_info' :	CmdParamNameTypeInfo,
				'tax_set_type_info' :	CmdParamTaxSetTypeInfo,
				'tree_set_type_info' :	CmdParamTreeSetTypeInfo,
				'phycas_impl' : CmdParamPhycasImpl,
				'distribution_type_info' : CmdParamDistributionTypeInfo,
				'mixed_type_info' : CmdParamMixedTypeInfo
				},
	multiEl = {	'cmd_param' : CmdParam, 
				'availability' : CommandRequireOrAvail,
				'requirement' : CommandRequireOrAvail})
#if writingNCL:
#CmdParam._CTH.singleElementDict['distribution_type_info'] = CmdParamDistributionTypeInfo
DistribRangeConstraint._CTH = ChildrenToHandle(
	attr = [	'constraint'], 
	singleEl = {'min_val' : GenericSingleVal, 
				'max_val' : GenericSingleVal})		
CmdParamDistributionTypeInfo._CTH = ChildrenToHandle(
	attr = [	'distrib_class'], 
	singleEl = {'default': GenericSingleVal,
				'num_variates': TextOnlyElement,
				'range_constraint': DistribRangeConstraint, 
				'phycas_impl' : PhycasCallback})
#else:
	#CmdParam._CTH.singleElementDict['distribution_type_info'] = CmdParamNxsStringTextField

FileOutputElement._CTH = ChildrenToHandle(attr = ['append', 'replace', 'path']) 
	
OutputValue._CTH = ChildrenToHandle(
	attr = [	'suppress','plottable'], 
	multiEl = {	'redirect': TextOnlyElement,
				'file':	FileOutputElement})
CmdParamBoolTypeInfo._CTH = ChildrenToHandle(
	singleEl = {'default' : GenericSingleVal})
CmdParamNumberTypeInfo._CTH = ChildrenToHandle(
	singleEl = {'default' : GenericSingleVal, 
				'min_val' : GenericSingleVal, 
				'max_val' : GenericSingleVal})
CmdParamStringTextField._CTH = ChildrenToHandle(
	singleEl = {'default' : GenericSingleVal})
CmdParamRestrictedStringTypeInfo._CTH = CmdParamStringTextField._CTH.getUnion(ChildrenToHandle(
	singleEl = {'default' : GenericSingleVal, 
	 			'disallowed_values' : PhycasStringList}))
CmdParamNameTypeInfo._CTH = CmdParamStringTextField._CTH.getUnion(ChildrenToHandle(
	singleEl = {'default' : GenericSingleVal, 
	 			'disallowed_values' : PhycasStringList}))
CmdParamChoiceTypeInfo._CTH = CmdParamStringTextField._CTH.getUnion(ChildrenToHandle(
	singleEl = {'default' : GenericSingleVal, 
				'choices' : PhycasStringList,
				'phycas_impl' : CmdParamPhycasImpl}))
CmdParamOutfileTypeInfo._CTH = CmdParamStringTextField._CTH.getUnion(ChildrenToHandle(
	singleEl = {'prefix' : TextOnlyElement}))
CmdParamOutputTypeInfo._CTH = ChildrenToHandle(
	singleEl = {'default' : OutputValue})
CmdParamSetTypeInfo._CTH = ChildrenToHandle(attr = ['check_indeces'],
	singleEl = {'default' : GenericSingleVal,
				'phycas_impl' : PhycasSetImpl})
CmdParamMixedTypeInfo._CTH = ChildrenToHandle(
	singleEl = {'default' : GenericSingleVal},
	multiEl = {	'cmd_param' : CmdParam})
CommandList._CTH = ChildrenToHandle(
	multiEl = {'command_set' : CommandSet})
CommandSet._CTH = ChildrenToHandle(
	attr = [	'name','user_interface'],
	multiEl = {'command' : Command})
CmdParamGroup._CTH = ChildrenToHandle(
	attr = [	'name'],
	multiEl = {	'cmd_param' : CmdParam,})
Command._CTH = ChildrenToHandle( 	
	attr = [	'gui_order', 
				'compilation_target',
				'label'],
	singleEl = {'description': TextOnlyElement,
				'phycas_impl' : CmdPhycasImpl},
	multiEl = {	'cmd_param_group' : CmdParamGroup,
				'test' : CommandTest,
				'requirement' : CommandRequireOrAvail,
				'availability' : CommandRequireOrAvail
				})
if writingNCL:
	CmdPhycasImpl._CTH = ChildrenToHandle(	
		singleEl = {'setup' : PhycasCmdSetupInfo,
					'handler': PhycasCallback})
	PhycasCmdSetupInfo._CTH = PhycasCallback._CTH.getUnion(ChildrenToHandle( 
		multiEl = {	'setup_header' : TextOnlyElement}))
else:
	CmdPhycasImpl._CTH = ChildrenToHandle()
