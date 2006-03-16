#! /usr/bin/python2.2
from sax_constructible import SAXConstructible, ChildrenToHandle, TextOnlyElement, IgnoreContentElement

class PhycasCallback(SAXConstructible):
	def getSetupRequiredHeader(self):
		headerList = [self.receiver_class == 'PhoFloor' and 'phycas/floor.hpp' or 'ncl/command/nxs_command_manager.hpp']
		headerList.extend(self.setup_header)   # automaticallyConvertToString = false -> headerList.extend([ str(i) for i in self.setup_header])
		return headerList
	def getReceiverClass(self):
		return self.receiver_class
	def getFuncName(self):
		return self.function
	def getPtrName(self):
		return self.ptr_name_in_setup
	def getBoostFuncDef(self):
		return 'boost::bind(&%s::%s, %s)' % (self.getReceiverClass(), self.getFuncName(), self.ptr_name_in_setup)

class PhycasLabileListImpl(SAXConstructible):
	def __init__(self, elName, parseContext):
		super(PhycasLabileListImpl, self).__init__(elName, parseContext, PhycasLabileListImpl._CTH)


class GenericValBase(SAXConstructible):
	def getValueAsString(self):
		if self.isVariable: # isVariable is set in endSelfElement of derived classes GenericSingleVal PhycasStringList
			if self.label != None: raise TypeError, 'labels are not supported in choices or restricted strings yet'
			return self.labile.phycas_impl.getBoostFuncDef()
		else:
			return '"' + '|'.join(self.constant_val) + '"'

class GenericSingleVal(GenericValBase):
	def endSelfElement(self, name):
		self.isVariable = (self.constant_val == None)
		if self.isVariable: # need to handle empty GenericSingleVal elements.  we'll consider them Not variable with '' as the constant val
			if self.__dict__.get('label') == None and self.__dict__.get('labile') == None:
				self.constant_val = ''
				self.isVariable = False
			
class GenericMultiVal(GenericValBase):
	def endSelfElement(self, name):
		self.isVariable = len(self.constant_val) == 0
		if self.isVariable: # need to handle empty GenericSingleVal elements.  we'll consider them Not variable with [] as the constant val
			if self.__dict__.get('label') == None and self.__dict__.get('labile') == None:
				self.constant_val = []
				self.isVariable = False

class CmdParamPhycasImpl(SAXConstructible):
	def getManipulatedVariableName(self): return self.manipulated_var # automaticallyConvertToString = fale -> self.manipulated_var.getRawChars()
	
class PhycasStringList(GenericMultiVal):pass
class PhycasSetImpl(SAXConstructible):pass

PhycasCallback._CTH = ChildrenToHandle(	
	singleEl = {'receiver_class' : TextOnlyElement,
				'function': TextOnlyElement,
				'ptr_name_in_setup' : TextOnlyElement})
PhycasLabileListImpl._CTH = ChildrenToHandle( 
	singleEl = {'phycas_impl' : PhycasCallback} )
GenericValBase._CTH = ChildrenToHandle(
	singleEl = {'labile' : PhycasLabileListImpl,
				'label' : TextOnlyElement} )
GenericSingleVal._CTH = GenericValBase._CTH.getUnion(ChildrenToHandle(
	singleEl = {'constant_val' : TextOnlyElement} ))
GenericMultiVal._CTH = GenericValBase._CTH.getUnion(ChildrenToHandle(
	multiEl = {	'constant_val' : TextOnlyElement} ))
CmdParamPhycasImpl._CTH = ChildrenToHandle(
	singleEl = {'manipulated_var' : TextOnlyElement})
PhycasSetImpl._CTH = ChildrenToHandle(
	singleEl = {'manager_ptr' : TextOnlyElement})

	
