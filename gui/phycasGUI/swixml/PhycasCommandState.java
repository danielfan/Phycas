package phycasGUI.swixml;

import org.phycas.commandState.*;
import org.phycas.commandLanguage.*;
import java.lang.reflect.*;
import phycasGUI.swixml.handler.*;

/**
 * Parses command_state updates from the main program to the gui.
 * It will need to be modified any time the XML schema's are modified.
 * @author Vanessa Jackson
 */
public class PhycasCommandState {
    
    private PhycasMain main = null;
    private org.phycas.commandState.CommandState commandState = null;
    
    /**
     * Creates a new instance of PhycasCommandState
     * @param m instance of main gui window
     */
    public PhycasCommandState(PhycasMain m) {
        main = m;
    }
    
    /**
     * Parse the command_state message and update the command language XMLBeans
     * @param d XML Document to parse
     */
    public void parse(org.w3c.dom.Document d) {
        String commandLabel = null;
        String paramLabel = null;
        try {
            // this class reflects the base tag of the xml doc
            org.phycas.commandState.CommandStateDocument commandStateDoc = org.phycas.commandState.CommandStateDocument.Factory.parse(d);
            commandState = commandStateDoc.getCommandState();
            
            // get all Commands in CommandState
            org.phycas.commandState.Command[] allCommands = commandState.getCommandArray();
            main.disableAllCommands();
            for(int i=0;i<allCommands.length;i++) {
                commandLabel = allCommands[i].getLabel();
                if(allCommands[i].getAvailable())
                    main.enableAvailableCommand(commandLabel);
                else
                    continue;

                // find method to get Handler class for command
                Class mainClass = main.getClass();
                String methodName = "get"+commandLabel+"Handler";
                Method method = null;
                try {
                    method = mainClass.getMethod(methodName,null);
                }
                catch (NoSuchMethodException e) {
                    continue; // skip if command does not exist
                }
                Object handlerInstance = method.invoke(main, null);
                Class handlerClass = handlerInstance.getClass();

                // get all Cmd Params in Command
                org.phycas.commandState.CmdParam[] allParams = allCommands[i].getCmdParamArray();
                for(int j=0;j<allParams.length;j++) {
                    paramLabel = allParams[j].getLabel();
                    
                    Object csTypeInfo = null;
                    Object clTypeInfo = null;
                    
                    // find appropriate TypeInfo in command state CmdParam
                    Class paramClass = allParams[j].getClass();
                    Method[] allMethods = paramClass.getDeclaredMethods();
                    String csParamMethodName = "";
                    boolean foundCSTI = false;
                    for(int k=0;k<allMethods.length;k++) {
                        if(foundCSTI) break;
                        String paramMethodName = allMethods[k].getName();
                        if(paramMethodName.startsWith("get") && paramMethodName.indexOf("TypeInfo")!=-1) {
                            Object csTI = allMethods[k].invoke(allParams[j], null);
                            if(csTI!=null) {
                                Method[] typeInfoMethods = csTI.getClass().getDeclaredMethods();
                                csTypeInfo = csTI;
                                csParamMethodName = paramMethodName;
                                foundCSTI=true;
                            }
                        }
                    }
                    
                    // find appropriate TypeInfo in command language CmdParam
                    Method getParams = handlerClass.getMethod("getAllParams", null);
                    org.phycas.commandLanguage.CmdParam[] clCmdParams = (org.phycas.commandLanguage.CmdParam[])getParams.invoke(handlerInstance,null);
                    // check if param should be disabled
                    boolean cmdParamAvail = true;
                    if(!allParams[j].getAvailable()) {
                        cmdParamAvail = false;
                    }
                    boolean processed = false;
                    for(int l=0;l<clCmdParams.length;l++) {
                        boolean ready = false;
                        String label = clCmdParams[l].getLabel();
                        String clPlacement = "";
                        String csPlacement = "";
                        if(label.equals("")) { // if no label, use placement instead
                            clPlacement = clCmdParams[l].getPlacement().toString();
                            csPlacement = allParams[j].getPlacement().toString();
                        }
                        if(clCmdParams[l].getMixedTypeInfo()!=null) {
                            org.phycas.commandLanguage.CmdParam[] mixedParams = clCmdParams[l].getMixedTypeInfo().getCmdParamArray();
                            for(int m=0;m<mixedParams.length;m++) {
                                label = mixedParams[m].getLabel();
                                if((!label.equals("") && label.equals(paramLabel)) || (label.equals("") && clPlacement.equals(csPlacement))) {
                                    ready = true;
                                }
                            }
                        }
                        else {
                            if((!label.equals("") && label.equals(paramLabel)) || (label.equals("") && clPlacement.equals(csPlacement))) {
                                ready = true;
                            }
                        }
                        if(ready) { // have the right command param
                            if(paramLabel.equals(""))
                                paramLabel = allParams[j].getPlacement().toString();
                            if(!cmdParamAvail)
                                disableParam(commandLabel, clCmdParams[l].getPhycasImpl().getManipulatedVar(), main, main.getClass());
                            paramClass = clCmdParams[l].getClass();
                            allMethods = paramClass.getDeclaredMethods();
                            boolean foundCLTI = false;
                            for(int k=0;k<allMethods.length;k++) {
                                if(foundCLTI) break;
                                String paramMethodName = allMethods[k].getName();
                                if(paramMethodName.startsWith("get") && paramMethodName.indexOf("TypeInfo")!=-1) {
                                    Object clTI = allMethods[k].invoke(clCmdParams[l], null);
                                    if(clTI!=null && paramMethodName.equals("getMixedTypeInfo")) {
                                        org.phycas.commandLanguage.CmdParam[] mixedParams = clCmdParams[l].getMixedTypeInfo().getCmdParamArray();
                                        for(int m=0;m<mixedParams.length;m++) {
                                            if(foundCLTI) break;
                                            paramClass = mixedParams[m].getClass();
                                            Method[] allMixedMethods = paramClass.getDeclaredMethods();
                                            for(int n=0;n<allMixedMethods.length;n++) {
                                                if(foundCLTI) break;
                                                String paramMixedMethodName = allMixedMethods[n].getName();
                                                if(paramMixedMethodName.startsWith("get") && paramMixedMethodName.indexOf("TypeInfo")!=-1) {
                                                    clTI = allMixedMethods[n].invoke(mixedParams[m], null);
                                                    if(clTI!=null && paramMixedMethodName.equals(csParamMethodName)) {
                                                        Method[] typeInfoMethods = clTI.getClass().getDeclaredMethods();
                                                        clTypeInfo = clTI;
                                                        foundCLTI = true;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else if(clTI!=null) {
                                        Method[] typeInfoMethods = clTI.getClass().getDeclaredMethods();
                                        clTypeInfo = clTI;
                                        foundCLTI = true;
                                    }
                                }
                            }
                            // set values from command state to command language
                            if(csTypeInfo!=null && clTypeInfo!=null) {
                                processed = setFields(csTypeInfo,clTypeInfo, commandLabel);
                            }
                        }
                    }
                    if(!processed) {
                        System.out.println("Command State Error: command "+commandLabel+", param "+paramLabel+" not processed.");
                    }
                }
            }
            // set PhycasSetManager values
            SetManager charSetMgr = commandState.getCharSetManager();
            if(charSetMgr!=null)
                setManager(charSetMgr, "char");
            SetManager taxSetMgr = commandState.getTaxSetManager();
            if(taxSetMgr!=null)
                setManager(taxSetMgr, "tax");
            SetManager treeSetMgr = commandState.getTreeSetManager();
            if(treeSetMgr!=null)
                setManager(treeSetMgr, "tree");
        }
        catch (Exception e) {
            main.taMsgs.append("Unable to parse new command state.\n");
            main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
            System.out.println("Command State Error: command "+commandLabel+", param "+paramLabel+" not processed.");
            e.printStackTrace();
        }
    }
    
    // move values from command state into command language objects
    private boolean setFields(Object csTypeInfo, Object clTypeInfo, String commandLabel) {
        Class csTypeInfoClass = csTypeInfo.getClass();
        String csClassName = csTypeInfoClass.getName();
        Class clTypeInfoClass = clTypeInfo.getClass();
        String clClassName = clTypeInfoClass.getName();
        Method[] csMethods = csTypeInfoClass.getDeclaredMethods();
        Method[] clMethods = clTypeInfoClass.getDeclaredMethods();
        
        boolean success = true;
        
        // must be handled separately because commandState and commandLanguage
        // structures are too different
        if(csTypeInfoClass.getName().equals("org.phycas.commandState.impl.ChoiceTypeInfoImpl")) {
            try {
                Method csMethod = csTypeInfoClass.getMethod("getChoiceArray",null);
                String setMethodName = "setChoices";
                Object csGetObj = csMethod.invoke(csTypeInfo,null);
                Method clMethod = clTypeInfoClass.getMethod("getChoices",null);
                Object clGetObj = clMethod.invoke(clTypeInfo,null);
                
                setCLtoCS(clMethods, setMethodName, clGetObj, csGetObj, clTypeInfo, csTypeInfo);
            }
            catch(Exception e) {
                e.printStackTrace();
                success = false;
            }
        }
        else if(csTypeInfoClass.getName().equals("org.phycas.commandState.impl.RestrictedStringTypeInfoImpl")) {
            try {
                Method csMethod = csTypeInfoClass.getMethod("getDisallowedValueArray",null);
                String setMethodName = "setDisallowedValues";
                Object csGetObj = csMethod.invoke(csTypeInfo,null);
                Method clMethod = clTypeInfoClass.getMethod("getDisallowedValues",null);
                Object clGetObj = clMethod.invoke(clTypeInfo,null);
                
                setCLtoCS(clMethods, setMethodName, clGetObj, csGetObj, clTypeInfo, csTypeInfo);
            }
            catch(Exception e) {
                e.printStackTrace();
                success = false;
            }
        }
        else { // all others
            success = setFieldsGeneral(csTypeInfo, clTypeInfo, commandLabel);
        }
        return success;
    }
    
    private boolean setFieldsGeneral(Object csTypeInfo, Object clTypeInfo, String commandLabel) {
        Class csTypeInfoClass = csTypeInfo.getClass();
        String csClassName = csTypeInfoClass.getName();
        Class clTypeInfoClass = clTypeInfo.getClass();
        String clClassName = clTypeInfoClass.getName();
        Method[] csMethods = csTypeInfoClass.getDeclaredMethods();
        Method[] clMethods = clTypeInfoClass.getDeclaredMethods();
        
        boolean success = true;
        for(int k=0;k<csMethods.length;k++) {
            String methodName = csMethods[k].getName();
            if(methodName.startsWith("get")) { 
                String setMethodName = "set"+methodName.substring(methodName.indexOf("get")+3);
                try {
                    Method clGetMethod = clTypeInfoClass.getMethod(methodName,null);
                    Object clGetObj = clGetMethod.invoke(clTypeInfo,null);
                    Method csGetMethod = csTypeInfoClass.getMethod(methodName,null);
                    Object csGetObj = csMethods[k].invoke(csTypeInfo,null);
                    
                    if(setMethodName.indexOf("Min")!=-1) {
                        // may have label for value, need to translate
                        boolean label = false;
                        try {
                            int test = new Integer(csGetObj.toString()).intValue();
                        }
                        catch (Exception e) {
                            label = true;
                        }
                        if(label)
                            csGetObj = translateLabelToValue("Min", csGetObj, commandLabel);
                    }
                    else if(setMethodName.indexOf("Max")!=-1) {
                        // may have label for value, need to translate
                        boolean label = false;
                        try {
                            int test = new Integer(csGetObj.toString()).intValue();
                        }
                        catch (Exception e) {
                            label = true;
                        }
                        if(label)
                            csGetObj = translateLabelToValue("Max", csGetObj, commandLabel);
                    }
                    setCLtoCS(clMethods, setMethodName, clGetObj, csGetObj, clTypeInfo, csTypeInfo);
                }
                catch(Exception e) {
                    System.out.println("Error calling method "+setMethodName+" on "+clClassName+" of command "+commandLabel);
                    e.printStackTrace();
                    success = false;
                }
            }
        }
        return success;
    }
    
    // move values from command state into command language objects
    private void setCLtoCS(Method[] clMethods, String setMethodName, Object cl, Object cs, Object clTypeInfo, Object csTypeInfo) throws Exception {
        Class csTypeInfoClass = csTypeInfo.getClass();
        Class clTypeInfoClass = clTypeInfo.getClass();
        
        if(cl==null || cs==null) //value not found
            return;
        Object[] parameter = new Object[1];
        for (int i = 0; i < clMethods.length; i++) {
            String methodString = clMethods[i].getName();
            if(methodString.equals(setMethodName)) {
                Class[] parameterTypes = clMethods[i].getParameterTypes();
                Method setMethod = clTypeInfoClass.getMethod(setMethodName,parameterTypes);
                // have to handle here and not in translateCL because Enum
                if(setMethodName.indexOf("DistribClass")!=-1) {
                    Method methodCS = csTypeInfoClass.getMethod("getDistribClass",null);
                    Object csValue = methodCS.invoke(csTypeInfo,null);
                    String csName = csValue.getClass().getName();
                    //org.phycas.commandState.DistribClass.Enum csValue = (org.phycas.commandState.DistribClass.Enum)methodCS.invoke(csTypeInfo,null);
                    Method methodCL = clTypeInfoClass.getMethod("getDistribClass",null);
                    org.phycas.commandLanguage.DistribClass.Enum clValue = (org.phycas.commandLanguage.DistribClass.Enum)methodCL.invoke(clTypeInfo,null);
                    clValue = clValue.forString(csValue.toString());
                    Class[] clValueClass = {clValue.getClass()};
                    Method method = clTypeInfoClass.getMethod("setDistribClass",clValueClass);
                    Object o = method.invoke(clTypeInfo,new Object[] {clValue});
                }
                else if(!cl.getClass().getName().equals(cs.getClass().getName())) {
                    // command state and command language return types are different 
                    translateCL(cl,cs);
                }
                else { // command state and command language return types are the same
                    parameter[0] = cs;
                    Object o = setMethod.invoke(clTypeInfo, parameter);
                }
            }
        }
        
    }
    
    // necessary because structures of cs and cl schemas are different
    private void translateCL(Object cl, Object cs) throws Exception {
        String className = cl.getClass().getName();
        String csName = cs.getClass().getName();
        if(className.equals("org.phycas.commandLanguage.impl.DefStringValueElementImpl")) {
            Class c = cl.getClass();
            Class[] paramClass = {cs.getClass()};
            Method method = c.getMethod("setConstantVal",paramClass);
            Object o = method.invoke(cl,new Object[] {cs});
        }
        else if(className.equals("org.phycas.commandLanguage.impl.DefBoolValueElementImpl")) {
            Class c = cl.getClass();
            Class[] paramClass = {boolean.class};
            Method method = c.getMethod("setConstantVal",paramClass);
            Object o = method.invoke(cl,new Object[] {cs});
        }
        else if(className.equals("org.phycas.commandLanguage.impl.DefIntegerValueElementImpl") || className.equals("org.phycas.commandLanguage.impl.IntegerValueElementImpl")) {
            Class c = cl.getClass();
            java.math.BigInteger csBI = new java.math.BigInteger((String)cs);
            Class[] paramClass = {csBI.getClass()};
            Method method = c.getMethod("setConstantVal",paramClass);
            Object o = method.invoke(cl,new Object[] {csBI});
        }
        else if(className.equals("org.phycas.commandLanguage.impl.DefDoubleValueElementImpl") || className.equals("org.phycas.commandLanguage.impl.DoubleValueElementImpl")) {
            Class c = cl.getClass();
            Class[] paramClass = {cs.getClass()};
            Method method = c.getMethod("setConstantVal",paramClass);
            Object o = method.invoke(cl,new Object[] {cs});
        }
        else if(className.equals("org.phycas.commandLanguage.impl.StringValueListImpl")) {
            Class c = cl.getClass();
            Class[] paramClass = {cs.getClass()};
            Method method = c.getMethod("setConstantValArray",paramClass);
            Object o = method.invoke(cl,new Object[] {cs});
        }
        else if(className.equals("org.phycas.commandLanguage.impl.DistribRangeImpl")) {
            Class c = cl.getClass();
            Class c2 = cs.getClass();

            // constraint
            Method methodCS = c2.getMethod("getConstraint",null);
            org.phycas.commandState.DistribRangeEnum.Enum csValue = (org.phycas.commandState.DistribRangeEnum.Enum)methodCS.invoke(cs,null);
            Method methodCL = c.getMethod("getConstraint",null);
            org.phycas.commandLanguage.DistribRangeEnum.Enum clValue = (org.phycas.commandLanguage.DistribRangeEnum.Enum)methodCL.invoke(cl,null);
            clValue = clValue.forString(csValue.toString());
            Class[] clValueClass = {clValue.getClass()};
            Method method = c.getMethod("setConstraint",clValueClass);
            Object o = method.invoke(cl,new Object[] {clValue});
            
            // minVal
            methodCS = c2.getMethod("getMinVal",null);
            java.math.BigInteger csValue2 = (java.math.BigInteger)methodCS.invoke(cs,null);
            methodCL = c.getMethod("getMinVal",null);
            org.phycas.commandLanguage.IntegerValueElement clValue2 = (org.phycas.commandLanguage.IntegerValueElement)methodCL.invoke(cl,null);
            if(csValue2!=null) {
                if(clValue2==null) {
                    methodCL = c.getMethod("addNewMinVal",null);
                    clValue2 = (org.phycas.commandLanguage.IntegerValueElement)methodCL.invoke(cl,null);
                }
                clValue2.setConstantVal(csValue2);
            }
            else if(csValue2==null && clValue2!=null) { //unset
                method = c.getMethod("unsetMinVal",null);
                o = method.invoke(cl,null);
            }
            
            // maxVal
            methodCS = c2.getMethod("getMaxVal",null);
            csValue2 = (java.math.BigInteger)methodCS.invoke(cs,null);
            methodCL = c.getMethod("getMaxVal",null);
            clValue2 = (org.phycas.commandLanguage.IntegerValueElement)methodCL.invoke(cl,null);
            if(csValue2!=null) {
                if(clValue2==null) {
                    methodCL = c.getMethod("addNewMaxVal",null);
                    clValue2 = (org.phycas.commandLanguage.IntegerValueElement)methodCL.invoke(cl,null);
                }
                clValue2.setConstantVal(csValue2);
            }
            else if(csValue2==null && clValue2!=null) { //unset
                method = c.getMethod("unsetMaxVal",null);
                o = method.invoke(cl,null);
            }
        }
        else if (className.equals("org.phycas.commandLanguage.impl.OutputTypeInfoImpl$DefaultImpl")) {
            Class c = cl.getClass();
            Class c2 = cs.getClass();

            // suppress
            Method methodCS = c2.getMethod("getSuppress",null);
            boolean csValue = ((Boolean)methodCS.invoke(cs,null)).booleanValue();
            Class[] csValueClass = {boolean.class};
            Method method = c.getMethod("setSuppress",csValueClass);
            Object o = method.invoke(cl,new Object[] {new Boolean(csValue)});
            
            // file
            methodCS = c2.getMethod("getFile",null);
            //org.phycas.commandState.impl.OutputTypeInfoImpl.DefaultImpl.FileImpl csValue2 = (org.phycas.commandState.impl.OutputTypeInfoImpl.DefaultImpl.FileImpl)methodCS.invoke(cs,null);
            org.phycas.commandState.OutputTypeInfo.Default.File csValue2 = (org.phycas.commandState.OutputTypeInfo.Default.File)methodCS.invoke(cs,null);
            Method methodCL = c.getMethod("getFile",null);
            org.phycas.commandLanguage.OutputTypeInfo.Default.File clValue2 = (org.phycas.commandLanguage.OutputTypeInfo.Default.File)methodCL.invoke(cl,null);
            if(csValue2!=null && clValue2!=null) {
                clValue2.setAppend(csValue2.getAppend());
                clValue2.setReplace(csValue2.getReplace());
                clValue2.setPath(csValue2.getPath());
            }
            else if(clValue2==null && csValue2!=null) {
                Class[] classes = {org.phycas.commandLanguage.OutputTypeInfo.Default.File.class};
                method = c.getMethod("addNewFile",null);
                o = method.invoke(cl,null);
                methodCL = c.getMethod("getFile",null);
                clValue2 = (org.phycas.commandLanguage.OutputTypeInfo.Default.File)methodCL.invoke(cl,null);
                clValue2.setAppend(csValue2.getAppend());
                clValue2.setReplace(csValue2.getReplace());
                clValue2.setPath(csValue2.getPath());
            }
            else if(csValue2==null && clValue2!=null) {
                Class[] classes = {org.phycas.commandLanguage.OutputTypeInfo.Default.File.class};
                method = c.getMethod("unsetFile",null);
                o = method.invoke(cl,null);
            }
            
            // redirect
            methodCS = c2.getMethod("getRedirectArray",null);
            org.phycas.commandState.OutputRedirectionEnum.Enum[] csValue3 = (org.phycas.commandState.OutputRedirectionEnum.Enum[])methodCS.invoke(cs,null);
            methodCL = c.getMethod("getRedirectArray",null);
            org.phycas.commandLanguage.OutputRedirectionEnum.Enum[] clValue3 = (org.phycas.commandLanguage.OutputRedirectionEnum.Enum[])methodCL.invoke(cl,null);
            if(csValue3!=null && clValue3!=null) {
                for(int i=0;i<csValue3.length;i++) {
                    if(clValue3.length>=csValue3.length)
                        clValue3[i] = clValue3[i].forString(csValue3[i].toString());
                    else {
                        Class[] classes = {org.phycas.commandLanguage.OutputRedirectionEnum.Enum.class};
                        method = c.getMethod("addRedirect",classes);
                        o = method.invoke(cl,new Object[] {org.phycas.commandLanguage.OutputRedirectionEnum.Enum.forString(csValue3[i].toString())});
                    }
                }
            }
        }
        else {
            throw new Exception("Translation to set commandState value to commandLanguage not handled.");
        }
    }
    
    // get the value in the command language for the label
    private Object translateLabelToValue(String type, Object cs, String commandLabel) throws Exception {
        org.phycas.commandState.Command[] allCommands = commandState.getCommandArray();
        for(int i=0;i<allCommands.length;i++) {
            String currentLabel = allCommands[i].getLabel();
            if(currentLabel.equals(commandLabel)) {
                // find method to get Handler class for command
                Class mainClass = main.getClass();
                String methodName = "get"+commandLabel+"Handler";
                Method method = null;
                try {
                    method = mainClass.getMethod(methodName,null);
                }
                catch (NoSuchMethodException e) {
                    continue; // skip if command does not exist
                }
                Object handlerInstance = method.invoke(main, null);
                Class handlerClass = handlerInstance.getClass();
                // find appropriate TypeInfo in command language CmdParam
                Method getParams = handlerClass.getMethod("getAllParams", null);
                org.phycas.commandLanguage.CmdParam[] clCmdParams = (org.phycas.commandLanguage.CmdParam[])getParams.invoke(handlerInstance,null);
                for(int l=0;l<clCmdParams.length;l++) {
                    String label = clCmdParams[l].getLabel();
                    if(label.equals(cs.toString())) {
                        Class paramClass = clCmdParams[l].getClass();
                        Method[] allMethods = paramClass.getDeclaredMethods();
                        for(int k=0;k<allMethods.length;k++) {
                            String paramMethodName = allMethods[k].getName();
                            if(paramMethodName.equals("getIntegerTypeInfo") || paramMethodName.equals("getDoubleTypeInfo")) {
                                Object clTI = allMethods[k].invoke(clCmdParams[l], null);
                                if(clTI!=null) {
                                    Method[] typeInfoMethods = clTI.getClass().getDeclaredMethods();
                                    for(int m=0;m<typeInfoMethods.length;m++) {
                                        String name = typeInfoMethods[m].getName();
                                        if(name.equals("get"+type+"Val")) {
                                            Object value = typeInfoMethods[m].invoke(clTI,null);
                                            String valueClass = value.getClass().getName();
                                            String v = value.toString();
                                            // get all Cmd Params in Command
                                            org.phycas.commandState.CmdParam[] allParams = allCommands[i].getCmdParamArray();
                                            for(int j=0;j<allParams.length;j++) {
                                                String paramLabel = allParams[j].getLabel();
                                                if(cs.toString().equals(paramLabel)) {
                                                    return value;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        
        return cs;
    }
    
    // handle sets
    private void setManager(SetManager setMgr, String type) {
        PhycasSetManager psm = null;
        if(type.equals("char"))
            psm = main.charSetMgr;
        else if(type.equals("tax"))
            psm = main.taxSetMgr;
        else if(type.equals("tree"))
            psm = main.treeSetMgr;
        if(psm==null)
            return;
        
        SetManager.IndexLabel[] indexLabels = setMgr.getIndexLabelArray();
        if(indexLabels!=null) {
            for(int i=0;i<indexLabels.length;i++) {
                String l = indexLabels[i].getLabel();
                int in = indexLabels[i].getIndex().intValue();
                psm.addIndexLabel(in,l);
            }
        }
        
        PhycasSet[] all = psm.getAllSets();
        KnownSet[] allSets = setMgr.getKnownSetArray();
        if(allSets!=null){
            for(int i=0;i<allSets.length;i++) {
                String label = allSets[i].getLabel();
                boolean exists = false;
                if(all!=null) {
                    for(int j=0;j<all.length;j++) {
                        if(all[j].getLabel().equals(label)) {
                            exists = true;
                            all[j].setMembers(allSets[i].getMembers());
                            break;
                        }
                    }
                }
                if(!exists) {
                    PhycasSet s = new PhycasSet(label);
                    s.setMembers(allSets[i].getMembers());
                    psm.addSet(s);
                }
            }
        }
        
        java.math.BigInteger maxIndex = setMgr.getMaxIndex();
        if(maxIndex!=null)
            psm.setMaxIndex(maxIndex.intValue());
    }
    
    // disable a command param
    private void disableParam(String commandLabel, String param, PhycasMain main, Class mainClass) throws Exception {
        Field[] allFields = mainClass.getFields();
        for (int i=0;i<allFields.length;i++) {
            String fieldName = allFields[i].getName();
            Class fieldClass = allFields[i].getType();
            if(fieldName.indexOf("cmd"+commandLabel)!=-1) {
                if(fieldName.indexOf(param)!=-1 && fieldName.indexOf("Action")==-1) {
                    Object fieldValue = allFields[i].get(main);
                    Method method = fieldClass.getMethod("setEnabled",new Class[] {boolean.class});
                    Boolean t = new Boolean(false);
                    Object o = method.invoke(fieldValue,new Object[] {t});
                }
            }
        }
    }
}
