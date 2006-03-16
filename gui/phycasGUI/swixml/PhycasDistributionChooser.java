package phycasGUI.swixml;

import javax.swing.*;
import java.awt.event.*;
import java.awt.*;
import java.util.Hashtable;
import org.phycas.commandLanguage.*;

/**
 * Dialog box for choosing a probability distribution type.
 * @author Vanessa Jackson
 */
public class PhycasDistributionChooser extends JDialog {
    
    // Need dynamic way to know which distribution types are available
    
    // info,allowedvalues,categories,ranges [i] corresponds to names[i]
    private static final String[] names = {"Bernoulli", "Binomial", "Exponential",
    "Gamma", "InvGamma", "Dirichlet", "Beta", "Uniform"};
    private static final String[][] info = {
        {"Probability", "OR", "Mean"},
        {"Probability", "AND", "NTrials"},
        {"Hazard", "OR", "Mean"},
        {"Shape", "AND", "Scale", "OR", "Mean", "AND", "Var", "OR", "Mean", "AND", "StdDev"},
        {"Shape", "AND", "Scale", "OR", "Mean", "AND", "Var", "OR", "Mean", "AND", "StdDev"},
        {"Values"},
        {"Alpha", "AND", "Beta", "OR", "Mean", "AND", "StdDev"},
        {"Lower", "AND", "Upper"}
    };
    // name, paramLabel, min, max
    private static final String[][] allowedValues = {
        {"Bernoulli", "Probability", "0", "1"},
        {"Bernoulli", "Mean", "0", "1"},
        {"Binomial", "Probability", "0", "1"},
        {"Binomial", "NTrials", "0", "Integer.MAX_VALUE"},
        {"Exponential", "Hazard", "0", "Double.MAX_VALUE"},
        {"Exponential", "Mean", "0", "Double.MAX_VALUE"},
        {"Gamma", "Shape", "0", "Double.MAX_VALUE"},
        {"Gamma", "Scale", "0", "Double.MAX_VALUE"},
        {"Gamma", "Mean", "0", "Double.MAX_VALUE"},
        {"Gamma", "Var", "0", "Double.MAX_VALUE"},
        {"Gamma", "StdDev", "0", "Double.MAX_VALUE"},
        {"InvGamma", "Shape", "0", "Double.MAX_VALUE"},
        {"InvGamma", "Scale", "0", "Double.MAX_VALUE"},
        {"InvGamma", "Mean", "0", "Double.MAX_VALUE"},
        {"InvGamma", "Var", "0", "Double.MAX_VALUE"},
        {"InvGamma", "StdDev", "0", "Double.MAX_VALUE"},
        {"Dirichlet", "Values", "0", "Double.MAX_VALUE"},
        {"Beta", "Alpha", "0", "Double.MAX_VALUE"},
        {"Beta", "Beta", "0", "Double.MAX_VALUE"},
        {"Beta", "Mean", "0", "Double.MAX_VALUE"},
        {"Beta", "StdDev", "0", "Double.MAX_VALUE"},
        {"Uniform", "Lower", "Double.MIN_VALUE", "Double.MAX_VALUE"},
        {"Uniform", "Upper", "Double.MIN_VALUE", "Double.MAX_VALUE"},
    };
    private static final String[] categories = {"Discrete", "Discrete", "Continuous",
    "Continuous", "Continuous", "Continuous", "Continuous", "Continuous"};
    private static final String[] ranges = {"Bounded", "Bounded", "NonNegative",
    "NonNegative", "NonNegative", "SumToOne", "SumToOne", "Bounded"};
    
    private String selectedName = "";
    private boolean madeSelection = false;
    private Hashtable textFields = null;
    private JRadioButton[] radioButtons = null;
    private StringBuffer errorMsgs = new StringBuffer();
    private DistributionFocusListener focusListener = null;
    private String category = "";
    private int numVariates = -1;
    private String rangeConstraint = "";
    private int rangeMin = -1;
    private int rangeMax = -1;
    
    /**
     * Creates a modal Dialog to choose a distribution
     * @param owner owner of the dialog box
     * @param allParams command param data from XMLBeans
     */
    public PhycasDistributionChooser(JFrame owner, CmdParam[] allParams) {
        super(owner,"Choose a Distribution",true);
        
        radioButtons = new JRadioButton[names.length];
        final ButtonGroup group = new ButtonGroup();
        
        for(int i=0;i<names.length;i++) {
            JRadioButton rb = new JRadioButton(names[i]);
            rb.setPreferredSize(new Dimension(100,20));
            rb.setMaximumSize(new Dimension(100,20));
            radioButtons[i] = rb;
            radioButtons[i].setActionCommand(names[i]);
        }
        
        for (int i = 0; i < radioButtons.length; i++) {
            group.add(radioButtons[i]);
        }
        focusListener = new DistributionFocusListener(radioButtons);
        
        radioButtons[0].setSelected(true);
        
        textFields = new Hashtable();
        for (int i=0;i<names.length;i++) {
            Object[] fields = new Object[info[i].length*2];
            int count = 0;
            for (int j=0;j<info[i].length;j++) {
                if(info[i][j].equals("AND") || info[i][j].equals("OR")) {
                    JLabel l = new JLabel(info[i][j]);
                    l.setPreferredSize(new Dimension(30,20));
                    l.setMaximumSize(new Dimension(30,20));
                    fields[count] = l;
                    count++;
                }
                else {
                    JLabel l = new JLabel(info[i][j]);
                    l.setPreferredSize(new Dimension(50,20));
                    l.setMaximumSize(new Dimension(50,20));
                    fields[count] = l;
                    JTextField tf = new JTextField();
                    tf.setPreferredSize(new Dimension(50,20));
                    tf.setMaximumSize(new Dimension(50,20));
                    tf.addFocusListener(focusListener);
                    tf.setName(names[i]+"_"+info[i][j]);
                    fields[count+1] = tf;
                    count+=2;
                }
            }
            textFields.put(names[i],fields);
        }
        
        JButton selectButton = new JButton("Select");
        selectButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                madeSelection = true;
                selectedName = group.getSelection().getActionCommand();
                if(validate(selectedName))
                    hide();
                else
                    madeSelection = false;
            }
        });
        
        JPanel box = new JPanel();
        box.setLayout(new BoxLayout(box, BoxLayout.PAGE_AXIS));
        for (int i = 0; i < radioButtons.length; i++) {
            JPanel line = new JPanel();
            line.setLayout(new FlowLayout(FlowLayout.LEFT,15, 5));
            line.add(radioButtons[i]);
            Object[] fields = (Object[])textFields.get(names[i]);
            for (int j=0;j<fields.length;j++) {
                JComponent jc = (JComponent)fields[j];
                if(jc!=null)
                    line.add(jc);
            }
            box.add(line);
        }
        
        JPanel bottom = new JPanel();
        bottom.add(selectButton);
        
        JPanel pane = new JPanel(new BorderLayout());
        pane.add(box, BorderLayout.PAGE_START);
        pane.add(bottom, BorderLayout.PAGE_END);
        
        getContentPane().add(pane);
        
        getDistributionRestrictions(allParams);
        disableOptions();
    }
    
    /**
     * Get the user's chosen distribution
     * @return distribution
     */    
    public String getChoice() {
        if(madeSelection)
            return format(selectedName);
        else
            return "";
    }
    
    // get the restrictions on the distribution type available from the XMLBeans
    private void getDistributionRestrictions(CmdParam[] allParams) {
        for(int i=0;i<allParams.length;i++) {
            DistributionTypeInfo dti = allParams[i].getDistributionTypeInfo();
            if(dti!=null) {
                numVariates = dti.getNumVariates().intValue();
                DistribRange dr = dti.getRangeConstraint();
                IntegerValueElement min = dr.getMinVal();
                if(min!=null)
                    rangeMin = min.getConstantVal().intValue();
                IntegerValueElement max = dr.getMaxVal();
                if(max!=null)
                    rangeMax = max.getConstantVal().intValue();
                rangeConstraint = dr.getConstraint().toString();
                category = dti.getDistribClass().toString();
                System.out.println("DistribRest = "+category+" "+rangeMin+" "+rangeMax+" "+rangeConstraint+" "+numVariates);
            }
        }
    }
    
    private void disableOptions() {
        for(int i=0;i<names.length;i++) {
            boolean enable = true;
            if(!category.equals("Any") && !categories[i].equals(category)) {
                enable = false;
            }
            else if (!rangeConstraint.equals("Any") && !ranges[i].equals(rangeConstraint)) {
                enable = false;
            }
            Object[] fields = (Object[])textFields.get(names[i]);
            for(int j=0;j<fields.length;j++) {
                if(fields[j]!=null)
                    ((JComponent)fields[j]).setEnabled(enable);
            }
            for(int j=0;j<radioButtons.length;j++) {
                String name = ((JRadioButton)radioButtons[i]).getText();
                if(name.equals(names[i]))
                    ((JRadioButton)radioButtons[i]).setEnabled(enable);
            }
        }
    }
    
    private boolean validate(String name) {
        boolean success = true;
        Object[] fields = (Object[])textFields.get(name);
        String value1 = "";
        String value2 = "";
        String label1 = "";
        String label2 = "";
        for (int i=0;i<fields.length;i++) {
            if(fields[i]!=null) {
                JLabel l = (JLabel)fields[i];
                String label = l.getText();
                if(label.equals("OR")) { // make sure only one filled in
                    if(value1.equals("") && value2.equals(""))
                        continue;
                    l = (JLabel)fields[i+1];
                    label = l.getText();
                    JTextField tf = (JTextField)fields[i+2];
                    String value = tf.getText();
                    value2 = value;
                    label2 = label;
                    if(!value1.equals("") && !value2.equals("")) {
                        errorMsgs.append("Please choose either field for the "+name+" distribution, but not both.\n");
                        success = false;
                        break;
                    }
                    i+=2;
                    continue;
                }
                else if(label.equals("AND")) { // validate both together
                    l = (JLabel)fields[i+1];
                    label = l.getText();
                    JTextField tf = (JTextField)fields[i+2];
                    String value = tf.getText();
                    value2 = value;
                    label2 = label;
                    if(!value1.equals("") && !value1.equals("")) {
                        boolean success1 = specialValidation(selectedName,label1,value1);
                        boolean success2 = specialValidation(selectedName,label2,value2);
                        if(!success1 || !success2) {
                            success = false;
                            break;
                        }
                    }
                    label1 = ""; value1 = "";
                    label2 = ""; value2 = "";
                    i+=2;
                    continue;
                }
                else {
                    if(!label1.equals("") && !value1.equals("")) {
                        success = specialValidation(selectedName,label1,value1);
                        label1 = "";
                        value1 = "";
                    }
                }
                
                JTextField tf = (JTextField)fields[i+1];
                String value = tf.getText();
                double val = -1;
                if(!value.equals("")) {
                    try {
                        val = new Double(value).doubleValue();
                    }
                    catch (Exception e) {
                        errorMsgs.append(label+" must be a number.\n");
                        success = false;
                    }
                }
                value1 = value;
                label1 = label;
                i++;
            }
            else if(success) {
                if(!label1.equals("") && !value1.equals("")) {
                    success = specialValidation(selectedName,label1,value1);
                    label1 = "";
                    value1 = "";
                }
            }
            
        }
        if(!success) {
            JOptionPane.showMessageDialog(this,errorMsgs.toString(),"Error",JOptionPane.ERROR_MESSAGE);
        }
        errorMsgs = new StringBuffer();
        return success;
    }
    
    private boolean specialValidation(String name, String label, String value) {
        boolean success = true;
        for(int i=0;i<allowedValues.length;i++) {
            String[] vals = (String[])allowedValues[i];
            int type = -1; // 0 for int, 1 for double
            if(vals[0].equals(name)) {
                if(vals[1].equals(label)) {
                    String minS = vals[2];
                    if(minS.indexOf("Integer")!=-1) type = 0;
                    else type = 1;
                    String maxS = vals[3];
                    if(maxS.indexOf("Integer")!=-1) type = 0;
                    else type = 1;
                    if(type==0) {
                        int min,max;
                        if(minS.indexOf("MIN_VAL")!=-1)
                            min = Integer.MIN_VALUE;
                        else
                            min = new Integer(minS).intValue();
                        if(maxS.indexOf("MAX_VAL")!=-1)
                            max = Integer.MAX_VALUE;
                        else
                            max = new Integer(maxS).intValue();
                        try {
                            int v = new Integer(value).intValue();
                            if(v<min || v>max) {
                                errorMsgs.append(label+" must be greater than or equal to "+min+" and less than or equal to "+max+".\n");
                                success = false;
                            }
                        }
                        catch(NumberFormatException nfe) {
                            errorMsgs.append(label+" must be an integer greater than or equal to "+min+" and less than or equal to "+max+".\n");
                            success = false;
                        }
                    }
                    else {
                        double min,max;
                        if(minS.indexOf("MIN_VAL")!=-1)
                            min = Long.MIN_VALUE;
                        else
                            min = new Double(minS).doubleValue();
                        if(maxS.indexOf("MAX_VAL")!=-1)
                            max = Double.MAX_VALUE;
                        else
                            max = new Double(maxS).doubleValue();
                        try {
                            double v = new Double(value).doubleValue();
                            if(v<min || v>max) {
                                errorMsgs.append(label+" must be greater than or equal to "+min+" and less than or equal to "+max+".\n");
                                success = false;
                            }
                        }
                        catch(NumberFormatException nfe) {
                            errorMsgs.append(label+" must be a number greater than or equal to "+min+" and less than or equal to "+max+".\n");
                            success = false;
                        }
                    }
                }
            }
        }
        return success;
    }
    
    private String format(String name) {
        StringBuffer sb = new StringBuffer(name);
        sb.append("(");
        if(textFields==null)
            return "";
        Object[] fields = (Object[])textFields.get(name);
        if(fields==null)
            return "";
        for (int i=0;i<fields.length;i++) {
            if(fields[i]!=null) {
                JLabel l = (JLabel)fields[i];
                String label = l.getText();
                if(label.equals("AND") || label.equals("OR")) {
                    continue;
                }
                JTextField tf = (JTextField)fields[i+1];
                String value = tf.getText();
                if(!value.equals(""))
                    sb.append(label+"="+value+" ");
                i++;
            }
        }
        sb.append(")");
        String r = sb.toString();
        return r;
    }
    
    private class DistributionFocusListener implements FocusListener {
        JRadioButton[] radioButtons = null;
        
        public DistributionFocusListener(JRadioButton[] rbs) {
            radioButtons = rbs;
        }
        
        public void focusGained(java.awt.event.FocusEvent focusEvent) {
            String fieldName = focusEvent.getComponent().getName();
            if(fieldName!=null) {
                System.out.println(fieldName);
                String distribName = fieldName.substring(0,fieldName.indexOf('_'));
                System.out.println(distribName);
                for(int i=0;i<radioButtons.length;i++) {
                    if(radioButtons[i].getText().equals(distribName)) {
                        radioButtons[i].setSelected(true);
                    }
                }
            }
        }
        
        public void focusLost(java.awt.event.FocusEvent focusEvent) {
        }
        
    }
    
}

