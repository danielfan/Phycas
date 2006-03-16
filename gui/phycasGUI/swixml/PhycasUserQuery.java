package phycasGUI.swixml;

import org.jdom.*;
import java.util.*;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JDialog;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeEvent;

/**
 * Handle user_query messages from the main program to the gui
 * @author Vanessa Jackson
 */
public class PhycasUserQuery {
    
    private JFrame frame = null;
    private String xml_title = null;
    private String xml_message = null;
    private String xml_choice = null;
    private ArrayList choices = null;
    private String xml_attribType = null;
    
    private  static final String TITLE = "title";
    private  static final String MSG = "message";
    private  static final String TYPE = "type";
    private  static final String CHOICE = "choice";
    private  static final String DEFAULT = "default";
    
    private  static final String FILE = "file";
    private  static final String STRING = "string";
    private  static final String ALERT = "alert";
    private  static final String CHOICES = "choices";
    private  static final String C_OK = "cancel_ok";
    private  static final String N_Y = "no_yes";
    
    private int defaultChoice = 0;
    
    /**
     * Creates a new instance of PhycasUserQuery
     * @param f main frame of Phycas gui
     */
    public PhycasUserQuery(JFrame f) {
        frame = f;
    }
    
    // parse user_query xml message
    private int assignChildren(Element current) {
        
        //System.out.println("UQ element name = "+current.getName());
        //System.out.println("UQ element text = "+current.getText());
        Attribute att = current.getAttribute(TYPE);
        System.out.println("UQ TYPE = "+att.getValue());
        xml_attribType = att.getValue();
        List children = current.getChildren();
        Iterator iterator = children.iterator();
        int numChoice = 0;
        while (iterator.hasNext()) {
            Element child = (Element) iterator.next();
            System.out.println("UQ CHILD = "+child.getName());
            if(child.getName().equals(TITLE))
                xml_title = child.getText().trim();
            else if(child.getName().equals(MSG))
                xml_message = child.getText().trim();
            else if(child.getName().equals(CHOICE)) {
                Attribute childAtt = child.getAttribute(DEFAULT);
                if(childAtt != null) {
                    String def = childAtt.getValue();
                    if(def.equals("true"))
                        defaultChoice = numChoice;
                }
                xml_choice = child.getText().trim();
                numChoice++;
            }
            else {
                System.out.println("INVALID XML TAG");
                return 1;
            }
            if(xml_choice!=null) {
                if(choices==null)
                    choices = new ArrayList();
                choices.add(xml_choice);
            }
        }
        return 0;
    }
    
    /**
     * Show dialog box with message
     * @param uq main element of xml message
     * @return user's answer to message
     */
    public String showDialog(Element uq) {
        assignChildren(uq);
        final JDialog dialog = new JDialog(frame,xml_title,true);
        dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
        if(xml_attribType.equals(FILE)) {
            //option pane with button to popup file chooser
            // VKJ 11/04 - easier to just refer user back to main window
            /*String s = (String)JOptionPane.showInputDialog(
            frame,
            xml_message,
            xml_title,
            JOptionPane.PLAIN_MESSAGE,
            null,
            null,
            "");*/
            JOptionPane.showMessageDialog(frame,
            xml_message,
            xml_title,
            JOptionPane.PLAIN_MESSAGE);
            System.out.println("RESPONSE = Cancel");
            return "cancel";
        }
        else if(xml_attribType.equals(STRING)) {
            // option pane with entry field
            String s = (String)JOptionPane.showInputDialog(
            frame,
            xml_message,
            xml_title,
            JOptionPane.PLAIN_MESSAGE,
            null,
            null,
            "");
            System.out.println("RESPONSE = "+s);
            if(s==null)
                return "";
            return s;
        }
        else if(xml_attribType.equals(ALERT)) {
            // option pane with OK button
            JOptionPane.showMessageDialog(frame,
            xml_message,
            xml_title,
            JOptionPane.PLAIN_MESSAGE);
            return "1";
        }
        else if(xml_attribType.equals(CHOICES)) {
            // option pane with combo box
            Object[] c = choices.toArray();
            String s = (String)JOptionPane.showInputDialog(
            frame,
            xml_message,
            xml_title,
            JOptionPane.PLAIN_MESSAGE,
            null,
            c,
            c[defaultChoice]);
            int n = -1;
            if(s==null) {
                // look for option that says cancel
                for(int i=0;i<c.length;i++) {
                    String lower = ((String)c[i]).toLowerCase();
                    if(lower.indexOf("cancel")!=-1)
                        n = i;
                }
            }
            else {
                // get int for chosen string
                for(int i=0;i<c.length;i++) {
                    if(((String)c[i]).equals(s))
                        n = i;
                }
            }
            System.out.println("RESPONSE = "+n);
            if(n!=-1)
                return Integer.toString(n);
            return "";
        }
        else if(xml_attribType.equals(C_OK)) {
            // option pane with ok and cancel buttons
            int n = JOptionPane.showConfirmDialog(
            frame,
            xml_message,
            xml_title,
            JOptionPane.OK_CANCEL_OPTION);
            System.out.println("RESPONSE = "+n);
            return Integer.toString(n);
        }
        else if(xml_attribType.equals(N_Y)) {
            // option pane with yes and no buttons
            int n = JOptionPane.showConfirmDialog(
            frame,
            xml_message,
            xml_title,
            JOptionPane.YES_NO_OPTION);
            System.out.println("RESPONSE = "+n);
            return Integer.toString(n);
        }
        else
            return null;
    }
    
    
}
