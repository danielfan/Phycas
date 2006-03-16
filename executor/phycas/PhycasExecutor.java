
package phycas;

import org.swixml.SwingEngine;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

/**
 * The main class of an application used to start the two components
 * of Phycas - the Main cpp program and the Java Swing GUI.
 * @author Vanessa Jackson
 */
public class PhycasExecutor {
    
    // these are required to be declared here and public since using Swixml
    public JFrame frame;
    public JLabel lbl_BaseDir;
    public JTextField tf_BaseDir;
    public JLabel lbl_Host;
    public JTextField tf_Host;
    public JLabel lbl_Port;
    public JTextField tf_Port;
    public JLabel lbl_Location;
    public JRadioButton rb_Host_Local;
    public JRadioButton rb_Host_Remote;
    public JTextArea taMsgs;
    
    private SwingEngine swix = null;
    public JFileChooser fc; // file selection dialog
    public PhycasExecutorHandler handler = new PhycasExecutorHandler(this);
    public Action okAction = handler;
    public Action exitAction = handler;
    public Action action_ChooseBaseDir = handler;
    public Action action_Local = handler;
    public Action action_Remote = handler;
    
    
    /**
     * Creates a new instance of PhycasExecutor
     * @throws Exception any Exception
     */
    public PhycasExecutor() throws Exception {
        swix = new SwingEngine(this);
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        fc = new JFileChooser();
        fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    }
    
    /**
     * The main method
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            PhycasExecutor pe = new PhycasExecutor();
            pe.display();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Render the Swixml file and Display the window.
     */    
    public void display() {
        try {
            System.out.println("Executor: running");
            Container cFrame = swix.render("xml/PhycasExecutorMain.xml");
            boolean gotProps = getProperties();
            frame.setDefaultLookAndFeelDecorated(true);
            frame.setLocationRelativeTo(null);
            if(gotProps)
                tf_BaseDir.setText(System.getProperty("phycas.home",""));
            frame.addWindowListener(handler);
            tf_BaseDir.getInputMap().put(KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, 0),"okAction");
            tf_BaseDir.getActionMap().put("okAction",handler.okAction);
            tf_Host.getInputMap().put(KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, 0),"okAction");
            tf_Host.getActionMap().put("okAction",handler.okAction);
            tf_Port.getInputMap().put(KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, 0),"okAction");
            tf_Port.getActionMap().put("okAction",handler.okAction);
            frame.setVisible(true);
            tf_BaseDir.requestFocus();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Get phycas properties from file.
     * @return <CODE>true</CODE> if successfully retrieved properties from file
     */    
    public boolean getProperties() {
        try {
            String userHome = System.getProperty("user.home");
            String fileSep = System.getProperty("file.separator");
            String filename = userHome+fileSep+"phycas"+fileSep+"phycas.properties";
            FileInputStream propFile = new FileInputStream(filename);
            Properties p = new Properties(System.getProperties());
            p.load(propFile);
            System.setProperties(p);
            propFile.close();
        }
        catch(Exception e) {
            taMsgs.append("Unable to load properties file.\n");
            taMsgs.setCaretPosition(taMsgs.getDocument().getLength());
            System.out.println("Unable to load properties file: "+e.getMessage());
            return false;
        }
        return true;
    }
    
    /**
     * Save phycas properties to file.
     */    
    public void saveProperties() {
        try {
            String userHome = System.getProperty("user.home");
            String fileSep = System.getProperty("file.separator");
            // check if user.home/phycas dir exists
            File dir = new File(userHome+fileSep+"phycas");
            if(!dir.exists())
                dir.mkdir();
            String filename = userHome+fileSep+"phycas"+fileSep+"phycas.properties";
            FileOutputStream propFile = new FileOutputStream(filename);
            Properties p = System.getProperties();
            Properties phycasP = new Properties();
            Enumeration keys = p.keys();
            while(keys.hasMoreElements()) {
                String key = (String)keys.nextElement();
                if(key.startsWith("phycas"))
                    phycasP.setProperty(key, p.getProperty(key));
            }
            phycasP.store(propFile,null);
            propFile.close();
        }
        catch(Exception e) {
            taMsgs.append("Unable to write properties file.\n");
            taMsgs.setCaretPosition(taMsgs.getDocument().getLength());
            System.out.println("Unable to write properties file: "+e.getMessage());
        }
    }
    
}
