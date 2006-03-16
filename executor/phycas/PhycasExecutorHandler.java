
package phycas;

import javax.swing.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;
import java.net.*;

/**
 * Handles the execution of the Main and GUI programs.
 * @author Vanessa Jackson
 */
public class PhycasExecutorHandler extends AbstractAction implements WindowListener {
    
    private PhycasExecutor main = null;
    private Process cppProcess = null;
    private Process guiProcess = null;
    private int port = -1;
    private String host = null;
    private boolean startPhycasProcessor = true;
    
    /**
     * Creates a new instance of PhycasExecutorHandler
     * @param m instance of PhycasExecutor which contains the window and it's components.
     */
    public PhycasExecutorHandler(PhycasExecutor m) {
        main = m;
    }
    
    /**
     * Handles any action on the PhycasExecutor window.
     * @param e event on window
     */
    public void actionPerformed(ActionEvent e) {
        
        // Exit menu item
        if(e.getActionCommand().equals("exitAction"))
            end();
        
        // Start Phycas button
        if(e.getActionCommand().equals("okAction"))
            handleStartButton();
        
        // Choose Base Directory button
        if(e.getActionCommand().equals("action_ChooseBaseDir"))
            chooseBaseDir();
        
        // Local machine radio button
        if(e.getActionCommand().equals("action_Local"))
            main.tf_Host.setEnabled(false);
        
        // Remote machine radio button
        if(e.getActionCommand().equals("action_Remote"))
            main.tf_Host.setEnabled(true);
        
    }
    
    /**
     * Validate the fields in the window.
     * @return <CODE>true</CODE> if all fields valid
     */
    private boolean validateAllFields() {
        boolean success = true;
        if(main.tf_BaseDir.getText().trim().equals("")) {
            main.tf_BaseDir.setText("");
            main.lbl_BaseDir.setForeground(java.awt.Color.RED);
            main.taMsgs.append("Please enter a Phycas Installation Directory.\n");
            main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
            success = false;
        }
        else
            main.lbl_BaseDir.setForeground(java.awt.Color.BLACK);
        
        if(main.rb_Host_Local.isSelected()) {
            startPhycasProcessor = true;
            String p = main.tf_Port.getText().trim();
            if(!p.equals("")) {
                try {
                    port = new Integer(p).intValue();
                    main.lbl_Port.setForeground(java.awt.Color.BLACK);
                }
                catch (NumberFormatException e) {
                    main.tf_Port.setText("");
                    main.lbl_Port.setForeground(java.awt.Color.RED);
                    main.taMsgs.append("Port must be a number.\n");
                    main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                    success = false;
                }
            }
            else
                main.lbl_Host.setForeground(java.awt.Color.BLACK);
        }
        else if(main.rb_Host_Remote.isSelected()) {
            startPhycasProcessor = false;
            host = main.tf_Host.getText().trim();
            if(host.equals("")) {
                main.tf_Host.setText("");
                main.lbl_Host.setForeground(java.awt.Color.RED);
                main.taMsgs.append("Please enter a Host ip address.\n");
                main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                success = false;
            }
            else
                main.lbl_Host.setForeground(java.awt.Color.BLACK);
            String p = main.tf_Port.getText().trim();
            if(p.equals("")) {
                main.tf_Port.setText("");
                main.lbl_Port.setForeground(java.awt.Color.RED);
                main.taMsgs.append("Please enter a Port number.\n");
                main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                success = false;
            }
            else {
                try {
                    port = new Integer(p).intValue();
                    main.lbl_Port.setForeground(java.awt.Color.BLACK);
                }
                catch (NumberFormatException e) {
                    main.tf_Port.setText("");
                    main.lbl_Port.setForeground(java.awt.Color.RED);
                    main.taMsgs.append("Port must be a number.\n");
                    main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                    success = false;
                }
            }
        }
        if(port==0) {
            main.tf_Port.setText("");
            main.lbl_Port.setForeground(java.awt.Color.RED);
            main.taMsgs.append("Port is invalid.\n");
            main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
            success = false;
        }
        return success;
    }
    
    /**
     * Shows file chooser for Phycas directory
     */
    private void chooseBaseDir() {
        System.out.println("Choose");
        int retVal = main.fc.showOpenDialog(main.tf_BaseDir);
        if (retVal == JFileChooser.APPROVE_OPTION) {
            File file = main.fc.getSelectedFile();
            main.tf_BaseDir.setText(file.getAbsolutePath());
            main.fc.setSelectedFile(new File(""));
        }
        
    }
    
    /**
     * Clean up processes
     */
    private void end() {
        killProcesses();
        System.exit( 0 );
    }
    
    /**
     * Main method to create processes
     */
    private void begin() {
        try {
            System.out.println("Start");
            // first check for an open port, if one has not been found previously
            if(port==-1) {
                ServerSocket s = null;
                try {
                    s = new ServerSocket(0); // use 0 for any free port
                    System.out.println("Found available port "+s.getLocalPort());
                    port = s.getLocalPort();
                }
                catch (Exception e) {
                    System.out.println("Unable to find available port.");
                    e.printStackTrace();
                    main.taMsgs.append("Unable to find available port.\n");
                    main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                    return;
                }
                try {
                    if(s!=null)
                        s.close();
                }
                catch (Exception e) {
                    System.out.println("Unable to close socket.");
                    e.printStackTrace();
                }
            }
            
            if(port==-1) {
                main.taMsgs.append("Unable to find available port.\n");
                main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                return;
            }
            
            // set necessary properties
            System.setProperty("phycas.home",main.tf_BaseDir.getText());
            String fileSep = System.getProperty("file.separator");
            String javaHome = System.getProperty("java.home");
            System.out.println("JavaHome = "+javaHome);
            String jreHome = System.getProperty("phycas.jre_home",javaHome);
            System.out.println("JRE = "+jreHome);
            String phycasHome = System.getProperty("phycas.home");
            System.out.println("Phycas = "+phycasHome);
            
            if(startPhycasProcessor) {
                main.taMsgs.append("Starting main Phycas process\n");
                main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                String cpp = phycasHome+fileSep+"PhycasSocket.exe "+port;
                System.out.println("Executing "+cpp);
                if(cppProcess==null)
                    cppProcess = Runtime.getRuntime().exec(cpp);
                else {
                    main.taMsgs.append("Phycas main process already running.\n");
                    main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                }
                PhycasReader cppReader = new PhycasReader(main,cppProcess.getInputStream(),0,0);
                cppReader.start();
                PhycasReader cppReaderErr = new PhycasReader(main,cppProcess.getErrorStream(),1,1);
                cppReaderErr.start();
            }
            
            if(host==null || host.equals(""))
                host = "localhost";
            // If there are spaces in the phycasHome directory, the program will not run
            // cannot put "" around the path because unix will not accept
            String javaExec = fileSep+"bin"+fileSep+"java -jar \""+phycasHome+fileSep+"PhycasGui.jar\" "+host+" "+port;
            
            main.taMsgs.append("Starting Phycas user interface\n");
            main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
            boolean javaHomeSuccess = false;
            try { //try javaHome first
                String gui = javaHome+javaExec;
                System.out.println("Executing "+gui);
                if(guiProcess==null)
                    guiProcess = Runtime.getRuntime().exec(gui);
                else {
                    main.taMsgs.append("Phycas user interface already running.\n");
                    main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                }
                System.setProperty("phycas.jre_home",javaHome);
                javaHomeSuccess = true;
            }
            catch(Exception e) {
                System.out.println("Java Home didn't work, will try jreHome");
            }
            if(!javaHomeSuccess) {
                try { //try jreHome next
                    String gui = jreHome+javaExec;
                    System.out.println("Executing "+gui);
                    if(guiProcess==null)
                        guiProcess = Runtime.getRuntime().exec(gui);
                    else {
                        main.taMsgs.append("Phycas user interface already running.\n");
                        main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                    }
                    System.setProperty("phycas.jre_home",jreHome);
                }
                catch(Exception e) {
                    System.out.println("Throwing exception");
                    throw(e);
                }
            }
            PhycasReader guiReader = new PhycasReader(main,guiProcess.getInputStream(),1,0);
            guiReader.start();
            PhycasReader guiReaderErr = new PhycasReader(main,guiProcess.getErrorStream(),1,1);
            guiReaderErr.start();
            
            new PhycasExitListener(main,cppProcess,0).start();
            new PhycasExitListener(main,guiProcess,1).start();
            
            main.saveProperties();
        }
        catch (Exception e) {
            main.taMsgs.append("Error: "+e.getMessage());
            main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
            e.printStackTrace();
        }
    }
    
    /**
     * Destroy appropriate processes
     */
    private void killProcesses() {
        if(cppProcess!=null) {
            main.taMsgs.append("Killing main Phycas process\n");
            main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
            System.out.println("Killing main Phycas process\n");
            cppProcess.destroy();
        }
        if(guiProcess!=null) {
            main.taMsgs.append("Killing Phycas user interface\n");
            main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
            System.out.println("Killing Phycas user interface\n");
            guiProcess.destroy();
        }
    }
    
    /**
     * Set appropriate process to null
     * @param i process id
     */
    public void resetProcess(int i) {
        if(i==0)
            cppProcess = null;
        else
            guiProcess = null;
    }
    
    public Action okAction = new AbstractAction() {
        public void actionPerformed(ActionEvent e) {
            handleStartButton();
        }
    };
    
    /**
     * Called when start button action event received
     */
    private void handleStartButton() {
        if(!validateAllFields())
            return;
        begin();
    }
    
    /**
     * Close program gracefully.
     * @param e window event
     */
    public void windowClosing(WindowEvent e) {
        end();
    }
    
    // All of these are required for implementing WindowListener
    // leave empty so generic handling will occur
    
    public void windowActivated(WindowEvent e) {
    }
    
    public void windowClosed(WindowEvent e) {
    }
    
    public void windowDeactivated(WindowEvent e) {
    }
    
    public void windowDeiconified(WindowEvent e) {
    }
    
    public void windowIconified(WindowEvent e) {
    }
    
    public void windowOpened(WindowEvent e) {
    }
    
}
