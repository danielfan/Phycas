/*
 * PhycasOutput.java
 *
 * Created on December 15, 2004, 3:15 PM
 */

package phycasGUI.swixml;

import java.awt.*;
import java.awt.event.*;
import java.awt.print.*;
import javax.swing.*;
import java.io.*;
/**
 * Handle the Output window
 * @author vjackson
 */
public class PhycasOutput extends AbstractAction implements ActionListener, MouseListener {
    
    private PhycasMain main = null;
    private JPopupMenu popup = null;
    
    /**
     * Creates a new instance of PhycasOutput
     * @param m Phycas Main object
     */
    public PhycasOutput(PhycasMain m) {
        main = m;
        popup = createPopupMenu("Output");
    }
    
    /**
     * Handle an event
     * @param e event
     */    
    public void actionPerformed(ActionEvent e) {
        if(e.getActionCommand().equals("showAction_Output")) {
            main.output.setVisible(true);
        }
        
        if(e.getActionCommand().equals("closeAction_Output")) {
            main.output.setVisible(false);
        }
        
        if(e.getActionCommand().equals("saveAction_Output")) {
            doSaveAs();
        }
        
        if(e.getActionCommand().equals("clearAction_Output")) {
            int confirm = JOptionPane.showConfirmDialog(
            main.output,
            "Note: This action will permanently erase the output text.\nPlease confirm that you have " +
            "saved and are ready to clear the output text.",
            "Confirm Clear",
            JOptionPane.OK_CANCEL_OPTION);
            System.out.println("RESPONSE = "+confirm);
            if(confirm==0)
                main.taOutput.setText(null);
        }
        
    }
    
    private void doSaveAs() {
        PrintWriter pw = null;
        File file = null;
        try {
            String out = main.taOutput.getText();
            int confirm = 1;
            int retVal = -1;
            while(confirm==1) {
                retVal = main.fc.showSaveDialog(main.output);
                if (retVal == JFileChooser.APPROVE_OPTION) {
                    file = main.fc.getSelectedFile();
                    if(file.exists()) {
                        confirm = JOptionPane.showConfirmDialog(
                        main.output,
                        "This file already exists.  Do you want to replace it?",
                        "File Exists",
                        JOptionPane.YES_NO_OPTION);
                        System.out.println("RESPONSE = "+confirm);
                    }
                    else
                        confirm = 0;
                }
                else
                    confirm = 0;
            }
            if(retVal==JFileChooser.APPROVE_OPTION) {
                pw = new PrintWriter(new FileWriter(file));
                pw.print(out);
                pw.flush();
                main.taMsgs.append("\nSaved output to file "+file.getAbsolutePath());
                main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
            main.taMsgs.append("\nFailed to write to file because ");
            main.taMsgs.append(ex.getMessage()+"\n");
            main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
        }
        finally {
            if(pw!=null)
                pw.close();
            main.fc.setSelectedFile(new File(""));
        }
        
    }
    
    /**
     * Creates a popup menu for a panel.
     *
     * @return The popup menu.
     */
    private JPopupMenu createPopupMenu(String name) {
        
        JPopupMenu result = new JPopupMenu();
        
        Font f = new Font("Tahoma",Font.PLAIN,11);
        
        JMenuItem saveItem = new JMenuItem("Save As");
        saveItem.setMnemonic('S');
        saveItem.setFont(f);
        saveItem.setBackground(Color.WHITE);
        saveItem.setActionCommand("saveAction_"+name);
        saveItem.addActionListener(this);
        result.add(saveItem);
        result.addSeparator();
        
        JMenuItem clearItem = new JMenuItem("Clear");
        clearItem.setMnemonic('L');
        clearItem.setFont(f);
        clearItem.setBackground(Color.WHITE);
        clearItem.setActionCommand("clearAction_"+name);
        clearItem.addActionListener(this);
        result.add(clearItem);
        
        return result;
        
    }
    
    /**
     * Handle mouse event
     * @param e event
     */    
    public void mousePressed(MouseEvent e) {
        maybeShowPopup(e);
    }
    
    /**
     * Handle mouse event
     * @param e event
     */    
    public void mouseReleased(MouseEvent e) {
        maybeShowPopup(e);
    }
    
    private void maybeShowPopup(MouseEvent e) {
        if (e.isPopupTrigger())
            popup.show(e.getComponent(),e.getX(), e.getY());
    }
    
    // These are not used, but must be here
    public void mouseClicked(MouseEvent e) {
    }
    
    public void mouseEntered(MouseEvent e) {
    }
    
    public void mouseExited(MouseEvent e) {
    }
    
    
}
