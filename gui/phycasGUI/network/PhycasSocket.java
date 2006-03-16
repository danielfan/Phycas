
package phycasGUI.network;

import phycasGUI.swixml.*;
import java.io.*;
import java.net.*;
import org.jdom.input.SAXBuilder;
import org.jdom.*;
import org.jdom.output.*;
import java.util.*;
import javax.swing.JOptionPane;

/**
 * Handles the socket connection of the gui to the main program
 * @author Vanessa Jackson
 */
public class PhycasSocket extends Thread {
    private Socket socket = null;
    private PrintWriter out = null;
    private BufferedReader in = null;
    private InputStream inputStream = null;
    private PhycasMain main = null;
    
    private String xml_hiddenQuery = null;
    private String xml_userQuery = null;
    private String xml_plot = null;
    private String xml_out = null;
    private String xml_error = null;
    private String xml_comment = null;
    private String xml_warning = null;
    
    private  static final String HQ = "hidden_query";
    private  static final String OUT = "out";
    private  static final String ERROR = "error";
    private  static final String COMMENT = "comment";
    private  static final String WARNING = "warning";
    private  static final String IDLE = "idle";
    private  static final String UQ = "user_query";
    private  static final String PLOT = "plot";
    private String host = null;
    private int port = -1;
    
    /**
     * Constructor
     * @param h host name
     * @param p port number
     */
    public PhycasSocket(String h, int p) {
        host = h;
        port = p;
    }
    
    /**
     * Set the instance of PhycasMain
     * @param m instance of main gui window
     */
    public void setMain(PhycasMain m) {
        main = m;
    }
    
    /**
     * Open a socket connection to the main program
     * @throws IOException unable to open socket
     * @throws UnknownHostException host cannot be found
     */
    public void openSocket() throws IOException,UnknownHostException {
        
        socket = new Socket(host, port);
        
        //POL 25-Sep-2004 don't want autoflushing
        //out = new PrintWriter(socket.getOutputStream(), true);
        out = new PrintWriter(socket.getOutputStream(), false);
        
        inputStream = socket.getInputStream();
        in = new BufferedReader(new InputStreamReader(inputStream));
    }
    
    /**
     * Write a string to the socket
     * @param s string to write
     * @throws Exception any exception
     */
    public synchronized void writeToSocket(String s) throws Exception {
        if(socket==null)
            openSocket();
        out.println(s);
        //VKJ 29-Sep-2004 Windows needs new line character when
        //  receiving message through socket
        //POL 25-Sep-2004 don't want line feed characters
        //out.print(s);
        out.flush();
        System.out.println("wrote - "+s);
    }
    
    /**
     * Close any open streams and the socket
     */
    public void closeSocket() {
        try {
            if(out!=null)
                out.close();
            if(in!=null)
                in.close();
            if(socket!=null && !socket.isClosed())
                socket.close();
        }
        catch (IOException io) {
            io.printStackTrace();
        }
    }
    
    /**
     * Handle xml from socket according to tag
     * @param current main tag
     * @return 0 - success, 1 - failure
     */
    public int assignChildren(Element current) throws Exception {
        //System.out.println("AC element name = "+current.getName());
        //System.out.println("AC element text = "+current.getText());
        List children = current.getChildren();
        Iterator iterator = children.iterator();
        while (iterator.hasNext()) {
            Element child = (Element) iterator.next();
            System.out.println("Tag = "+child.getName());
            if(child.getName().equals(HQ))
                xml_hiddenQuery = child.getText();
            else if(child.getName().equals(UQ))
                xml_userQuery = child.getText();
            else if(child.getName().equals(PLOT))
                xml_plot = child.getText();
            else if(child.getName().equals(OUT))
                xml_out = child.getText();
            else if(child.getName().equals(ERROR))
                xml_error = child.getText();
            else if(child.getName().equals(COMMENT))
                xml_out = child.getText();
            else if(child.getName().equals(WARNING))
                xml_warning = child.getText();
            else if(child.getName().equals(IDLE))
                System.out.println("IDLE");
            else {
                System.out.println("INVALID XML TAG");
                return 1;
            }
            if(xml_out!=null) {
                if(!xml_out.trim().equals("")) {
                    main.taOutput.append("\n"+xml_out);
                    main.taOutput.setCaretPosition(main.taOutput.getDocument().getLength());
                    main.output.requestFocus();
                }
                xml_out = null;
            }
            if(xml_error!=null){
                main.goToLastPageViewed();
                JOptionPane.showMessageDialog(main.pnl,xml_error, "ERROR",
                JOptionPane.ERROR_MESSAGE);
                xml_error = null;
            }
            if(xml_warning!=null){
                if(main.mi_settings_displayWarnings.isSelected()) {
                    main.goToLastPageViewed();
                    JOptionPane.showMessageDialog(main.pnl,xml_warning, "WARNING",
                    JOptionPane.WARNING_MESSAGE);
                }
                else if(main.mi_settings_logWarnings.isSelected() && !xml_warning.trim().equals("")) {
                    main.taOutput.append("\n WARNING > "+xml_warning);
                    main.taOutput.setCaretPosition(main.taOutput.getDocument().getLength());
                    main.output.requestFocus();
                }
                xml_warning = null;
            }
            if(xml_plot!=null){
                main.phycasChart.assignChildren(child);
            }
            if(xml_userQuery!=null){
                main.goToLastPageViewed();
                PhycasUserQuery uq = new PhycasUserQuery(main.frame);
                String response = uq.showDialog(child);
                try {
                    if(response!=null) {
                        main.socket.writeToSocket(response);
                    }
                    else {
                        main.taMsgs.append("\nUnable to write response.");
                        main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                    }
                }
                catch (Exception e) {
                    main.taMsgs.append("\nUnable to write response.");
                    main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                    e.printStackTrace();
                }
                xml_userQuery = null;
            }
            /** ~=================================================================== */
            if(xml_hiddenQuery!=null) { //@Use when command_state ready
                System.out.println("Found hidden_query");
                // cut off hidden_query tag before parsing
                Element cl = child.getChild("command_state");
                cl.setNamespace(Namespace.getNamespace("http://commandState.phycas.org"));
                System.out.println("Namespace in Element = "+cl.getNamespace());
                Document doc = new Document((Element)cl.clone());
                System.out.println("Command_State XML:");
                try {
                    XMLOutputter outputter = new XMLOutputter();
                    outputter.output(doc, System.out);
                }
                catch (IOException e) {
                    e.printStackTrace();
                }
                DOMOutputter domOut = new DOMOutputter();
                org.w3c.dom.Document domDoc = null;
                try {
                    domDoc = domOut.output(doc);
                }
                catch (JDOMException e) {
                    e.printStackTrace();
                }
                PhycasCommandState cs = new PhycasCommandState(main);
                cs.parse(domDoc);
                xml_hiddenQuery = null;
            }
            /* //@~======================================================================*/
        }
        return 0;
    }
    
    /**
     * Read from socket
     */
    public void run() {
        try{
            openSocket();
            
            while(true){
                String ret = "";
                char[] cbuf = new char[1000];
                int countNulls = 0;
                int numBytesRead = -1;
                StringBuffer data = new StringBuffer();
                int count = 1;
                boolean readAgain = true;
                // need base tag in order to properly parse xml message
                data.append("<doc>");
                while(readAgain) {
                    //System.out.println("Reading group "+count);
                    cbuf = new char[1000];
                    numBytesRead = in.read(cbuf);
                    //System.out.println("read "+numBytesRead+" chars with "+countNulls+" nulls");
                    if(numBytesRead < 0)
                        break; //@ this means that the socket connection has died, we should reopen or throw and exception.
                    // move to new buffer in case read less than 1000 chars
                    char[] dbuf = new char[numBytesRead];
                    System.arraycopy(cbuf,0,dbuf,0,numBytesRead);
                    String s = new String(dbuf);
                    // read until either idle tag or end of out or user_query tag
                    boolean foundEnd = false;
                    if(s.endsWith("</out>") || s.endsWith("<idle/>") || s.endsWith("</user_query>")  || s.endsWith("</plot>"))
                        foundEnd = true;
                    else if(s.endsWith("<idle/>\n") || s.endsWith("<idle/>\r") || s.endsWith("</out>\n") || s.endsWith("</out>\r") || s.endsWith("<out/>\n") || s.endsWith("<out/>\r") || s.endsWith("</plot>\n") || s.endsWith("</plot>\r") )
                        foundEnd = true;
                    if(numBytesRead<1000 && foundEnd)
                        readAgain = false;
                    data.append(s);
                    //System.out.println("Appended: "+s);
                    //System.out.println("GOT: "+data);
                    count++;
                }
                data.append("</doc>");
                char[] newBuf = new char[data.length()];
                data.getChars(0,data.length(),newBuf,0);
                System.out.println("Buffer START\n"+new String(newBuf)+"\nBuffer END");
                //System.out.println("NewBuf size = "+newBuf.length);
                
                // parse XML message
                Reader r = new CharArrayReader(newBuf);
                SAXBuilder parser = new SAXBuilder();
                //System.out.println("got parser");
                Document doc = parser.build(r);
                //System.out.println("built doc");
                Element root = doc.getRootElement();
                //System.out.println("got root");
                //System.out.println("current element name = "+root.getName());
                //System.out.println("current element text = "+root.getText());
                int retVal = assignChildren(root);
                if(retVal!=0)
                    throw new Exception("Error reading children");
            }
        }
        catch (UnknownHostException e) {
            System.err.println("Don't know about host: "+host);
            showErrMsgAndExit();
        }
        catch (JDOMException e) {
            System.err.println("XML is not well-formed.");
            System.err.println(e.getMessage());
            showErrMsgAndExit();
        }
        catch (IOException e) {
            System.err.println("IO Error");
            System.err.println(e.getMessage());
            showErrMsgAndExit();
        }
        catch (Exception e) {
            e.printStackTrace();
            showErrMsgAndExit();
        }
    }
    
    private void showErrMsgAndExit() {
        JOptionPane.showMessageDialog(main.pnl,"An error has occurred in the socket, so the program Phycas must exit.", "ERROR",
        JOptionPane.ERROR_MESSAGE);
        closeSocket();
        System.exit(0);
    }
    
    /**
     * Get the xml in the hidden_query tag
     * Need to deprecate once command_state handling is working
     * @return xml inside tag
     */
    public synchronized String getHiddenQuery() {
        System.out.println("Waiting for hidden query ...");
        while (true) {
            if(xml_hiddenQuery!=null)
                break;
        }
        System.out.println("Returning "+xml_hiddenQuery);
        String ret = xml_hiddenQuery.toString();
        xml_hiddenQuery = null;
        return ret;
    }
    
}
