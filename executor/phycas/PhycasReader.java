
package phycas;

import java.io.*;

/**
 * Connect output and error streams of processes to output of Executor
 * @author Vanessa Jackson
 */
public class PhycasReader extends BufferedReader implements Runnable {
    
    private Thread thread = null;
    private String type = "";
    private String processType = "";
    private PhycasExecutor main = null;
    
    /**
     * Creates a new instance of PhycasReader
     * @param is InputStream of process
     * @param i id of process - 0 = Main, 1 = GUI
     * @param t type of process stream - 0 = output, 1 = error
     */
    public PhycasReader(PhycasExecutor m, InputStream is, int i, int t) {
        super(new InputStreamReader(is));
        if(i==0)
            processType = "Main";
        else
            processType = "GUI";
        if(t==0)
            type="Out";
        else
            type="Err";
        main = m;
    }
    
    /**
     * Start thread
     */    
    public void start() {
        if (thread == null) {
            thread = new Thread(this);
            thread.start();
        }
    }
    
    /**
     * Forward output and errors from process to output of Executor
     */    
    public void run() {
        try {
            System.out.println("Reader thread for process "+processType+" starting");
            String line = readLine();
            while (line!=null) {
                System.out.println(type+"_"+processType+" > "+line);
                if(type.equals("Err") || line.toUpperCase().indexOf("EXCEPTION")!=-1 || line.toUpperCase().indexOf("ERROR")!=-1) {
                    main.taMsgs.append(processType+" > "+line+"\n");
                    main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                }
                line = readLine();
            }
            System.out.println("Reader thread for process "+processType+" ending");
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    
}
