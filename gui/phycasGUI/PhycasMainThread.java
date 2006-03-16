package phycasGUI;

import phycasGUI.network.*;
import phycasGUI.swixml.*;

/**
 * Starts the socket and gui program threads for the PhycasGUI
 * @author Vanessa Jackson
 */
public class PhycasMainThread {
    
    /** Creates a new instance of PhycasMainThread */
    public PhycasMainThread() {
    }
    
    /**
     * Spawn threads for the socket and gui
     * Expects length of argc to be:
     *  -   0 -> (host = localhost, port = 4444)
     *  -   1 -> (host = localhost port =argc[0])
     *  -   2 -> (host = argc[0] port = argc[1])
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            String host = "localhost";// mth added default arguments that match phycas_server defaults
            int port = 4444; 
            if (args.length > 2) {
                System.out.println("The only accepted arguments are host-IP and port.");
                return;
            }
            else if (args.length == 2){
                host = args[0];
                port = new Integer(args[1]).intValue();
            }
            else if (args.length == 1){
                port = new Integer(args[0]).intValue();
            }
            PhycasSocket s = new PhycasSocket(host,port);
            PhycasMain m = new PhycasMain(s);
            s.start();
            m.start();
            s.setMain(m);
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
