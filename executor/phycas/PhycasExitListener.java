
package phycas;

/**
 * A separate thread to wait for process to end
 * @author Vanessa Jackson
 */
public class PhycasExitListener extends Thread {
    
        int exitVal = -1;
        int processNum = -1;
        Process process = null;
        PhycasExecutor main = null;
        
        /**
         * Constructor
         * @param m instance of PhycasExecutor containing window and it's components
         * @param p Process to monitor
         * @param i id of Process to monitor
         */        
        public PhycasExitListener(PhycasExecutor m, Process p, int i) {
            main = m;
            process = p;
            processNum = i;
        }
        
        /**
         * Run the thread, waiting for process to end
         */        
        public void run() {
            try {
                if(processNum==0) 
                    System.out.println("Phycas main process starting");
                else 
                    System.out.println("Phycas user interface starting");
                int exitVal = process.waitFor();
                if(exitVal != 0) {
                    if(processNum==0) {
                        System.out.println("Phycas main process has exited with error "+exitVal);
                        main.taMsgs.append("Phycas main process has exited with error\n");
                        main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                    }
                    else {
                        System.out.println("Phycas user interface has exited with error "+exitVal);
                        main.taMsgs.append("Phycas user interface has exited with error\n");
                        main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                    }
                }
                else {
                    if(processNum==0) {
                        System.out.println("Phycas main process has exited successfully ");
                        main.taMsgs.append("Phycas main process and user interface have exited.\n");
                        main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                    }
                    else {
                        System.out.println("Phycas user interface has exited successfully");
                        main.taMsgs.append("Phycas user interface has exited, but the main process may still running.\n");
                        main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                    }
                }
                main.handler.resetProcess(processNum);
            }
            catch(Exception e) {
                e.printStackTrace();
            }
        }
    
}
