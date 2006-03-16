package phycasGUI.swixml;

import java.util.StringTokenizer;

/**
 * Contains set data
 * @author Vanessa Jackson
 */

public class PhycasSet {
    
    private String members = "";
    private String label = "";
    
    /**
     * Constructor
     * @param l label (name) of set
     */
    public PhycasSet(String l) {
        label = l;
    }
    
    /**
     * Mutator for members of set
     * @param m members of set
     */
    public void setMembers(String m) {
        members = m;
    }
    
    /**
     * Accessor for members of set
     * @return members of set
     */
    public String[] getMembers() {
        String[] s = null;
        StringTokenizer st = new StringTokenizer(members, "-\\ ", true);
        int countTokens = 0;
        while(st.hasMoreTokens()) {
            String t = st.nextToken();
            if(!t.equals(" ")) {
                countTokens++;
            }
        }
        st = new StringTokenizer(members, "-\\ ", true);
        String[] tokens = new String[countTokens];
        int k=0;
        int count = 0;
        while(st.hasMoreTokens()) {
            String t = st.nextToken();
            if(!t.equals(" ")) {
                tokens[k] = t;
                if(!tokens[k].equals("-") && !tokens[k].equals("\\"))
                    count++;
                k++;
            }
        }
        // check if space separated list of numbers or nexus format (1-3,or 1-10\2)
        if(members.indexOf('-')!=-1) {
            s = new String[0];
            String begin = tokens[0];
            // see if followed by -
            String test = "";
            if(tokens.length>1)
                test = tokens[1];
            if(!test.equals("-")) {
                String[] temp = new String[1];
                temp[0] = begin;
                s = mergeArrays(s,temp);
            }
            for(int i=1;i<tokens.length;i++) {
                String next = tokens[i];
                String end = "", stride = "", extra = "";
                if(next.equals("-")) {
                    end = tokens[i+1];
                    int f = new Integer(begin).intValue();
                    int l = new Integer(end).intValue();
                    if(i+2>=tokens.length) {
                        String[] temp = new String[(l-f)+1];
                        int c = 0;
                        for(int r=f;r<=l;r++) {
                            temp[c] = new Integer(r).toString();
                            c++;
                        }
                        s = mergeArrays(s,temp);
                        i++;
                    }
                    else if(tokens[i+2].equals("\\")) {
                        stride = tokens[i+3];
                        int j = 0;
                        int strideInt = new Integer(stride).intValue();
                        String[] temp = new String[((l-f)+strideInt)/strideInt];
                        for(int r=f;r<=l;r+=strideInt) {
                            temp[j] = new Integer(r).toString();
                            j++;
                        }
                        s = mergeArrays(s,temp);
                        i+=3;
                    }
                    else {
                        String[] temp = new String[(l-f)+1];
                        int c = 0;
                        for(int r=f;r<=l;r++) {
                            temp[c] = new Integer(r).toString();
                            c++;
                        }
                        s = mergeArrays(s,temp);
                        i++;
                    }
                }
                else {
                    // see if followed by -
                    test = "";
                    if(i+1<tokens.length)
                        test = tokens[i+1];
                    if(test.equals("-")) {
                        begin = next;
                    }
                    else {
                        String[] temp = new String[1];
                        temp[0] = next;
                        s = mergeArrays(s,temp);
                    }
                }
            }
        }
        else {
            st = new StringTokenizer(members, " ");
            s = new String[st.countTokens()];
            int i = 0;
            while(st.hasMoreTokens()) {
                s[i] = st.nextToken();
                i++;
            }
        }
        return s;
    }
    
    /**
     * Accessor for label of set
     * @return label (name) of set
     */
    public String getLabel() {
        return label;
    }
    
    // merge two arrays into one
    private String[] mergeArrays(String[] pa, String[] pb) {
        String[] arr = new String[pa.length+pb.length];
        for (int x=0; x < pa.length; x++) {
            arr[x] = pa[x];
        }
        for (int x=0; x < pb.length; x++) {
            arr[x+pa.length] = pb[x];
        }
        return arr;
    }
    
}

