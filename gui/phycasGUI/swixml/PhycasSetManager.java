package phycasGUI.swixml;

import java.util.*;

/**
 * Manage set data
 * @author Vanessa Jackson
 */
public class PhycasSetManager {
    
    private PhycasSet[] allSets = null;
    private int maxIndex = -1;
    private String type = "";
    private Hashtable indexLabel = null;
    
    /**
     * Constructor
     * @param t type - either Char, Tree, or Tax
     */    
    public PhycasSetManager(String t) {
        type = t;
    }

    /**
     * List of possible components of set
     * @return list of data
     */    
    public String[] getListData() {
        String[] data = null;
        if(maxIndex!=-1)
            data = new String[maxIndex];
        else
            data = new String[indexLabel.size()];
        for(int i=0;i<data.length;i++) {
            int j = i+1;
            if(indexLabel!=null && indexLabel.get(new Integer(j))!=null) 
                data[i] = (String)indexLabel.get(new Integer(j));
            else
                data[i] = new Integer(j).toString();
        }
        return data;
    }
    
    /**
     * Add a set to this manager
     * @param s set to add
     */    
    public void addSet(PhycasSet s) {
        if(allSets == null) {
            allSets = new PhycasSet[1];
            allSets[0] = s;
        }
        else {
            allSets = mergeArrays(allSets,new PhycasSet[]{s});
        }
    }
    
    /**
     * Accessor for set in this manager
     * @param i index of set
     * @return set
     */    
    public PhycasSet getSet(int i) {
        if(allSets.length<i)
            return allSets[i];
        return null;
    }
    
    /**
     * Accessor for set in this manager
     * @param label label (name) of set
     * @return set
     */    
    public PhycasSet getSet(String label) {
        for(int i=0;i<allSets.length;i++) {
            if(label.equals(allSets[i].getLabel()))
                return allSets[i];
        }
        return null;
    }
    
    /**
     * Accessor for all sets in manager
     * @return all sets
     */    
    public PhycasSet[] getAllSets() {
        return allSets;
    }
    
    /**
     * Accessor for maximum index of possible components of set
     * @return maximum index
     */    
    public int getMaxIndex() {
        return maxIndex;
    }
    
    /**
     * Mutator for maximum index of possible components of set
     * @param i maximum index
     */    
    public void setMaxIndex(int i) {
        maxIndex = i;
    }
    
    /**
     * Associate a label to an index of a possible component of set
     * @param index index
     * @param label label
     */    
    public void addIndexLabel(int index, String label) {
        if(indexLabel == null)
            indexLabel = new Hashtable();
        indexLabel.put(new Integer(index),label);
    }

    /**
     * Get the index of a certain label
     * @param label label to get index for
     * @return index
     */    
    public int getIndexForLabel(String label) {
        Enumeration en = indexLabel.keys();
        while(en.hasMoreElements()) {
            Integer key = (Integer)en.nextElement();
            String element = (String)indexLabel.get(key);
            if(element.equals(label)) {
                return key.intValue();
            }
        }
        return -1;
    }
    
    /**
     * Get the Hashtable of all index/label pairs
     * @return all index/label pairs
     */    
    public Hashtable getIndexLabel() {
        return indexLabel;
    }

    private PhycasSet[] mergeArrays(PhycasSet[] pa, PhycasSet[] pb) {
        PhycasSet[] arr = new PhycasSet[pa.length+pb.length];
        for (int x=0; x < pa.length; x++) {
            arr[x] = pa[x];
        }
        for (int x=0; x < pb.length; x++) {
            arr[x+pa.length] = pb[x];
        }
        return arr;
    }
    
}

