/* 
 VKJ 12/04 This file is no longer being used, as it was a component
 of PhycasGraph.  Use PhycasChart instead.
 */ 

package phycasGUI.swixml;

import java.io.File;
import javax.swing.*;
import javax.swing.filechooser.*;

/**
 * File Filter to allow only jpeg,jpg,or png files for saving images
 * @author Vanessa Jackson
 */

public class PhycasImageFilter  extends FileFilter {
    
    public final static String JPEG = "jpeg";
    public final static String JPG = "jpg";
    public final static String PNG = "png";
    
    /** Creates a new instance of PhycasImageFilter */
    public PhycasImageFilter() {
    }
    
    /**
     * Get the file extension
     * @param f instance of file
     * @return extension
     */    
    public static String getExtension(File f) {
        String ext = null;
        String s = f.getName();
        int i = s.lastIndexOf('.');

        if (i > 0 &&  i < s.length() - 1) {
            ext = s.substring(i+1).toLowerCase();
        }
        return ext;
    }
    
    /**
     * Restrict extensions of files that are acceptable
     * @param f instance of file
     * @return <CODE>true</CODE> if acceptable extension
     */    
    public boolean accept(File f) {
        if (f.isDirectory()) {
            return true;
        }

        String extension = getExtension(f);
        if (extension != null) {
            if (extension.equals(JPEG) ||
                extension.equals(JPG) ||
                extension.equals(PNG)) {
                    return true;
            } else {
                return false;
            }
        }

        return false;
    }

    
    /**
     * The description of this filter
     * @return description
     */    
    public String getDescription() {
        return "Image File Type JPEG or PNG";
    }
}