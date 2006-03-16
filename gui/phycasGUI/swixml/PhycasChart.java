
package phycasGUI.swixml;

import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;

import org.jdom.input.SAXBuilder;
import org.jdom.*;

import org.jfree.chart.axis.*;
import org.jfree.chart.renderer.xy.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.ui.*;
import org.jfree.chart.*;
import org.jfree.data.xy.*;
import org.jfree.ui.*;


/**
 * Trace Plot Graphing
 * @author  Vanessa Jackson
 */
public class PhycasChart implements ActionListener, ItemListener {
    
    private PhycasMain main = null;
    private String xAxisLabel = "";
    private String title = "";
    private XYSeries[] dataSeries = null;
    private Hashtable data = null; //key = label, value = double[]
    private JPanel choices = null;
    private String[] labels = null;
    private ArrayList selectedLabels = null;
    ButtonGroup xGroup = null;
    
    /**
     * Creates a new instance of PhycasChart
     * @param m instance of PhycasMain
     */
    public PhycasChart(PhycasMain m) {
        main = m;
    }
    
    /**
     * Add the main graph tab to the output window.
     * Contains options for graphs.
     * @throws Exception any exception
     */
    private void addGraphMain() throws Exception {
        if(labels==null || labels.length==0)
            return;
        
        JPanel main = new JPanel();
        main.setLayout(new BoxLayout(main,BoxLayout.X_AXIS));
        
        JPanel x = new JPanel();
        x.setLayout(new BoxLayout(x,BoxLayout.Y_AXIS));
        JLabel xTitle = new JLabel("X AXIS");
        xTitle.setFont(new Font("Arial", Font.BOLD, 12));
        xTitle.setAlignmentX(JLabel.LEFT_ALIGNMENT);
        x.add(xTitle);
        xGroup = new ButtonGroup();
        // add radio buttons for x axis since can only choose one
        for(int i=0;i<labels.length;i++) {
            JRadioButton rb = new JRadioButton(labels[i]);
            rb.setAlignmentX(JRadioButton.LEFT_ALIGNMENT);
            rb.setActionCommand(labels[i]);
            xGroup.add(rb);
            x.add(rb);
            if(i==0)
                rb.setSelected(true);
        }
        
        JPanel y = new JPanel();
        y.setLayout(new BoxLayout(y,BoxLayout.Y_AXIS));
        JLabel yTitle = new JLabel("Y AXIS");
        yTitle.setFont(new Font("Arial", Font.BOLD, 12));
        yTitle.setAlignmentX(JLabel.LEFT_ALIGNMENT);
        y.add(yTitle);
        selectedLabels = new ArrayList();
        // add check boxes for y axis since can choose one or more
        for(int i=0;i<labels.length;i++) {
            JCheckBox cb = new JCheckBox(labels[i]);
            cb.setAlignmentX(JCheckBox.LEFT_ALIGNMENT);
            cb.addItemListener(this);
            y.add(cb);
        }
        
        main.add(x);
        main.add(y);
        
        JLabel choose = new JLabel("Graphing Options");
        choose.setFont(new Font("Arial", Font.BOLD, 12));
        choose.setAlignmentX(JLabel.CENTER_ALIGNMENT);
        
        JButton submit = new JButton("Submit");
        submit.setActionCommand("Submit");
        submit.addActionListener(this);
        
        choices = new JPanel();
        choices.setLayout(new BoxLayout(choices,BoxLayout.PAGE_AXIS));
        choices.add(choose);
        choices.add(Box.createRigidArea(new Dimension(20, 20)));
        choices.add(main);
        choices.add(Box.createRigidArea(new Dimension(20, 20)));
        choices.add(submit);
    }

    /**
     * NOT YET IMPLEMENTED!
     * Update the main graph tab on the output window.
     * Contains options for other graphs.
     * @throws Exception any exception
     */
    private void updateGraphMain() throws Exception {
        // Need to implement!
        // if there are new labels, add radioButtons & checkBoxes to choices panel
    }
    
    /**
     * Create a graph
     * @param name Name of the series of data and the graph, corresponds to y-axis label.
     * @return JFreeChart ChartPanel
     */
    private ChartPanel createGraph(String name) {
        NumberAxis xAxis = new NumberAxis(xAxisLabel);
        NumberAxis yAxis = new NumberAxis(); 
        
        yAxis.setAutoRangeIncludesZero(false);
        xAxis.setAutoRangeIncludesZero(false);
        
        XYItemRenderer renderer = new StandardXYItemRenderer(StandardXYItemRenderer.LINES);
        int current = -1;
        XYSeriesCollection xyc = new XYSeriesCollection();
        dataSeries = new XYSeries[labels.length];
        // get data for currently selected x and y values and add to XYDataSeries
        for(int i=0;i<selectedLabels.size();i++) {
            ArrayList x = (ArrayList)data.get(xAxisLabel);
            String s = (String)selectedLabels.get(i);
            ArrayList y = (ArrayList)data.get(s);
            addDataSeries(xAxisLabel, x.toArray(), s,y.toArray());
        }
        // add all XYDataSeries to XYSeriesCollection
        for(int i=0;i<selectedLabels.size();i++)
            xyc.addSeries(dataSeries[i]);
        
        XYPlot myplot=new XYPlot(xyc, xAxis, yAxis, renderer);
        myplot.setOrientation(PlotOrientation.VERTICAL);
        
        JFreeChart graph = new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT, myplot, true);
        graph.setBorderVisible(true);
        graph.setBorderPaint(Color.BLACK);
        graph.setBackgroundPaint(Color.white);
        
        ChartPanel panel = new ChartPanel(graph,true,true,true,false,false);
        panel.setAutoscrolls(true);
        panel.setVisible(true);
        JPopupMenu popup = panel.getPopupMenu();
        JMenuItem closeItem = new JMenuItem("Close");
        closeItem.setActionCommand("closeAction_Graph_"+name);
        closeItem.addActionListener(this);
        popup.addSeparator();
        popup.add(closeItem);
        
        return panel;
    }
    
    /**
     * Reads an xml file containing data to graph
     * @param filename name of file containing data
     * @throws Exception any exception
     */
    public void readFile(String filename) throws Exception {
        File file = new File(filename);
        SAXBuilder parser = new SAXBuilder();
        Document doc = parser.build(file);
        Element root = doc.getRootElement();
        int retVal = assignChildren(root);
        if(retVal!=0)
            throw new Exception("Error reading children");
    }
    
    /**
     * Add data to a DataSeries or create one
     * @param xName Name of the x axis
     * @param xData X Axis Data to add to graph
     * @param yName Name of the series, corresponding to name of y-axis label
     * @param yData Y Axis Data to add to graph
     */
    private void addDataSeries(String xName, Object[] xData, String yName, Object[] yData){
        int current = -1;
        int last = -1;
        for(int i=0;i<dataSeries.length;i++){
            if(dataSeries[i]==null) {
                last = i;
                break;
            }
            else if(dataSeries[i].getName().equals(yName)) {
                current = i;
            }
        }
        XYSeries series = null;
        if(current==-1)
            series = new XYSeries(yName);
        else
            series = (XYSeries)dataSeries[current];
        for(int i=0;i<yData.length;i++){
            double d0 = new Double((String)xData[i]).doubleValue();
            double d1 = new Double((String)yData[i]).doubleValue();
            series.add(d0,d1);
        }
        if(current==-1) {
            dataSeries[last] = series;
        }
    }
    
    /**
     * Parse MCMC graphing data from xml containing labels and/or entries
     * @param current base tag
     * @throws Exception any exception
     * @return 0 - success, 1 - failure
     */
    public int assignChildren(Element current) throws Exception {
        /* format of xml for this graph is as follows:
         a tab-separated list of labels of the values to plot
         <plot><label>Gen    LnL	TL	pi(A)	pi(C)	pi(G)	pi(T)</label></plot>
         a tab-separated list of values to plot, in order corresponding to labels
         <plot><entry>1	-23560.466	10.500	0.250000	0.250000	0.250000	0.250000	</entry></plot>
         Each of these lines will come in separately as updates to the gui from the cpp program
         */
        
        java.util.List children = current.getChildren();
        Iterator iterator = children.iterator();
        String xml_label = null;
        String xml_entry = null;
        StringTokenizer st = null;
        while (iterator.hasNext()) {
            Element child = (Element) iterator.next();
            if(child.getName().equals("label"))
                xml_label = child.getText();
            else if(child.getName().equals("entry"))
                xml_entry = child.getText();
            else {
                System.out.println("INVALID XML TAG");
                return 1;
            }
            if(xml_label!=null) {
                st = new StringTokenizer(xml_label);
                int num = st.countTokens();
                labels = new String[num];
                for(int i=0;i<num;i++) {
                    labels[i] = st.nextToken();
                }
                boolean alreadyExists = false;
                for(int i=0;i<main.tpOutput.getTabCount();i++) {
                    String title = main.tpOutput.getTitleAt(i);
                    if(title.equals("Graph - Main"))
                        alreadyExists = true;
                }
                if(!alreadyExists && choices==null)
                    addGraphMain();
                else
                    updateGraphMain();
                // clear out old data from previous plot
                // assumes one plot finishes before another begins
                dataSeries = new XYSeries[labels.length];
                data = new Hashtable(); 
            }
            else if(xml_entry!=null) {
                boolean alreadyExists = false;
                for(int i=0;i<main.tpOutput.getTabCount();i++) {
                    String title = main.tpOutput.getTitleAt(i);
                    if(title.equals("Graph - Main"))
                        alreadyExists = true;
                }
                if(!alreadyExists)
                    main.tpOutput.addTab("Graph - Main", choices);
                if(labels!=null) {
                    st = new StringTokenizer(xml_entry);
                    int num = st.countTokens();
                    String[] tokens = new String[num];
                    for(int i=0;i<num;i++) { 
                        tokens[i] = st.nextToken();
                    }
                    for(int j = 0;j<labels.length;j++) {
                        ArrayList d = null;
                        if(data!=null)
                            d = (ArrayList)data.get(labels[j]);
                        else
                            data = new Hashtable();
                        if(d==null)
                            d = new ArrayList();
                        d.add(tokens[j]);
                        data.put(labels[j],d);
                    }
                    // check to see if chart of currently selected values exists 
                    xAxisLabel = xGroup.getSelection().getActionCommand();
                    if(selectedLabels.size()==0)
                        break;
                    StringBuffer y = new StringBuffer();
                    for(int j=0;j<selectedLabels.size();j++) {
                        y.append((String)selectedLabels.get(j));
                        if(j!=selectedLabels.size()-1)
                            y.append(",");
                    }
                    String yAxisLabel = y.toString();
                    String graphLabel = "Graph - "+xAxisLabel+"/"+yAxisLabel;
                    alreadyExists = false;
                    int index = -1;
                    for(int i=0;i<main.tpOutput.getTabCount();i++) {
                        String title = main.tpOutput.getTitleAt(i);
                        if(title.equals(graphLabel)) {
                            alreadyExists = true;
                            index = i;
                            break;
                        }
                    }
                    // if chart already exists, add recent data to it
                    if(alreadyExists) {
                        ChartPanel cp = (ChartPanel)main.tpOutput.getComponentAt(index);
                        XYSeriesCollection xyc = (XYSeriesCollection)cp.getChart().getXYPlot().getDataset();
                        for(int i=0;i<xyc.getSeriesCount();i++) {
                            XYSeries s = xyc.getSeries(i);
                            for(int j=0;j<selectedLabels.size();j++) {
                                if(s.getName().equals(selectedLabels.get(j))) {
                                    ArrayList xal = (ArrayList)data.get(xAxisLabel);
                                    ArrayList yal = (ArrayList)data.get(selectedLabels.get(j));
                                    XYDataItem xydi = s.getDataItem(s.getItemCount()-1);
                                    double lastX = xydi.getX().doubleValue();
                                    double currentX = new Double((String)xal.get(xal.size()-1)).doubleValue();
                                    // only add most recent data, do not add if already there
                                    if(lastX<currentX) {
                                        for(int k=0;k<yal.size();k++) {
                                            double xd = new Double((String)xal.get(k)).doubleValue();
                                            if(xd>lastX) {
                                                double yd = new Double((String)yal.get(k)).doubleValue();
                                                s.add(xd,yd);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        main.tpOutput.setSelectedIndex(index);
                    }
                }
                else {
                    throw new Exception("Labels not yet read!");
                }
                xml_entry = null;
            }
            assignChildren(child);
        }
        return 0;
    }
    
    /**
     * Handles selection of data to plot from main graph tab
     * @param e event in graph tab
     */
    public void actionPerformed(ActionEvent e) {
        
        if(e.getActionCommand().equals("Submit")) {
            xAxisLabel = xGroup.getSelection().getActionCommand();
            if(selectedLabels.size()==0)
                return;
            StringBuffer y = new StringBuffer();
            for(int j=0;j<selectedLabels.size();j++) {
                y.append((String)selectedLabels.get(j));
                if(j!=selectedLabels.size()-1)
                    y.append(",");
            }
            String yAxisLabel = y.toString();
            String graphLabel = "Graph - "+xAxisLabel+"/"+yAxisLabel;
            boolean alreadyExists = false;
            int index = -1;
            for(int i=0;i<main.tpOutput.getTabCount();i++) {
                String title = main.tpOutput.getTitleAt(i);
                if(title.equals(graphLabel)) {
                    alreadyExists = true;
                    index = i;
                    break;
                }
            }
            if(alreadyExists)
                main.tpOutput.setSelectedIndex(index);
            else {
                ChartPanel graph = createGraph(graphLabel);
                main.tpOutput.addTab(graphLabel, graph);
                for(int i=0;i<main.tpOutput.getTabCount();i++) {
                    String title = main.tpOutput.getTitleAt(i);
                    if(title.equals(graphLabel)) {
                        index = i;
                        break;
                    }
                }
                main.tpOutput.setSelectedIndex(index);
                createMenuItems(graphLabel);
            }
        }
        
        if(e.getActionCommand().startsWith("closeAction_Graph_")) {
            String name = e.getActionCommand().substring(18);
            for(int i=0;i<main.tpOutput.getTabCount();i++) {
                String title = main.tpOutput.getTitleAt(i);
                if(title.equals(name)) {
                    main.tpOutput.remove(i);
                    destroyMenuItems(name);
                    break;
                }
            }
        }
        
        if(e.getActionCommand().startsWith("printAction_Graph_")) {
            String name = e.getActionCommand().substring(18);
            for(int i=0;i<main.tpOutput.getTabCount();i++) {
                String title = main.tpOutput.getTitleAt(i);
                if(title.equals(name)) {
                    ChartPanel p = (ChartPanel)main.tpOutput.getComponentAt(i);
                    p.createChartPrintJob();
                    break;
                }
            }
        }
        
        if(e.getActionCommand().startsWith("propAction_Graph_")) {
            String name = e.getActionCommand().substring(17);
            for(int i=0;i<main.tpOutput.getTabCount();i++) {
                String title = main.tpOutput.getTitleAt(i);
                if(title.equals(name)) {
                    ChartPanel p = (ChartPanel)main.tpOutput.getComponentAt(i);
                    ChartPropertyEditPanel panel = new ChartPropertyEditPanel(p.getChart());
                    int result = JOptionPane.showConfirmDialog(main.output, panel,
                    "Chart_Properties", JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);
                    if (result == JOptionPane.OK_OPTION) {
                        panel.updateChartProperties(p.getChart());
                    }
                    break;
                }
            }
        }
        
        if(e.getActionCommand().startsWith("saveAction_Graph_")) {
            String name = e.getActionCommand().substring(17);
            for(int i=0;i<main.tpOutput.getTabCount();i++) {
                String title = main.tpOutput.getTitleAt(i);
                if(title.equals(name)) {
                    ChartPanel p = (ChartPanel)main.tpOutput.getComponentAt(i);
                    try {
                        doSaveAs(p);
                    }
                    catch (IOException ex) {
                        main.taMsgs.append("\nUnable to save "+name);
                        main.taMsgs.setCaretPosition(main.taMsgs.getDocument().getLength());
                        ex.printStackTrace();
                    }
                    break;
                }
            }
        }
    }
    
    /**
     * Handle selection/deselection of checkbox
     * @param e event
     */    
    public void itemStateChanged(ItemEvent e) {
        for(int i=0;i<labels.length;i++) {
            if(((JCheckBox)e.getItem()).getText().equals(labels[i])) {
                if(e.getStateChange()==e.SELECTED) {
                    selectedLabels.add(labels[i]);
                }
                else if(e.getStateChange()==e.DESELECTED) {
                    for(int j=0;j<selectedLabels.size();j++) {
                        String s = (String)selectedLabels.get(j);
                        if(s.equals(labels[i]))
                            selectedLabels.remove(j);
                    }
                }
            }
        }
    }
    
    private void createMenuItems(String graphLabel) {
        JMenuItem miSave = new JMenuItem(graphLabel);
        JMenuItem miPrint = new JMenuItem(graphLabel);
        JMenuItem miClose = new JMenuItem(graphLabel);
        JMenuItem miProp = new JMenuItem(graphLabel);
        
        miSave.setMnemonic(graphLabel.charAt(0));
        miSave.setActionCommand("saveAction_Graph_"+graphLabel);
        miSave.addActionListener(this);
        
        miPrint.setMnemonic(graphLabel.charAt(0));
        miPrint.setActionCommand("printAction_Graph_"+graphLabel);
        miPrint.addActionListener(this);
        
        miClose.setMnemonic(graphLabel.charAt(0));
        miClose.setActionCommand("closeAction_Graph_"+graphLabel);
        miClose.addActionListener(this);
        
        miProp.setMnemonic(graphLabel.charAt(0));
        miProp.setActionCommand("propAction_Graph_"+graphLabel);
        miProp.addActionListener(this);
        
        main.menu_save.add(miSave);
        if(!main.menu_print.isEnabled())
            main.menu_print.setEnabled(true);
        main.menu_print.add(miPrint);
        main.menu_close.add(miClose);
        if(!main.menu_prop.isEnabled())
            main.menu_prop.setEnabled(true);
        main.menu_prop.add(miProp);
    }
    
    private void destroyMenuItems(String graphLabel) {
        for(int i=0;i<main.menu_save.getItemCount();i++) {
            JMenuItem mi = main.menu_save.getItem(i);
            if(mi.getText().equals(graphLabel)) {
                main.menu_save.remove(i);
                break;
            }
        }
        for(int i=0;i<main.menu_print.getItemCount();i++) {
            JMenuItem mi = main.menu_print.getItem(i);
            if(mi.getText().equals(graphLabel)) {
                main.menu_print.remove(i);
                break;
            }
        }
        for(int i=0;i<main.menu_close.getItemCount();i++) {
            JMenuItem mi = main.menu_close.getItem(i);
            if(mi.getText().equals(graphLabel)) {
                main.menu_close.remove(i);
                break;
            }
        }
        for(int i=0;i<main.menu_prop.getItemCount();i++) {
            JMenuItem mi = main.menu_prop.getItem(i);
            if(mi.getText().equals(graphLabel)) {
                main.menu_prop.remove(i);
                break;
            }
        }
        if(main.menu_print.getItemCount()==0)
            main.menu_print.setEnabled(false);
        if(main.menu_prop.getItemCount()==0)
            main.menu_prop.setEnabled(false);
    }
    
    /**
     * Opens a file chooser and gives the user an opportunity to save the chart
     * in PNG format.  If file exists, ask the user to replace or select another.
     *
     * @throws IOException if there is an I/O error.
     */
    private void doSaveAs(ChartPanel cp) throws IOException {
        
        ExtensionFileFilter filter1 = new ExtensionFileFilter("PNG_Image_Files", ".png");
        main.fc.addChoosableFileFilter(filter1);
        ExtensionFileFilter filter2 = new ExtensionFileFilter("All files", "");
        main.fc.addChoosableFileFilter(filter2);
        
        File file = null;
        int confirm = 1;
        int option = -1;
        while(confirm==1) {
            option = main.fc.showSaveDialog(cp);
            if (option == JFileChooser.APPROVE_OPTION) {
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
        if (option == JFileChooser.APPROVE_OPTION) {
            String filename = main.fc.getSelectedFile().getPath();
            if (!filename.endsWith(".png")) {
                filename = filename + ".png";
            }
            ChartUtilities.saveChartAsPNG(new File(filename), cp.getChart(), cp.getWidth(), cp.getHeight());
        }
        main.fc.removeChoosableFileFilter(filter1);
        main.fc.removeChoosableFileFilter(filter2);
    }
    
    
}
