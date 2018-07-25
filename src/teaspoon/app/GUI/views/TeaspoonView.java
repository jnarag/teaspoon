/**
 * 
 */
package teaspoon.app.GUI.views;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GraphicsConfiguration;
import java.awt.GridLayout;
import java.awt.HeadlessException;
import java.awt.Toolkit;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.event.ChangeListener;

import teaspoon.app.GUI.models.TeaspoonMaskModel;
import teaspoon.app.GUI.models.TeaspoonModel;
import teaspoon.app.utils.RateEstimationBehaviour;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 5 Dec 2017
 * @version 0.0.1
 */
public class TeaspoonView extends JFrame {

	/*
	 * GUI components
	 */
	private JButton runAnalysis, addMasksTable, removeMasks, combineMasks, removeAlignment, guessDates, selectAncestral, selectBins, showSpectrum;
	private JCheckBox doSlidingWindow, doBootstraps;
	private JLabel maskLabel, fileTableLabel, maskTableLabel, progressLabel;
	private JProgressBar taskBar;
	private JTextField numberBootstraps, slidingWindowSize;
	private JTable filesTable, analysisMasksTable;
	private JPanel mainPanel, maskDisplayPanel, fileListPanel, maskListPanel, tablesPanel, controlsPanel;
	private MaskDisplayPanel maskContentsDisplay;
	private JScrollPane maskPane, maskListPane, filesPane;
	private JSlider bootstrapSlider;
	private JMenuBar menuBar;
	private JMenu menu;
	private JMenuItem menuAbout, menuHelp, menuClear, menuQuit, menuOpen, menuOpenSingle, menuRemoveSingle;
	private JFileChooser setWorkdirLocationChooser, fileChooser;
	static final int BS_MIN = 0;
	static final int BS_MAX = 100;
	static final int BS_INIT = 30;  
	
	/**
	 * Default no-arg constructor
	 * @throws HeadlessException
	 */
	public TeaspoonView() throws HeadlessException {
		super();
		
		// First make the menu
		menuBar = new JMenuBar();
		menu = new JMenu("Menu");
		menuAbout = new JMenuItem("About TEASPOON");
		menuHelp = new JMenuItem("Help");
		menuOpen = new JMenuItem("Open directory");
		menuOpen.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()+1));
		menuOpenSingle = new JMenuItem("Open single");
		menuOpenSingle.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		menuRemoveSingle = new JMenuItem("Remove single");
		menuClear = new JMenuItem("Reset all fields");
		menuQuit = new JMenuItem("Quit TEASPOON");
		menuQuit.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		
		menu.add(menuAbout);
		menu.addSeparator();
		menu.add(menuOpen);
		menu.add(menuOpenSingle);
		menu.add(menuRemoveSingle);
		menu.add(menuClear);
		menu.addSeparator();
		menu.add(menuHelp);
		menu.add(menuQuit);
		
		menuBar.add(menu);


		this.setJMenuBar(menuBar);
		
		
		// Make the buttons etc
		runAnalysis = new JButton("Run analysis");
		addMasksTable = new JButton("Add new mask");
		removeMasks = new JButton("Remove mask");
		combineMasks = new JButton("Combine 2 masks (union)");
		removeAlignment = new JButton("Remove alignment");
		guessDates = new JButton("Guess dates");
		selectAncestral = new JButton("Select ancestral");
		selectBins = new JButton("Select bin intervals");
		showSpectrum = new JButton("Show site-freq spectrum");
		doSlidingWindow = new JCheckBox("Sliding window");
		doBootstraps = new JCheckBox("Bootstrap analysis");
		maskLabel = new JLabel("Alignment mask display");
		fileTableLabel = new JLabel("Alignments");
		maskTableLabel = new JLabel("Masks list");
		numberBootstraps = new JTextField("100");
		slidingWindowSize = new JTextField("50 bp");

		bootstrapSlider = new JSlider(JSlider.HORIZONTAL, BS_MIN, BS_MAX, BS_INIT);
		//Turn on labels at major tick marks.
		bootstrapSlider.setMajorTickSpacing(10);
		bootstrapSlider.setMinorTickSpacing(5);
		bootstrapSlider.setPaintTicks(true);
		bootstrapSlider.setPaintLabels(false);		

		// taskbar
		progressLabel = new JLabel("(inactive)");
		taskBar = new JProgressBar();
		
		controlsPanel = new JPanel();
		controlsPanel.setLayout(new FlowLayout());
		controlsPanel.add(removeAlignment);
		controlsPanel.add(addMasksTable);
		controlsPanel.add(removeMasks);
		controlsPanel.add(combineMasks);
		controlsPanel.add(guessDates);
		controlsPanel.add(selectAncestral);
		controlsPanel.add(selectBins);
		controlsPanel.add(showSpectrum);
		controlsPanel.add(doBootstraps);
		controlsPanel.add(numberBootstraps);
		controlsPanel.add(bootstrapSlider);
		controlsPanel.add(doSlidingWindow);
		controlsPanel.add(slidingWindowSize);
		controlsPanel.add(runAnalysis);
		controlsPanel.add(taskBar);
		controlsPanel.add(progressLabel);

		
		// Make the masks scrollpane
		maskDisplayPanel = new JPanel();
		maskContentsDisplay = new MaskDisplayPanel();
		maskContentsDisplay.setPreferredSize(new Dimension(650,310));
		maskDisplayPanel.setLayout(new GridLayout(2,1));
		maskPane = new JScrollPane(maskContentsDisplay,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		maskPane.setSize(710, 100);
		maskDisplayPanel.add(maskLabel);
		maskDisplayPanel.add(maskPane);
		
		
		// Make the files list
		fileListPanel = new JPanel();
		fileListPanel.setLayout(new GridLayout(2,1));
		filesTable = new JTable();
		//filesTable.setPreferredSize(new Dimension(650,200));
		filesTable.setFillsViewportHeight(true);
		filesTable.setRowSelectionAllowed(true);
		filesTable.setColumnSelectionAllowed(true);
		filesTable.setCellSelectionEnabled(true);
		filesTable.setAutoCreateRowSorter(true);
		filesTable.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
		filesPane = new JScrollPane(filesTable,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		filesPane.setSize(450, 250);
		fileListPanel.add(fileTableLabel);
		fileListPanel.add(filesPane);

		
		// Make the masks list
		maskListPanel = new JPanel();
		maskListPanel.setLayout(new GridLayout(2,1));
		analysisMasksTable = new JTable();
		//filesTable.setPreferredSize(new Dimension(650,200));
		analysisMasksTable.setFillsViewportHeight(true);
		analysisMasksTable.setRowSelectionAllowed(true);
		analysisMasksTable.setColumnSelectionAllowed(true);
		analysisMasksTable.setCellSelectionEnabled(true);
		analysisMasksTable.setAutoCreateRowSorter(true);
		analysisMasksTable.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
		maskListPane = new JScrollPane(analysisMasksTable,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		maskListPane.setSize(350, 200);
		maskListPanel.add(maskTableLabel);
		maskListPanel.add(maskListPane);

		
		// Join the two tables
		tablesPanel = new JPanel();
		tablesPanel.setLayout(new GridLayout(1,2));
		tablesPanel.add(fileListPanel);
		tablesPanel.add(maskListPanel);
		tablesPanel.setSize(700,210);

		
		// Add component sets to main panel
		mainPanel = new JPanel();
		mainPanel.setLayout(new GridLayout(3,1));
		mainPanel.add(maskDisplayPanel);
		mainPanel.add(tablesPanel);
		mainPanel.add(controlsPanel);
		mainPanel.setVisible(true);
		
		add(mainPanel);
		pack();
		setTitle("Teaspoon");
		setDefaultCloseOperation(EXIT_ON_CLOSE);
		setSize(800, 600);
		setVisible(true);

		setWorkdirLocationChooser = new JFileChooser("Select working directory");
		setWorkdirLocationChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
	}

	/**
	 * @param gc
	 */
	public TeaspoonView(GraphicsConfiguration gc) {
		super(gc);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param title
	 * @throws HeadlessException
	 */
	public TeaspoonView(String title) throws HeadlessException {
		super(title);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param title
	 * @param gc
	 */
	public TeaspoonView(String title, GraphicsConfiguration gc) {
		super(title, gc);
		// TODO Auto-generated constructor stub
	}


	/**
	 * Constructs a dialog box to add a new TeaspoonMask.
	 * TODO at the moment there is little error handling for parse issues
	 * See https://stackoverflow.com/questions/1313390/is-there-any-way-to-accept-only-numeric-values-in-a-jtextfield

	 * @param alignmentLength
	 * @return
	 *  Object values:
	 *  [0] RateEstimationBehaviour 
	 *  [1] int start
	 *  [2] int end
	 *  [3] int length
	 *  [4] float ratio
	 */
	public Object[] showTeaspoonMaskDialog(int alignmentLength) {
		// first show a warning since we've called this method without a set length - user will have to input it
		JOptionPane.showMessageDialog(new JFrame(), 
				"You must input integer values for mask_mid start/end/length,"+
				" and floating-point (decimal) for neutral ratio.",
				"Alignment length warning!", JOptionPane.WARNING_MESSAGE);
		/*
		 *  Object values:
		 *  [0] RateEstimationBehaviour 
		 *  [1] int start
		 *  [2] int end
		 *  [3] int length
		 *  [4] float ratio
		 */
		Object[] retArr = {
				RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,
				0,
				0,
				0,
				0.0
		};

		JTextField startField = new JTextField(5);
		startField.setText("0");
		JTextField endField = new JTextField(5);
		endField.setText((alignmentLength-1)+"");
		JLabel lengthField = new JLabel(""+alignmentLength);
		JTextField ratioField = new JTextField(5);
		ratioField.setText("0.0");
		@SuppressWarnings({ "rawtypes", "unchecked" })
		JComboBox methodSelection = new JComboBox(RateEstimationBehaviour.values());

		JPanel myPanel = new JPanel();
		myPanel.add(new JLabel("Method for neutral ratio estimation:"));
		myPanel.add(methodSelection);
		myPanel.add(Box.createHorizontalStrut(15)); // a spacer
		myPanel.add(new JLabel("Start\n(NB: index 0->[k-1]):"));
		myPanel.add(startField);
		myPanel.add(Box.createHorizontalStrut(15)); // a spacer
		myPanel.add(new JLabel("End\n(NB: index 0->[k-1]):"));
		myPanel.add(endField);
		myPanel.add(Box.createHorizontalStrut(15)); // a spacer
		myPanel.add(new JLabel("Length:"));
		myPanel.add(lengthField);
		myPanel.add(Box.createHorizontalStrut(15)); // a spacer
		myPanel.add(new JLabel("Ratio:"));
		myPanel.add(ratioField);

		// show the dialog and get the result
		int result = JOptionPane.showConfirmDialog(null, myPanel, 
				"Please Enter start and end Values", JOptionPane.OK_CANCEL_OPTION);
		// get result and attempt to parse
		if (result == JOptionPane.OK_OPTION) {
			try {
				System.out.println("method value: " + methodSelection.getSelectedIndex());
				retArr[0] = RateEstimationBehaviour.values()[methodSelection.getSelectedIndex()];
				System.out.println("start value: " + startField.getText());
				retArr[1] = Integer.parseInt(startField.getText());
				System.out.println("end value: " + endField.getText());
				retArr[2] = Integer.parseInt(endField.getText());
				System.out.println("length value: " + alignmentLength);
				retArr[3] = alignmentLength;
				System.out.println("length value: " + ratioField.getText());
				retArr[4] = Double.parseDouble(ratioField.getText());
			} catch (NumberFormatException e) {
				JOptionPane.showMessageDialog(new JFrame(), "You must input integer values for mask_mid start/end/length, and floating-point (decimal) for neutral ratio.", "Number Format Error!", JOptionPane.ERROR_MESSAGE);
				e.printStackTrace();
				return null;
			}
		}
		
		return retArr;
	}
	
	/**
	 * Constructs a dialog box to add a new TeaspoonMask.
	 * TODO at the moment there is little error handling for parse issues
	 * See https://stackoverflow.com/questions/1313390/is-there-any-way-to-accept-only-numeric-values-in-a-jtextfield
	 * 
	 * @return
	 *  Object values:
	 *  [0] RateEstimationBehaviour 
	 *  [1] int start
	 *  [2] int end
	 *  [3] int length
	 *  [4] float ratio
	 */
	public Object[] showTeaspoonMaskDialog(){
		// first show a warning since we've called this method without a set length - user will have to input it
		JOptionPane.showMessageDialog(new JFrame(), 
				"You must input integer values for mask_mid start/end/length,"+
				" and floating-point (decimal) for neutral ratio."+
				"\n<b>Don't forget</b> mask_mid start must be < mask_mid end must be <= alignement length!",
				"Alignment length warning!", JOptionPane.WARNING_MESSAGE);
		/*
		 *  Object values:
		 *  [0] RateEstimationBehaviour 
		 *  [1] int start
		 *  [2] int end
		 *  [3] int length
		 *  [4] float ratio
		 */
		Object[] retArr = {
				RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED,
				0,
				0,
				0,
				0.0
		};

		JTextField startField = new JTextField(5);
		startField.setText("0");
		JTextField endField = new JTextField(5);
		endField.setText("0");
		JTextField lengthField = new JTextField(5);
		lengthField.setText("0");
		JTextField ratioField = new JTextField(5);
		ratioField.setText("0.0");
		@SuppressWarnings({ "rawtypes", "unchecked" })
		JComboBox methodSelection = new JComboBox(RateEstimationBehaviour.values());

		JPanel myPanel = new JPanel();
		myPanel.add(new JLabel("Method for neutral ratio estimation:"));
		myPanel.add(methodSelection);
		myPanel.add(Box.createHorizontalStrut(15)); // a spacer
		myPanel.add(new JLabel("Start\n(NB: index 0->[k-1]):"));
		myPanel.add(startField);
		myPanel.add(Box.createHorizontalStrut(15)); // a spacer
		myPanel.add(new JLabel("End\n(NB: index 0->[k-1]):"));
		myPanel.add(endField);
		myPanel.add(Box.createHorizontalStrut(15)); // a spacer
		myPanel.add(new JLabel("Length:"));
		myPanel.add(lengthField);
		myPanel.add(Box.createHorizontalStrut(15)); // a spacer
		myPanel.add(new JLabel("Ratio:"));
		myPanel.add(ratioField);

		// show the dialog and get the result
		int result = JOptionPane.showConfirmDialog(null, myPanel, 
				"Please Enter start and end Values", JOptionPane.OK_CANCEL_OPTION);
		// get result and attempt to parse
		if (result == JOptionPane.OK_OPTION) {
			try {
				System.out.println("method value: " + methodSelection.getSelectedIndex());
				retArr[0] = RateEstimationBehaviour.values()[methodSelection.getSelectedIndex()];
				System.out.println("start value: " + startField.getText());
				retArr[1] = Integer.parseInt(startField.getText());
				System.out.println("end value: " + endField.getText());
				retArr[2] = Integer.parseInt(endField.getText());
				System.out.println("length value: " + lengthField.getText());
				retArr[3] = Integer.parseInt(lengthField.getText());
				System.out.println("length value: " + ratioField.getText());
				retArr[4] = Double.parseDouble(ratioField.getText());
			} catch (NumberFormatException e) {
				JOptionPane.showMessageDialog(new JFrame(), "You must input integer values for mask_mid start/end/length, and floating-point (decimal) for neutral ratio.", "Number Format Error!", JOptionPane.ERROR_MESSAGE);
				e.printStackTrace();
				return null;
			}
		}
		
		return retArr;
	}

	/**
	 * Provides a dialog box to prompt user for a runID. Appends System.currentTimeMillis()
	 * @return
	 */
	public String showTeaspoonRunNameDialog() {

		String returnRunID = System.currentTimeMillis()+"";
		JTextField runIDField = new JTextField(20);
		runIDField.setText(returnRunID);
		JPanel myPanel = new JPanel();
		myPanel.add(new JLabel("Specify an ID for this run:"));
		myPanel.add(runIDField);

		// show the dialog and get the result
		int result = JOptionPane.showConfirmDialog(null, myPanel, 
				"Please Enter a Run ID", JOptionPane.OK_CANCEL_OPTION);
		// get result and attempt to parse
		if (result == JOptionPane.OK_OPTION) {
			try {
				System.out.println("run ID value: " + runIDField.getText());
				returnRunID = runIDField.getText();
			} catch (NumberFormatException e) {
				JOptionPane.showMessageDialog(new JFrame(), "Specify an ID for this run:", "Set run ID...", JOptionPane.ERROR_MESSAGE);
				e.printStackTrace();
				return null;
			}
		}
		return returnRunID;
	}

	/**
	 * @param globalAppModel
	 */
	public void addTables(TeaspoonModel globalAppModel) {
		// TODO Auto-generated method stub
		// FIXME implement
		// this.filesTable.setModel(dataModel);
		this.filesTable.setModel(globalAppModel);
		this.analysisMasksTable.setModel(globalAppModel.getMaskTracksModel());
		this.maskContentsDisplay.setModel(globalAppModel.getMaskTracksModel());
	}

	/**
	 * @param teaspoonCustomGUIaddMaskTrackListener
	 */
	public void addAddMaskListener(ActionListener teaspoonCustomGUIaddMaskTrackListener) {
		// TODO Auto-generated method stub
		this.addMasksTable.addActionListener(teaspoonCustomGUIaddMaskTrackListener);
	}

	/**
	 * @param teaspoonCustomGUIremoveMaskTrackListener
	 */
	public void addRemoveMaskListener(ActionListener teaspoonCustomGUIremoveMaskTrackListener) {
		// TODO Auto-generated method stub
		this.removeMasks.addActionListener(teaspoonCustomGUIremoveMaskTrackListener);
	}

	/**
	 * @param teaspoonCustomGUIcombineMaskTrackListener
	 */
	public void addCombineMaskListener(ActionListener teaspoonCustomGUIcombineMaskTrackListener) {
		// TODO Auto-generated method stub
		this.combineMasks.addActionListener(teaspoonCustomGUIcombineMaskTrackListener);
	}
	
	/**
	 * @param teaspoonCustomGUIaddWorkingDirectoryListener
	 */
	public void addAddWorkingDirectoryListener(ActionListener teaspoonCustomGUIaddWorkingDirectoryListener) {
		// TODO Auto-generated method stub
		this.menuOpen.addActionListener(teaspoonCustomGUIaddWorkingDirectoryListener);
	}

	/**
	 * @param teaspoonCustomGUIaddWorkingDirectoryListener
	 */
	public void addGuessDatesListener(ActionListener teaspoonCustomGUIguessDatesListener) {
		// TODO Auto-generated method stub
		this.guessDates.addActionListener(teaspoonCustomGUIguessDatesListener);
	}

	/**
	 * @param teaspoonCustomGUIclearAllSettingsListener
	 */
	public void addClearAllSettingsListener(ActionListener teaspoonCustomGUIclearAllSettingsListener) {
		// TODO Auto-generated method stub
		this.menuClear.addActionListener(teaspoonCustomGUIclearAllSettingsListener);
	}

	/**
	 * @param teaspoonCustomGUIrunAnalysisListener
	 */
	public void addRunAnalysisListener(ActionListener teaspoonCustomGUIrunAnalysisListener) {
		// TODO Auto-generated method stub
		this.runAnalysis.addActionListener(teaspoonCustomGUIrunAnalysisListener);
	}

	/**
	 * @param teaspoonCustomGUIselectAncestralListener
	 */
	public void addSelectAncestralListener(ActionListener teaspoonCustomGUIselectAncestralListener) {
		// TODO Auto-generated method stub
		this.selectAncestral.addActionListener(teaspoonCustomGUIselectAncestralListener);
	}

	/**
	 * @param teaspoonCustomGUIselectBinsListener
	 */
	public void addSelectBinsListener(ActionListener teaspoonCustomGUIselectBinsListener) {
		// TODO Auto-generated method stub
		this.selectBins.addActionListener(teaspoonCustomGUIselectBinsListener);
	}

	/**
	 * @param teaspoonCustomGUIshowSpectrumListener
	 */
	public void addShowSpectrumListener(ActionListener teaspoonCustomGUIshowSpectrumListener) {
		// TODO Auto-generated method stub
		this.showSpectrum.addActionListener(teaspoonCustomGUIshowSpectrumListener);
	}

	/**
	 * @return
	 */
	public JSlider getBootstrapSlider() {
		// TODO Auto-generated method stub
		return this.bootstrapSlider;
	}

	/**
	 * @param teaspoonCustomGUIBootstrapSliderListener
	 */
	public void addBootstrapSliderListener(ChangeListener teaspoonCustomGUIBootstrapSliderListener) {
		// TODO Auto-generated method stub
		this.bootstrapSlider.addChangeListener(teaspoonCustomGUIBootstrapSliderListener);
	}

	/**
	 * @param teaspoonCustomGUIRemoveAlignmentListener
	 */
	public void addRemoveAlignmentListener(ActionListener teaspoonCustomGUIRemoveAlignmentListener) {
		this.removeAlignment.addActionListener(teaspoonCustomGUIRemoveAlignmentListener);	
	}

	/**
	 * @param bs
	 */
	public void setBootstrapValueDisplay(int bs) {
		// TODO Auto-generated method stub
		this.numberBootstraps.setText(bs+"");
	}

	/**
	 * @return
	 */
	public JFileChooser getDirectoryChooser() {
		// TODO Auto-generated method stub
		return this.setWorkdirLocationChooser;
	}

	/**
	 * @return
	 */
	public JFileChooser getFileChooser() {
		// TODO Auto-generated method stub
		return this.fileChooser;
	}

	/**
	 * @return
	 */
	public JLabel getTaskLabel() {
		// TODO Auto-generated method stub
		return this.progressLabel;
	}

	/**
	 * @return
	 */
	public JProgressBar getTaskBar() {
		// TODO Auto-generated method stub
		return this.taskBar;
	}

	/**
	 * @return
	 */
	public JTable getFilesTable() {
		// TODO Auto-generated method stub
		return this.filesTable;
	}

	/**
	 * @return
	 */
	public JTable getMasksTable() {
		// TODO Auto-generated method stub
		return this.analysisMasksTable;
	}

	/**
	 * Constructs a dialog box to add custom bin intervals.
	 * TODO at the moment there is little error handling for parse issues
	 * See https://stackoverflow.com/questions/1313390/is-there-any-way-to-accept-only-numeric-values-in-a-jtextfield
	 * 
	 * @return
	 *  double[2][3] representing site-frequency bins' intervals
	 */
	public double[][] showCustomBinDialog(double[][] existingBins){
	
		double[][] newBins = existingBins;
		JTextField startLoFreqField = new JTextField(3);
		startLoFreqField.setText(existingBins[0][0]+"");
		JTextField startMidFreqField = new JTextField(3);
		startMidFreqField.setText(existingBins[0][1]+"");
		JTextField startHiFreqField = new JTextField(3);
		startHiFreqField.setText(existingBins[0][2]+"");

		JTextField endLoFreqField = new JTextField(3);
		endLoFreqField.setText(existingBins[1][0]+"");
		JTextField endMidFreqField = new JTextField(3);
		endMidFreqField.setText(existingBins[1][1]+"");
		JTextField endHiFreqField = new JTextField(3);
		endHiFreqField.setText(existingBins[1][2]+"");

		JPanel myPanel = new JPanel();
		myPanel.setLayout(new BoxLayout(myPanel,BoxLayout.Y_AXIS));
		JPanel lowPanel = new JPanel();
		JPanel midPanel = new JPanel();
		JPanel hiPanel = new JPanel();
		
		myPanel.add(Box.createHorizontalStrut(2)); // a spacer
		
		lowPanel.add(new JLabel("Low-freq bin, lower bound, upper:"));
		lowPanel.add(startLoFreqField);
		lowPanel.add(Box.createHorizontalStrut(2)); // a spacer
		lowPanel.add(endLoFreqField);
		myPanel.add(lowPanel);

		myPanel.add(Box.createVerticalStrut(15)); // a spacer
		myPanel.add(Box.createHorizontalStrut(2)); // a spacer

		midPanel.add(new JLabel("Mid-freq bin, lower bound, upper:"));
		midPanel.add(startMidFreqField);
		midPanel.add(Box.createHorizontalStrut(2)); // a spacer
		midPanel.add(endMidFreqField);
		myPanel.add(midPanel);

		myPanel.add(Box.createVerticalStrut(15)); // a spacer
		myPanel.add(Box.createHorizontalStrut(2)); // a spacer

		hiPanel.add(new JLabel("High-freq bin, lower bound, upper:"));
		hiPanel.add(startHiFreqField);
		hiPanel.add(Box.createHorizontalStrut(2)); // a spacer
		hiPanel.add(endHiFreqField);
		myPanel.add(hiPanel);
		
		// show the dialog and get the result
		int result = JOptionPane.showConfirmDialog(null, myPanel, 
				"Please bin start and end Values", JOptionPane.OK_CANCEL_OPTION);
		// get result and attempt to parse
		if (result == JOptionPane.OK_OPTION) {
			try {
				newBins[0][0] = Double.parseDouble(startLoFreqField.getText());
				newBins[0][1] = Double.parseDouble(startMidFreqField.getText());
				newBins[0][2] = Double.parseDouble(startHiFreqField.getText());
				newBins[1][0] = Double.parseDouble(endLoFreqField.getText());
				newBins[1][1] = Double.parseDouble(endMidFreqField.getText());
				newBins[1][2] = Double.parseDouble(endHiFreqField.getText());
			} catch (NumberFormatException e) {
				JOptionPane.showMessageDialog(new JFrame(), "You must input integer values for mask_mid start/end/length, and floating-point (decimal) for neutral ratio.", "Number Format Error!", JOptionPane.ERROR_MESSAGE);
				e.printStackTrace();
				return null;
			}
		}
		
		return newBins;
	}

	/**
	 * @param maskModel
	 */
	public void setMaskPanelModel(TeaspoonMaskModel maskModel) {
		this.maskContentsDisplay.setModel(maskModel);
		
	}

}
