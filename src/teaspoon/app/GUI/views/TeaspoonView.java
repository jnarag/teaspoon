/**
 * 
 */
package teaspoon.app.GUI.views;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GraphicsConfiguration;
import java.awt.GridLayout;
import java.awt.HeadlessException;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ComboBoxModel;
import javax.swing.GroupLayout;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
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
import javax.swing.JSeparator;
import javax.swing.JSlider;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.event.ChangeListener;

import teaspoon.app.TEASPOONVersion;
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

	/**
	 * 
	 * <b>TEASPOON:<b>
	 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
	 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
	 * University of Oxford, 2010-2018.
	 * 
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 2 Aug 2018
	 * @version 0.1
	 */
	public class AboutFrame extends JFrame{

		private static final long serialVersionUID = -4338113095121220945L;
	
		public AboutFrame(){
			super("About TEASPOON (v"+version.getVersion()+")");
			JPanel panel = new JPanel(new FlowLayout());
			/*
			panel.add(new JLabel("Phylogenomic Dataset Browser - alpha version."));
			panel.add(new JLabel("This is a development-only private alpha: use at your own risk."));
			panel.add(new JLabel("(c) Joe Parker / Queen Mary University of London, 2013-5."));
			 */
			panel.add(new JLabel("<html><center>"+version.getHTMLCredits()+"</html>"));
			add(panel);
			setSize(650,700);
			setLocationRelativeTo(null);
			setVisible(true);
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		}
	}

	private static TEASPOONVersion version;
	
	/*
	 * GUI components
	 */
	private JButton runAnalysis, addMasksTable, removeMasks, combineMasks, removeAlignment, guessDates, selectAncestral, selectBins, showSpectrum;
	private JCheckBox doSlidingWindow, doBootstraps;
	private JLabel maskLabel, fileTableLabel, maskTableLabel, progressLabel, numberBootstraps, windowSize;
	private JProgressBar taskBar;
	private JTable filesTable, analysisMasksTable;
	private JPanel mainPanel, maskDisplayPanel, alignmentsPanel, maskPanel, tablesPanel, controlsPanel, alignmentControls, maskControls;
	private MaskDisplayPanel maskContentsDisplay;
	private JScrollPane maskPane, maskListPane, filesPane;
	private JSlider bootstrapSlider, windowSizeSlider;
	private JMenuBar menuBar;
	private JMenu menu, masks, run, window, help;
	private JMenuItem menuClear, menuQuit, menuOpen, menuOpenSingle, menuRemoveSingle, maskMenuAdd, maskMenuDelete, maskMenuCombine, runmenuFastSpectrum, runmenuFullAnalysis;
	public JMenuItem menuAbout, menuHelp, helpOnline, helpBugReport;
	public  JCheckBoxMenuItem windowToggleScatter, windowToggleHistogram;
	private JFileChooser setWorkdirLocationChooser, fileChooser;
	public AboutFrame aboutFrame;
	static final int BS_MIN = 0;
	static final int BS_MAX = 100;
	static final int BS_INIT = 30;  
	
	/**
	 * Default no-arg constructor
	 * @throws HeadlessException
	 */
	public TeaspoonView(TEASPOONVersion versions) throws HeadlessException {
		super();
		
		TeaspoonView.version = versions;
		
		/*-- GUI application menu bar items --*/
		
		// First make the menu and submenus
		menuBar = new JMenuBar();
		menu = new JMenu("TEASPOON"); 	// main application menu
		masks = new JMenu("Masks");		// manipulate masks etc
		run = new JMenu("Run");			// run analyses
		window = new JMenu("Window");	// toggle view options
		help = new JMenu("Help");		// User help
		
		// populate 'TEASPOON' application menu, 
		menuAbout = new JMenuItem("About TEASPOON");
		menuOpen = new JMenuItem("Open directory of alignments");
		menuOpen.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()+1));
		menuOpenSingle = new JMenuItem("Open single alignment");
		menuOpenSingle.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		menuRemoveSingle = new JMenuItem("Remove alignment(s)");
		menuRemoveSingle.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_DELETE, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		menuClear = new JMenuItem("Reset all fields");
		// TODO implement parameter flush
		menuClear.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_R, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask())); // does nothing
		menuQuit = new JMenuItem("Quit TEASPOON");
		menuQuit.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		menuQuit.addActionListener(new ExplicitCloseApplicationListener());
		
		menu.add(menuAbout);
		menu.addSeparator();
		menu.add(menuOpen);
		menu.add(menuOpenSingle);
		menu.add(menuRemoveSingle);
		menu.add(menuClear);
		menu.addSeparator();
		menu.add(menuQuit);

		// populate 'Masks' menu
		maskMenuAdd = new JMenuItem("Add mask track");
		maskMenuAdd.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_M, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		maskMenuDelete = new JMenuItem("Remove selected mask track");
		maskMenuDelete.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_M, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()+1));
		maskMenuCombine = new JMenuItem("Combine two or more masks");
		maskMenuCombine.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		
		masks.add(maskMenuAdd);
		masks.add(maskMenuDelete);
		masks.add(maskMenuCombine);
		
		// populate 'Run' menu
		runmenuFastSpectrum = new JMenuItem("Calculate fast site-frequency spectrum");
		runmenuFastSpectrum.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()+1));
		runmenuFullAnalysis = new JMenuItem("Run full analysis");
		runmenuFullAnalysis.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_E, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		
		run.add(runmenuFastSpectrum);
		run.add(runmenuFullAnalysis);
		
		// populate 'Window' menu
		windowToggleHistogram = new JCheckBoxMenuItem("Show/hide approximate site-frequency plot");
		windowToggleHistogram.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_2, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		windowToggleScatter = new JCheckBoxMenuItem("Show/hide adaptation : time scatterplot");
		windowToggleScatter.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_3, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		
		window.add(windowToggleHistogram);
		window.add(windowToggleScatter);
		
		// populate 'Help' menu 
		menuHelp = new JMenuItem("Help");
		menuHelp.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_H, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		helpOnline = new JMenuItem("Online resources");
		helpBugReport = new JMenuItem("Report a bug...");
		
		help.add(menuHelp);
		help.add(helpOnline);
		help.add(helpBugReport);
		
		// add menus to menu bar
		menuBar.add(menu);
		menuBar.add(masks);
		menuBar.add(run);
		menuBar.add(window);
		menuBar.add(help);


		this.setJMenuBar(menuBar);
		
		
		/*-- GUI Application components for windows --*/
		
		// Make the buttons etc
		// Top (alignment section) group
		addMasksTable = new JButton("Add new mask");
		addMasksTable.setToolTipText("Adds a single mask (analysis region). Use 'TEASPOON>Add Alignments..>Directory' to add whole directories in batch mode.");
		removeMasks = new JButton("Remove mask");
		removeMasks.setToolTipText("Removes a single mask.");
		combineMasks = new JButton("Combine 2 masks (union)");
		combineMasks.setToolTipText("Combines two selected masks into a contiguous analysis region.");

		// Middle (mask section) group
		removeAlignment = new JButton("Remove alignment");
		removeAlignment.setToolTipText("Removes the currently selected alignment. Use 'TEASPOON>Add Alignments..>Directory' to add whole directories in batch mode.");
		guessDates = new JButton("Guess dates");
		guessDates.setToolTipText("Tries to infer dates from the last numeric (decimal) group in the filename. Dates can be edited directly in the Alignments pane above."); 
		selectAncestral = new JButton("Select ancestral");
		selectAncestral.setToolTipText("Defines the currently-selected alignment as the ancestral (outgroup) timepoint. Exactly ONE alignment must be selected.");

		// Bottom (controls section) group
		selectBins = new JButton("Select bin intervals");
		selectBins.setToolTipText("Allows you to define custom boundaries for the Low/Mid/High frequency bins.");
		runAnalysis = new JButton("Run analysis");
		runAnalysis.setToolTipText("Attempts to run an analysis using these parameters and inputs.");
		showSpectrum = new JButton("Show site-freq spectrum");
		showSpectrum.setToolTipText("Calculates an approximate site-frequency spectrum for these alignments (in aggregste mode) using 100 discretised bins from (0->1]."); 
		doSlidingWindow = new JCheckBox("Sliding window size:");
		doSlidingWindow.setToolTipText("(currently disabled feature)");
		doSlidingWindow.setSelected(false);
		doBootstraps = new JCheckBox("Bootstrap replicates:");
		doBootstraps.setSelected(true);
		maskLabel = new JLabel("Alignment masks");
		fileTableLabel = new JLabel("Alignments");
		maskTableLabel = new JLabel("Masks list");
		numberBootstraps = new JLabel("30");
		windowSize = new JLabel("n/a");
		windowSizeSlider = new JSlider(JSlider.HORIZONTAL, 10, 1000, 50);
		windowSizeSlider.setMajorTickSpacing(100);
		windowSizeSlider.setMinorTickSpacing(20);
		windowSizeSlider.setPaintTicks(true);
		windowSizeSlider.setPaintLabels(false);	
		windowSizeSlider.setMaximumSize(new Dimension(300,20));

		bootstrapSlider = new JSlider(JSlider.HORIZONTAL, BS_MIN, BS_MAX, BS_INIT);
		//Turn on labels at major tick marks.
		bootstrapSlider.setMajorTickSpacing(10);
		bootstrapSlider.setMinorTickSpacing(5);
		bootstrapSlider.setPaintTicks(true);
		bootstrapSlider.setPaintLabels(false);	
		bootstrapSlider.setMaximumSize(new Dimension(300,20));

		// taskbar
		progressLabel = new JLabel("(Analysis: inactive)");
		taskBar = new JProgressBar();
		taskBar.setMaximumSize(new Dimension(300,20));
						
		
		// Make the files list
		alignmentsPanel = new JPanel();
		alignmentsPanel.setLayout(new BoxLayout(alignmentsPanel,BoxLayout.PAGE_AXIS));
		alignmentControls = new JPanel();
		alignmentControls.setLayout(new BoxLayout(alignmentControls,BoxLayout.LINE_AXIS));
		filesTable = new JTable();
		//filesTable.setPreferredSize(new Dimension(650,200));
		filesTable.setFillsViewportHeight(true);
		filesTable.setRowSelectionAllowed(true);
		filesTable.setColumnSelectionAllowed(true);
		filesTable.setCellSelectionEnabled(true);
		filesTable.setAutoCreateRowSorter(true);
		filesTable.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
		filesPane = new JScrollPane(filesTable,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		//filesPane.setSize(450, 250);
		alignmentControls.add(removeAlignment);
		alignmentControls.add(guessDates);
		alignmentControls.add(selectAncestral);
		// finish it
		fileTableLabel.setHorizontalAlignment((int) Component.LEFT_ALIGNMENT);
		alignmentsPanel.add(fileTableLabel);
		alignmentsPanel.add(filesPane);
		alignmentsPanel.add(alignmentControls);

		
		// Make the masks list
		maskPanel = new JPanel();
		maskPanel.setLayout(new BoxLayout(maskPanel,BoxLayout.PAGE_AXIS));
		analysisMasksTable = new JTable();
		filesTable.setPreferredSize(new Dimension(650,400));
		analysisMasksTable.setFillsViewportHeight(true);
		analysisMasksTable.setRowSelectionAllowed(true);
		analysisMasksTable.setColumnSelectionAllowed(true);
		analysisMasksTable.setCellSelectionEnabled(true);
		analysisMasksTable.setAutoCreateRowSorter(true);
		analysisMasksTable.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
		maskListPane = new JScrollPane(analysisMasksTable,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		maskListPane.setSize(650, 300);
		//maskTableLabel.setHorizontalAlignment((int) Component.LEFT_ALIGNMENT);
		//maskPanel.add(maskTableLabel);
		maskPanel.add(maskListPane);

		// Make the masks scrollpane and panel
		maskDisplayPanel = new JPanel();
		maskContentsDisplay = new MaskDisplayPanel();
		maskContentsDisplay.setPreferredSize(new Dimension(650,410));
		
		maskContentsDisplay.setSize(650,410);
		maskDisplayPanel.setLayout(new BoxLayout(maskDisplayPanel,BoxLayout.PAGE_AXIS));
		maskPane = new JScrollPane(maskContentsDisplay,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		// Make the mask command buttons panel
		maskControls = new JPanel();
		maskControls.add(addMasksTable);
		maskControls.add(removeMasks);
		maskControls.add(combineMasks);
		maskControls.setLayout(new BoxLayout(maskControls,BoxLayout.LINE_AXIS));
		
		// finish the mask panel off
		maskLabel.setHorizontalAlignment((int) Component.LEFT_ALIGNMENT);
		maskDisplayPanel.add(maskLabel);
		maskDisplayPanel.add(maskPane);
		maskDisplayPanel.add(maskPanel);
		maskDisplayPanel.add(maskControls);
		maskDisplayPanel.setSize(650, 200);
		maskDisplayPanel.setMinimumSize(new Dimension(650,200));
		
		
		// Join the two tables
//		tablesPanel = new JPanel();
//		tablesPanel.setLayout(new GridLayout(1,2));
//		tablesPanel.add(alignmentsPanel);
//		tablesPanel.add(maskPanel);
//		tablesPanel.setSize(700,210);

		controlsPanel = new JPanel();
		//controlsPanel.add(new JSeparator(SwingConstants.HORIZONTAL)); // doesn't seem to work? TODO FIXME
		GroupLayout controlsLayout = new GroupLayout(controlsPanel);
		controlsPanel.setLayout(controlsLayout);
		controlsLayout.setAutoCreateGaps(true);
		controlsLayout.setAutoCreateContainerGaps(true);
		controlsLayout.setHorizontalGroup(
				controlsLayout.createParallelGroup()
				.addGroup(controlsLayout.createSequentialGroup()
						.addComponent(doBootstraps)
						.addComponent(bootstrapSlider)
						.addComponent(numberBootstraps)
						)				
				.addGroup(controlsLayout.createSequentialGroup()
						.addComponent(doSlidingWindow)
						.addComponent(windowSizeSlider)
						.addComponent(windowSize)
						)				
				.addGroup(controlsLayout.createSequentialGroup()
						.addComponent(selectBins)
						.addComponent(showSpectrum)
						.addComponent(runAnalysis)
						)				
				.addGroup(controlsLayout.createSequentialGroup()
						.addComponent(taskBar)
						.addComponent(progressLabel)
						)				
		);
		controlsLayout.setVerticalGroup(
				controlsLayout.createSequentialGroup()
				.addGroup(controlsLayout.createParallelGroup()
						.addComponent(doBootstraps)
						.addComponent(bootstrapSlider)
						.addComponent(numberBootstraps)
						)				
				.addGroup(controlsLayout.createParallelGroup()
						.addComponent(doSlidingWindow)
						.addComponent(windowSizeSlider)
						.addComponent(windowSize)
						)				
				.addGroup(controlsLayout.createParallelGroup()
						.addComponent(selectBins)
						.addComponent(showSpectrum)
						.addComponent(runAnalysis)
						)				
				.addGroup(controlsLayout.createParallelGroup()
						.addComponent(taskBar)
						.addComponent(progressLabel)
						)				
		);
		/*
		controlsPanel.setLayout(new FlowLayout());
		controlsPanel.add(selectBins);
		controlsPanel.add(showSpectrum);
		controlsPanel.add(doBootstraps);
		controlsPanel.add(numberBootstraps);
		controlsPanel.add(bootstrapSlider);
		controlsPanel.add(doSlidingWindow);
		controlsPanel.add(windowSizeSlider);
		controlsPanel.add(runAnalysis);
		controlsPanel.add(taskBar);
		controlsPanel.add(progressLabel);
		 */
		
		// toy icons
		/*
		controlsPanel.add(new JLabel(UIManager.getIcon("OptionPane.informationIcon"))); //boring, java logo in OSX
		controlsPanel.add(new JLabel(UIManager.getIcon("OptionPane.errorIcon"))); // useful
		controlsPanel.add(new JLabel(UIManager.getIcon("OptionPane.warningIcon")));	// useful
		controlsPanel.add(new JLabel(UIManager.getIcon("OptionPane.questionIcon"))); //boring, java logo in OSX
		 */
		
		// Add component sets to main panel
		mainPanel = new JPanel();
		mainPanel.setLayout(new BoxLayout(mainPanel,BoxLayout.PAGE_AXIS));
		mainPanel.add(Box.createRigidArea(new Dimension(10,10)));
		mainPanel.add(alignmentsPanel);
		mainPanel.add(Box.createRigidArea(new Dimension(10,10)));
		mainPanel.add(maskDisplayPanel);
		mainPanel.add(Box.createRigidArea(new Dimension(10,10)));
		mainPanel.add(controlsPanel);
		mainPanel.setVisible(true);
		
		add(mainPanel);
		
		
		/*-- Finalise application rendering and show to user --*/

		this.aboutFrame = new AboutFrame();
		aboutFrame.setVisible(false);
		pack();
		setTitle("Teaspoon: "+version.getVersion());
		setDefaultCloseOperation(EXIT_ON_CLOSE);
		setSize(800, 700);
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
		this.maskMenuAdd.addActionListener(teaspoonCustomGUIaddMaskTrackListener);
	}

	/**
	 * @param teaspoonCustomGUIremoveMaskTrackListener
	 */
	public void addRemoveMaskListener(ActionListener teaspoonCustomGUIremoveMaskTrackListener) {
		// TODO Auto-generated method stub
		this.removeMasks.addActionListener(teaspoonCustomGUIremoveMaskTrackListener);
		this.maskMenuDelete.addActionListener(teaspoonCustomGUIremoveMaskTrackListener);
	}

	/**
	 * @param teaspoonCustomGUIcombineMaskTrackListener
	 */
	public void addCombineMaskListener(ActionListener teaspoonCustomGUIcombineMaskTrackListener) {
		// TODO Auto-generated method stub
		this.combineMasks.addActionListener(teaspoonCustomGUIcombineMaskTrackListener);
		this.maskMenuCombine.addActionListener(teaspoonCustomGUIcombineMaskTrackListener);
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
		this.runmenuFullAnalysis.addActionListener(teaspoonCustomGUIrunAnalysisListener);
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
		this.runmenuFastSpectrum.addActionListener(teaspoonCustomGUIshowSpectrumListener);
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
		this.menuRemoveSingle.addActionListener(teaspoonCustomGUIRemoveAlignmentListener);
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

	/**
	 * @param teaspoonCustomGUItoggleHistogramViewListener
	 */
	public void addToggleHistogramListener(ActionListener teaspoonCustomGUItoggleHistogramViewListener) {
		this.windowToggleHistogram.addActionListener(teaspoonCustomGUItoggleHistogramViewListener);		
	}

	/**
	 * @param teaspoonCustomGUItoggleScatterplotViewListener
	 */
	public void addToggleScatterplotListener(ActionListener teaspoonCustomGUItoggleScatterplotViewListener) {
		this.windowToggleScatter.addActionListener(teaspoonCustomGUItoggleScatterplotViewListener);		
	}

	/**
	 * Shows a generic information dialog with the specified message.
	 * @param string
	 */
	public void showGenericDialog(String message) {
		JPanel myPanel = new JPanel();
		myPanel.add(new JLabel(message));

		// show the dialog and get the result
		JOptionPane.showConfirmDialog(null, myPanel,"Information", JOptionPane.OK_CANCEL_OPTION);	
	}
	
	/**
	 * <b>TEASPOON:<b>
	 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
	 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
	 * University of Oxford, 2010-2018.
	 * 
	 * Teensy private listener to close the application for the JMenu 'quit' option.
	 * 
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 6 Aug 2018
	 * @version 0.1
	 */
	private class ExplicitCloseApplicationListener implements ActionListener{
		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			Runtime.getRuntime().exit(EXIT_ON_CLOSE);
		}
	}
}
