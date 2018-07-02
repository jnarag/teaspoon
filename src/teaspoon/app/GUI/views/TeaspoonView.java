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

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.event.ChangeListener;

import teaspoon.app.GUI.models.TeaspoonModel;

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
	private JButton runAnalysis, addMasksTable, removeMasks, guessDates, selectAncestral;
	private JCheckBox doSlidingWindow, doBootstraps;
	private JLabel maskLabel, fileTableLabel, maskTableLabel, progressLabel;
	private JProgressBar taskBar;
	private JTextField numberBootstraps, slidingWindowSize;
	private JTable filesTable, analysisMasksTable;
	private JPanel mainPanel, maskDisplayPanel, fileListPanel, maskListPanel, tablesPanel, controlsPanel;
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
		guessDates = new JButton("Guess dates");
		selectAncestral = new JButton("Select ancestral");
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
		controlsPanel.add(addMasksTable);
		controlsPanel.add(removeMasks);
		controlsPanel.add(guessDates);
		controlsPanel.add(selectAncestral);
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
		maskDisplayPanel.setLayout(new GridLayout(2,1));
		maskPane = new JScrollPane();
		maskPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		maskPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
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
		maskListPane = new JScrollPane();
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
	 * @param globalAppModel
	 */
	public void addTables(TeaspoonModel globalAppModel) {
		// TODO Auto-generated method stub
		// FIXME implement
		// this.filesTable.setModel(dataModel);
		this.filesTable.setModel(globalAppModel);
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

}
