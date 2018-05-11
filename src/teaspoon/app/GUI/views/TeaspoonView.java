/**
 * 
 */
package teaspoon.app.GUI.views;

import java.awt.FlowLayout;
import java.awt.GraphicsConfiguration;
import java.awt.GridLayout;
import java.awt.HeadlessException;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;

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
	private JButton runAnalysis, addMasksTable, removeMasks, guessDates;
	private JCheckBox doSlidingWindow, doBootstraps;
	private JLabel maskLabel, fileTableLabel, maskTableLabel;
	private JTextField numberBootstraps, slidingWindowSize;
	private JTable filesTable, analysisMasksTable;
	private JPanel mainPanel, maskDisplayPanel, fileListPanel, maskListPanel, tablesPanel, controlsPanel;
	private JScrollPane maskPane, maskListPane, filesPane;
	private JMenuBar menuBar;
	private JMenu menu;
	private JMenuItem menuAbout, menuHelp, menuClear, menuQuit, menuOpen;
	private JFileChooser setWorkdirLocationChooser;
	
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
		menuClear = new JMenuItem("Reset all fields");
		menuQuit = new JMenuItem("Quit TEASPOON");
		
		menu.add(menuAbout);
		menu.addSeparator();
		menu.add(menuOpen);
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
		doSlidingWindow = new JCheckBox("Sliding window");
		doBootstraps = new JCheckBox("Bootstrap analysis");
		maskLabel = new JLabel("Alignment mask display");
		fileTableLabel = new JLabel("Alignments");
		maskTableLabel = new JLabel("Masks list");
		numberBootstraps = new JTextField("100 reps");
		slidingWindowSize = new JTextField("50 bp");
		
		controlsPanel = new JPanel();
		controlsPanel.setLayout(new FlowLayout());
		controlsPanel.add(addMasksTable);
		controlsPanel.add(removeMasks);
		controlsPanel.add(guessDates);
		controlsPanel.add(doBootstraps);
		controlsPanel.add(numberBootstraps);
		controlsPanel.add(doSlidingWindow);
		controlsPanel.add(slidingWindowSize);
		controlsPanel.add(runAnalysis);

		
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
		filesPane = new JScrollPane();
		filesPane.setSize(350, 200);
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
		setSize(720, 320);
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

}
