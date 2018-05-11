package teaspoon.app.GUI.controllers;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Date;

import teaspoon.app.GUI.models.TeaspoonModel;
import teaspoon.app.GUI.views.TeaspoonView;

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
public class TeaspoonController {

	TeaspoonView appView;
	TeaspoonModel appModel;
	
	/**
	 * No-arg constructor. Deprecated
	 */
	@Deprecated
	public TeaspoonController(){}
	
	/**
	 * Main contstructor, requires a view and model
	 * @param appView
	 * @param appModel
	 */
	public TeaspoonController(TeaspoonView globalAppView, TeaspoonModel globalAppModel) {
		// basic setup
		this.appView = globalAppView;
		this.appModel = globalAppModel;
		
		// now the complicated bits - first bind tables to view
		this.appView.addTables(globalAppModel);
		
		// now add actionListeners to GUI view
		this.appView.addAddMaskListener(new TeaspoonCustomGUIaddMaskTrackListener());
		this.appView.addRemoveMaskListener(new TeaspoonCustomGUIremoveMaskTrackListener());
		this.appView.addAddWorkingDirectoryListener(new TeaspoonCustomGUIaddWorkingDirectoryListener());
		this.appView.addGuessDatesListener(new TeaspoonCustomGUIguessDatesListener());
		this.appView.addClearAllSettingsListener(new TeaspoonCustomGUIclearAllSettingsListener());
		this.appView.addRunAnalysisListener(new TeaspoonCustomGUIrunAnalysisListener());
	}

	/**
	 * Opens a directory, crawls for FASTA and parses them into the files list
	 */
	public void addDirectory(){
		// FIXME implement
		guessDates();
	}
	
	/**
	 * Guesses the date for each alignment file
	 * @return - Date[] of each alignment
	 */
	public Date[] guessDates(){
		// FIXME implement
		return null;
	}
	
	/**
	 * The user sets a date for a single alignment manually
	 */
	public void userSetDate(){
		// FIXME implement
	}
	
	/**
	 * Runs analysis via masker/CLI
	 */
	public void runAnalysis(){
		// FIXME implement
		// don't forget to take bootstrap / sliding-window settings etc.
	}
	
	// ------------------------------------------
	
	/*
	 * Section for a load of GUI action listeners
	 */
	
	/**
	 * Need custom listeners
	 * FIXME implement
	 * <b>TEASPOON:<b>
	 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
	 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
	 * University of Oxford, 2010-2018.
	 * 
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 10 May 2018
	 * @version 0.1
	 * 
	 * Add a new mask track to the analysis.
	 */
	private class TeaspoonCustomGUIaddMaskTrackListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			System.out.println("Action Event: Add a mask track");
		}		
	}

	/**
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 11 May 2018
	 * @version 0.1
	 * 
	 * Remove the selected analysis masking track.
	 */
	private class TeaspoonCustomGUIremoveMaskTrackListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			System.out.println("Action Event: Remove a mask track");
		}		
	}

	/**
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 11 May 2018
	 * @version 0.1
	 * 
	 * Guess alignment file sampling dates
	 */
	private class TeaspoonCustomGUIguessDatesListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			System.out.println("Action Event: Guess alignment dates.");
		}		
	}

	/**
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 11 May 2018
	 * @version 0.1
	 * 
	 * Select a new working directory.
	 */
	private class TeaspoonCustomGUIaddWorkingDirectoryListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			System.out.println("Action Event: Add a working directory.");
		}		
	}

	/**
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 11 May 2018
	 * @version 0.1
	 * 
	 * Clears all settings from the GUI.
	 */
	private class TeaspoonCustomGUIclearAllSettingsListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			System.out.println("Action Event: Clear all settings");
		}		
	}
	
	/**
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 11 May 2018
	 * @version 0.1
	 * 
	 * Runs an analysis with the specified settings.
	 */
	private class TeaspoonCustomGUIrunAnalysisListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			System.out.println("Action Event: Run analysis");
			runAnalysis();
		}		
	}
}
