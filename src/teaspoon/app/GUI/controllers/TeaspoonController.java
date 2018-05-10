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
		this.appView = globalAppView;
		this.appModel = globalAppModel;
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
	 */
	private class TeaspoonCustomGUIListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub			
		}		
	}
}
