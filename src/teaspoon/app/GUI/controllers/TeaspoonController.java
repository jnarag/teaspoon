package teaspoon.app.GUI.controllers;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.JSlider;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import teaspoon.app.TeaspoonMask;
import teaspoon.app.GUI.models.TeaspoonMaskModel;
import teaspoon.app.GUI.models.TeaspoonModel;
import teaspoon.app.GUI.views.TeaspoonView;
import teaspoon.app.standalone.TeaspoonCommandLineApp;
import teaspoon.app.standalone.TeaspoonMaskFactory;
import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.MainAlignmentParser;
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
public class TeaspoonController {

	TeaspoonView appView;
	TeaspoonModel appModel;
	TeaspoonMaskModel maskModel;
	
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
		this.appView.addSelectAncestralListener(new TeaspoonCustomGUIselectAncestralListener());
		this.appView.addSelectBinsListener(new TeaspoonCustomGUIselectBinsListener());
		this.appView.addBootstrapSliderListener(new TeaspoonCustomGUIBootstrapSliderListener());
		this.appView.addRemoveAlignmentListener(new TeaspoonCustomGUIRemoveAlignmentListener());
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
	public void runAnalysis(BhattAdaptationParameters parameters){
		// FIXME implement
		// don't forget to take bootstrap / sliding-window settings etc.
		try {
			new TeaspoonCommandLineApp(parameters);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
			if(appModel.getAlignmentLength()>0 && appModel.hasAncestralSequenceAlignmentBeenSet()){
				Object[] getMaskSpecification = appView.showTeaspoonMaskDialog(appModel.getAlignmentLength());
				if(getMaskSpecification != null){
					RateEstimationBehaviour behaviour = (RateEstimationBehaviour) getMaskSpecification[0];
					int start = (int) getMaskSpecification[1];
					int end = (int) getMaskSpecification[2];
					int length = (int) getMaskSpecification[3];
					double ratio = (double) getMaskSpecification[4];
					TeaspoonMask newMaskTrack;
					if(behaviour == RateEstimationBehaviour.NEUTRAL_RATE_FIXED){
						newMaskTrack = TeaspoonMaskFactory.initialiseMask(ratio, start, end, length);
					}else{
						newMaskTrack = TeaspoonMaskFactory.initialiseMask(behaviour, start, end, length);
					}
					appModel.addMaskRow(newMaskTrack);
				}
				
			}else{
				Object[] getMaskSpecification = appView.showTeaspoonMaskDialog();
				if(getMaskSpecification != null){
					RateEstimationBehaviour behaviour = (RateEstimationBehaviour) getMaskSpecification[0];
					int start = (int) getMaskSpecification[1];
					int end = (int) getMaskSpecification[2];
					int length = (int) getMaskSpecification[3];
					double ratio = (double) getMaskSpecification[4];
					TeaspoonMask newMaskTrack;
					if(behaviour == RateEstimationBehaviour.NEUTRAL_RATE_FIXED){
						newMaskTrack = TeaspoonMaskFactory.initialiseMask(ratio, start, end, length);
					}else{
						newMaskTrack = TeaspoonMaskFactory.initialiseMask(behaviour, start, end, length);
					}
					appModel.addMaskRow(newMaskTrack);
				}
				
			}
		}		
	}

	/**
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 11 May 2018
	 * @version 0.1
	 * 
	 * Remove the selected analysis masking track.
	 */
	private class TeaspoonCustomGUIRemoveAlignmentListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			System.out.println("Action Event: Remove an input alignment");
			/*
			 * 1. check exactly one row is selected
			 * 2. toggle this row's ancestral property
			 * 3. add ancestral to model data 
			 */
			try{
				if(appView.getFilesTable().getSelectedRows().length == 1){
					// toggle state
					appModel.removeFileWithSelectedRow(appView.getFilesTable().getSelectedRow());
				}else{
					// invalid
					throw new Exception();
				}
			}catch (Exception ex){
				System.err.println("exactly one alignment must be selected");
				ex.printStackTrace();
			}
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
			/*
			 * 1. check exactly one row is selected
			 * 2. toggle this row's ancestral property
			 * 3. add ancestral to model data 
			 */
			try{
				if(appView.getMasksTable().getSelectedRows().length == 1){
					// toggle state
					appModel.removeMaskWithSelectedRow(appView.getMasksTable().getSelectedRow());
				}else{
					// invalid
					throw new Exception();
				}
			}catch (Exception ex){
				System.err.println("exactly one alignment must be selected");
				ex.printStackTrace();
			}
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
	
	private class TeaspoonCustomGUIselectAncestralListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			System.out.println("Action Event: select ancestral alignment");
			/*
			 * 1. check exactly one row is selected
			 * 2. toggle this row's ancestral property
			 * 3. add ancestral to model data 
			 */
			try{
				if(appView.getFilesTable().getSelectedRows().length == 1){
					// toggle state
					appModel.updateAncestralWithSelectedRow(appView.getFilesTable().getSelectedRow());
				}else{
					// invalid
					throw new Exception();
				}
			}catch (Exception ex){
				System.err.println("exactly one alignment must be selected");
				ex.printStackTrace();
			}
		}		
	}

	/**
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 11 May 2018
	 * @version 0.1
	 * 
	 * Select a new working directory.
	 */
	private class TeaspoonCustomGUIaddWorkingDirectoryListener implements ActionListener, PropertyChangeListener{
		AlignmentsImportTask task;
		JLabel taskLabel;
		JProgressBar taskBar;
		public final String completeText = "Done ";
		public final int completeInt = 100;

	    /**
	     * Invoked when task's progress property changes.
	     */
	    public void propertyChange(PropertyChangeEvent evt) {
	        if ("progress" == evt.getPropertyName()) {
	        	int progress = task.getProgress();
	            String message = "Adding alignments ("+progress+";%)...";
	            taskLabel.setText(message);
	            taskBar.setValue(progress);
	        } 
	    }

	    /* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent ev) {
			System.out.println("Action Event: Add a working directory.");
			int returnVal = appView.getDirectoryChooser().showOpenDialog(appView);
			if(returnVal == JFileChooser.APPROVE_OPTION){
				File alignmentDirectory = appView.getDirectoryChooser().getSelectedFile();
				if(alignmentDirectory.isDirectory()){
					taskLabel = appView.getTaskLabel();
					taskBar = appView.getTaskBar();
					task = new AlignmentsImportTask(taskLabel, taskBar);
			        task.addPropertyChangeListener((PropertyChangeListener) this);
			        task.execute();
			        taskLabel.setText(completeText);
			        taskBar.setValue(completeInt);

			        /*
			         * Non-swingworker, non-progress-bar method:
			         * 
					File[] files = alignmentDirectory.listFiles();
					for(File file:files){
						if(file.exists()){
							System.out.println(file);
							appModel.addFileRow(file);
						}
					}
			         */
				}
			}
		}
	}

	/**
	 * 
	 * <b>TEASPOON:<b>
	 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
	 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
	 * University of Oxford, 2010-2018.
	 * 
	 * Get the bootstrap count from slider 
	 * 
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 30 Jun 2018
	 * @version 0.1
	 */
	private class TeaspoonCustomGUIBootstrapSliderListener implements ChangeListener{
		
		@Override
		public void stateChanged(ChangeEvent e) {
			JSlider source = (JSlider)e.getSource();
			if (!source.getValueIsAdjusting()) {
				int bs = (int)source.getValue();
				appModel.setBootstraps(bs);
				appView.setBootstrapValueDisplay(bs);
				System.out.println("Update bootstraps: "+bs);
			}
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
			if(!appModel.hasAncestralSequenceAlignmentBeenSet()){
				// there is no ancestral alignment set, unlikely we can proceed
				JOptionPane.showMessageDialog(new JFrame(), "No sequence selected as ancestral!", "Analysis Specification Error!", JOptionPane.ERROR_MESSAGE);
			}else if(appModel.getMainSequenceAlignments().length<1){
				// there is no main sequence alignment set, unlikely we can proceed
				JOptionPane.showMessageDialog(new JFrame(), "No main alignments available!", "Analysis Specification Error!", JOptionPane.ERROR_MESSAGE);
			}else if(!appModel.areAlignmentsOfEqualLength()){
				// main alignments not same lengths, unlikely we can proceed
				JOptionPane.showMessageDialog(new JFrame(), "Alignments should all be same length", "Analysis Specification Error!", JOptionPane.ERROR_MESSAGE);
			}else{
				//good to go...
				System.out.println("Action Event: Run analysis");
				// set run prefix and name
				String runID = appView.showTeaspoonRunNameDialog();
				//double ratio = 0.7186788;
				File maskFile = new File("./HCV_data/sub_053/"+runID+".mask");
				//File ancestral = new File("./HCV_data/sub_053/FP7_05301_0.fasta");
				File output = new File("./HCV_data/sub_053/"+runID+".out");
				/*
				File[] inputList = new File[8];
				inputList[0] = new File("./HCV_data/sub_053/FP7_05302_0.3644.fasta");
				inputList[1] = new File("./HCV_data/sub_053/FP7_05303_0.6137.fasta");
				inputList[2] = new File("./HCV_data/sub_053/FP7_05304_0.8438.fasta");
				inputList[3] = new File("./HCV_data/sub_053/FP7_05305_1.3699.fasta");
				inputList[4] = new File("./HCV_data/sub_053/FP7_05306_1.7836.fasta");
				inputList[5] = new File("./HCV_data/sub_053/FP7_05307_3.8986.fasta");
				inputList[6] = new File("./HCV_data/sub_053/FP7_05308_6.8429.fasta");
				inputList[7] = new File("./HCV_data/sub_053/FP7_05309_7.6849.fasta");
				BhattAdaptationFullSiteMatrix alignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(ancestral).readFASTA());
				ArrayList<TeaspoonMask> masks = new ArrayList<TeaspoonMask>();
				int[] maskStartEnd = {0, alignment.alignmentLength()-1};
				ArrayList<int[]> maskRanges = new ArrayList<int[]>();
				maskRanges.add(maskStartEnd);
				masks.add(TeaspoonMaskFactory.initialiseMask(ratio, 0, alignment.alignmentLength()-1, alignment.alignmentLength()));
				*/
				// try programmatic set
				ArrayList<TeaspoonMask> masks = appModel.getMaskTracks();
				try {
					TeaspoonMaskFactory.writeMaskFile(maskFile, masks);
				} catch (IOException ex) {
					// TODO Auto-generated catch block
					ex.printStackTrace();
				}
				BhattAdaptationParameters parameters = appModel.getParametersSnapshot();
				try {
					// hardcoded here
					parameters.setMaskFile(maskFile);
					parameters.setOutputFile(output);
					/*
					 * From GUI model directly:
					 * parameters.setBootstrapReplicates(numBootstraps);
					 * parameters.setAncestralFile(ancestral);
					 * parameters.setInputFileList(inputList);
					 * parameters.setNeutralRate(ratio);
					 */
				} catch (IOException ioEx) {
					// TODO Auto-generated catch block
					ioEx.printStackTrace();
				}
				runAnalysis(parameters);
			}
		}		
	}

	class AlignmentsImportTask extends SwingWorker<Void, Void> {
		JLabel taskLabel;
		JProgressBar taskBar;
		File forceOpen = null;
		public final String completeText = "Done ";
		public final int completeInt = 100;
		
		/**
		 * Set an optional tasklabel/bar
		 * @param label
		 * @param bar
		 */
		public AlignmentsImportTask(JLabel label, JProgressBar bar){
			taskLabel = label;
			taskBar = bar;
		}
	
		/**
		 * No-arg constructor - risky since taskLabel and taskBar will not be instantiated.
		 */
		public AlignmentsImportTask() {
			// TODO Auto-generated constructor stub
		}
	
		/**
		 * Force the task to open a certain File
		 * @param file
		 */
		public void setForceOpen(File file){
			forceOpen = file;
		}
		/*
	     * Main task for adding alignment files. Executed in background thread.
	     */
	    @Override
	    public Void doInBackground() {
	        int progress = 0;
	        //Initialize progress property.
	        setProgress(0);
	        File[] files;
	        if(forceOpen != null && forceOpen.isDirectory()){
	        	files = forceOpen.listFiles();
	        }else{
				files  = appView.getDirectoryChooser().getSelectedFile().listFiles();
	        }
			int totalFiles = files.length;
			int filesTried = 0;
			for(File alignmentFile:files){
				// Attempt to check we have a valid filename.
				if(
					alignmentFile.exists() && 
					!alignmentFile.isDirectory() &&
					( 	alignmentFile.getName().endsWith("fa") || 
						alignmentFile.getName().endsWith("fasta")	
					) 
					){
					System.out.println(alignmentFile+"File looks like .fa or .fasta from filename. Parsing...: ");
					try {
						appModel.addFileRow(alignmentFile);
					} catch (NullPointerException ex) {
						ex.printStackTrace();
					}
				}else{
					System.out.println("File doesn't look like .fa or .fasta from filename, skipping: "+alignmentFile);
				}
				filesTried++;
				progress = Math.round(((float)filesTried / (float)totalFiles)*100f);
				String message = "Adding alignments ("+filesTried+"; "+progress+"%)...";
				taskLabel.setText(message);
				taskBar.setValue(progress);
	            setProgress(Math.min(progress, 100));
				System.out.println("Adding alignments ("+filesTried+"; "+progress+"%)...");
			}
	      return null;
	    }
	
	    /*
	     * Executed in event dispatching thread
	     */
	    @Override
	    public void done() {
	    	taskLabel.setText(completeText);
	    	taskBar.setValue(completeInt);
	    }
	}

	/**
	 * <b>TEASPOON:<b>
	 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
	 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
	 * University of Oxford, 2010-2018.
	 * 
	 * Sets custom site-frequency bin intervals (double[2][3])
	 * Bin intervals may or may not be overlapping, and from the range [0,1] in any case
	 * But there must be exactly 3, and the middle (ie index[1], not necessarily
	 * numerically central) bin will be presumed to be the 'mid-frequency' 
	 * substitutions - used in estimation for the neutral ratio.
	 * 
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 11 Jul 2018
	 * @version 0.1
	 */
	private class TeaspoonCustomGUIselectBinsListener implements ActionListener{
	
		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			System.out.println("Action Event: select bin intervals");
			/*
			 * 1. check exactly one row is selected
			 * 2. toggle this row's ancestral property
			 * 3. add ancestral to model data 
			 */

			double[][] defaultSensibleBins = {
						{0.0,0.15,0.75},
						{0.15,0.75,1.0}
				};
			appModel.setCustomBinIntervals(
				appView.showCustomBinDialog(defaultSensibleBins)
			);
		}		
	}
}
