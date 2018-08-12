package teaspoon.app.GUI.controllers;

import java.awt.Desktop;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.JSlider;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import teaspoon.app.TEASPOONVersion;
import teaspoon.app.TeaspoonMask;
import teaspoon.app.GUI.models.TeaspoonMaskModel;
import teaspoon.app.GUI.models.TeaspoonModel;
import teaspoon.app.GUI.views.SimpleHistogramPlottingFrame;
import teaspoon.app.GUI.views.SimpleRegressionPlottingFrame;
import teaspoon.app.GUI.views.TeaspoonView;
import teaspoon.app.standalone.TeaspoonCommandLineApp;
import teaspoon.app.standalone.TeaspoonFastSiteFreqApp;
import teaspoon.app.standalone.TeaspoonMaskFactory;
import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.BhattAdaptationResults;
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

	private static TEASPOONVersion version;
	TeaspoonView appView;
	TeaspoonModel appModel;
	SimpleRegressionPlottingFrame plotter;
	SimpleHistogramPlottingFrame histogram;

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
	public TeaspoonController(TeaspoonView globalAppView, TeaspoonModel globalAppModel, TEASPOONVersion version) {
		// basic setup
		TeaspoonController.version = version;
		this.appView = globalAppView;
		this.appModel = globalAppModel;
		this.plotter = new SimpleRegressionPlottingFrame();
		this.histogram = new SimpleHistogramPlottingFrame();

		// now the complicated bits - first bind tables to view
		this.appView.addTables(globalAppModel);

		// now add actionListeners to GUI view
		this.appView.addAddMaskListener(new TeaspoonCustomGUIaddMaskTrackListener());
		this.appView.addRemoveMaskListener(new TeaspoonCustomGUIremoveMaskTrackListener());
		this.appView.addCombineMaskListener(new TeaspoonCustomGUIcombineMaskTrackListener());
		this.appView.addAddWorkingDirectoryListener(new TeaspoonCustomGUIaddWorkingDirectoryListener());
		this.appView.addGuessDatesListener(new TeaspoonCustomGUIguessDatesListener());
		this.appView.addClearAllSettingsListener(new TeaspoonCustomGUIclearAllSettingsListener());
		this.appView.addRunAnalysisListener(new TeaspoonCustomGUIrunAnalysisListener());
		this.appView.addSelectAncestralListener(new TeaspoonCustomGUIselectAncestralListener());
		this.appView.addSelectBinsListener(new TeaspoonCustomGUIselectBinsListener());
		this.appView.addShowSpectrumListener(new TeaspoonCustomGUIshowSpectrumListener());
		this.appView.addBootstrapSliderListener(new TeaspoonCustomGUIBootstrapSliderListener());
		this.appView.addRemoveAlignmentListener(new TeaspoonCustomGUIRemoveAlignmentListener());
		this.appView.addToggleHistogramListener(new TeaspoonCustomGUItoggleHistogramViewListener());
		this.appView.addToggleScatterplotListener(new TeaspoonCustomGUItoggleScatterplotViewListener());

		// add actionlisteners for menu items
		this.appView.menuAbout.addActionListener(new AboutMenuListener());
		this.appView.menuHelp.addActionListener(new OpenURLListener("https://github.com/jnarag/teaspoon/blob/master/Usage.md"));
		this.appView.helpBugReport.addActionListener(new OpenURLListener("https://github.com/jnarag/teaspoon/issues"));
		this.appView.helpOnline.addActionListener(new OpenURLListener("https://github.com/jnarag/teaspoon"));
}

	public void fitAndShowRegression(HashMap<File,BhattAdaptationResults> results, String name){
		if(this.appModel != null){
			HashMap<File,Float> inputDates = appModel.getDates();
			Iterator<File> resultsIterator = results.keySet().iterator();
			ArrayList<Float[]> points = new ArrayList<Float[]>();
			while(resultsIterator.hasNext()){
				File mainFile = resultsIterator.next();
				if(inputDates.containsKey(mainFile)){
					BhattAdaptationResults mainFileResult = results.get(mainFile);
					double nonNeutralAdaptations = mainFileResult.getBhattSiteCounter().getNonNeutralSubstitutions()[2];
					float  inputFileDate = inputDates.get(mainFile);
					// check for bootstrap results
					if(!mainFileResult.hasBootstraps){
						// no bootstraps - just send the list of {x, y} points to regression plotting frame
						points.add(new Float[]{
								inputFileDate, 
								(float)nonNeutralAdaptations
								});
						System.out.println("results returned from "+mainFile+" :: "+nonNeutralAdaptations+"\t"+inputFileDate);
					}else{
						// print the BS results to console and send the {x, y, (y_error = y_hi-y_low)} points to regression plotting frame
						System.out.println("\t"+nonNeutralAdaptations+" bs lo "+mainFileResult.getBootstrapAdaptationEstimates().getMin());
						System.out.println("\t"+nonNeutralAdaptations+" bs mean "+mainFileResult.getBootstrapAdaptationEstimates().getMean());
						System.out.println("\t"+nonNeutralAdaptations+" bs median "+mainFileResult.getBootstrapAdaptationEstimates().getPercentile(50));
						System.out.println("\t"+nonNeutralAdaptations+" bs hi "+mainFileResult.getBootstrapAdaptationEstimates().getMax());
						points.add(new Float[]{
								inputFileDate, 
								(float)nonNeutralAdaptations,
								(float)mainFileResult.getBootstrapAdaptationEstimates().getPercentile(97) - (float)mainFileResult.getBootstrapAdaptationEstimates().getPercentile(3)
								});
						System.out.println("results returned from "+mainFile+" :: "+nonNeutralAdaptations+"\t"+inputFileDate);				
					}
				}
			} 

			// show the plot
			this.plotter.updateScatterChart(name, points);
			toggleRegressionPlotVisible();
		}
	}

	/**
	 * Runs analysis via masker/CLI
	 */
	public void runAnalysis(BhattAdaptationParameters parameters){
		// FIXME implement
		// don't forget to take bootstrap / sliding-window settings etc.
		try {
			TeaspoonCommandLineApp analysis = new TeaspoonCommandLineApp(parameters);
			HashMap<File,BhattAdaptationResults> results = analysis.getResults();
			fitAndShowRegression(results, parameters.getRunID());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Runs analysis via masker/CLI as above, but overloaded to accept UI task as well
	 */
	public void runAnalysis(BhattAdaptationParameters parameters, MainBhattAdaptationAnalysisTask task){
		// FIXME implement
		// don't forget to take bootstrap / sliding-window settings etc.
		try {
			TeaspoonCommandLineApp analysis = new TeaspoonCommandLineApp(parameters, task);
			HashMap<File,BhattAdaptationResults> results = analysis.getResults();
			fitAndShowRegression(results, parameters.getRunID());
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
	 * Runs analysis via masker/CLI
	 * @param siteFreqPlottingTask 
	 * @since 24 Jul 2018
	 */
	public void runFastSiteFreqAnalysis(BhattAdaptationParameters parameters, int numBins, SiteFreqPlottingTask siteFreqPlottingTask){
		// FIXME implement
		// don't forget to take bootstrap / sliding-window settings etc.
		try {
			TeaspoonFastSiteFreqApp analysis = new TeaspoonFastSiteFreqApp(parameters, numBins, siteFreqPlottingTask);
			HashMap<File,float[][]> siteFreqSpectrum = analysis.getResults();
			float[][] spectrum = siteFreqSpectrum.get(new File("fast"));
			ArrayList<Float[]> spectrumPlottingData = new ArrayList<Float[]>();
			for(float[] bin:spectrum){
				Float[] binAsFloat = new Float[2];
				binAsFloat[0] = bin[0];
				binAsFloat[1] = bin[1];
				spectrumPlottingData.add(binAsFloat);
			}
			histogram.updateHistogram("SiteFreq", spectrumPlottingData);
			//histogram.setVisible(true);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * toggles site-frequency historgram visibility on/off
	 */
	private void toggleSiteFreqHistogramVisible() {
		if(histogram.isVisible()){
			histogram.setVisible(false);
			appView.windowToggleHistogram.setSelected(false);
		}else{
			histogram.setVisible(true);
			appView.windowToggleHistogram.setSelected(true);

		}
	}

	/**
	 * toggles regression plot visibility on/off
	 */
	private void toggleRegressionPlotVisible() {
		if(plotter.isVisible()){
			plotter.setVisible(false);
		}else{
			plotter.setVisible(true);
			appView.windowToggleScatter.setSelected(true);
		}	
	}
	
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
	 * Add a new mask_mid track to the analysis.
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
					appView.repaint();
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
			/*
			 * 1. check exactly one row is selected
			 * 2. toggle this row's ancestral property
			 * 3. add ancestral to model data 
			 */
			try{
				if(appView.getFilesTable().getSelectedRows().length == 1){
					// toggle state
					System.out.println("Action Event: Remove an input alignment (row "+appView.getFilesTable().getSelectedRows()[0]+")");
					appModel.removeFileWithSelectedRow(appView.getFilesTable().getSelectedRow());
					appModel.fireTableDataChanged();
					appView.repaint();
				}else if(appView.getFilesTable().getSelectedRows().length > 1){
					// more than one row
					System.out.print("Action Event: Remove "+
							appView.getFilesTable().getSelectedRows().length+
							" input alignment (rows ");
					for(int row:appView.getFilesTable().getSelectedRows()){
						System.out.print(row+", ");
					}
					System.out.println(")");
					appModel.removeFileWithSelectedRows(appView.getFilesTable().getSelectedRows());
					appModel.fireTableDataChanged();
					appView.repaint();
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
					appView.repaint();
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
	 * <b>TEASPOON:<b>
	 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
	 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
	 * University of Oxford, 2010-2018.
	 * 
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 2 Aug 2018
	 * @version 0.1
	 * 
	 * Toggles visibility of the scatterplot frame
	 */
	private class TeaspoonCustomGUItoggleScatterplotViewListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			System.out.println("Action Event: toggle scatterplot view.");
			toggleRegressionPlotVisible();
		}		
	}

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
	 * 
	 * Toggles visibility of the site-freq histogram frame
	 */
	private class TeaspoonCustomGUItoggleHistogramViewListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			System.out.println("Action Event: toggle histogram view.");
			toggleSiteFreqHistogramVisible();
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
			appModel.inferDates();
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
	private class TeaspoonCustomGUIrunAnalysisListener implements ActionListener, PropertyChangeListener{
		MainBhattAdaptationAnalysisTask task;
		JLabel taskLabel;
		JProgressBar taskBar;
		public final String completeText = "Done";
		public final int numberOfsiteFreqBins = 100;

		/**
		 * Invoked when task's progress property changes.
		 */
		public void propertyChange(PropertyChangeEvent evt) {
			if ("progress" == evt.getPropertyName()) {
				int progress = task.getProgress();
				String message = "Calculating site-frequency spectrum ("+progress+";%)...";
				taskLabel.setText(message);
				taskBar.setValue(progress);
			} 
		}

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			taskLabel = appView.getTaskLabel();
			taskBar = appView.getTaskBar();

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
					parameters.setRunID(runID);
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
				//runAnalysis(parameters);
				task = new MainBhattAdaptationAnalysisTask(taskLabel, taskBar, parameters, numberOfsiteFreqBins);
				task.addPropertyChangeListener((PropertyChangeListener) this);
				task.execute();
				taskLabel.setText(completeText);
				taskBar.setValue(numberOfsiteFreqBins);
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
			// guess dates
			appModel.inferDates();
			System.out.println("date 0"+(float)appModel.getValueAt(0, 3));
			System.out.println("date 2"+(float)appModel.getValueAt(2, 3));
			
			// if guessdates worked, set ancestral to lowest date
			appModel.setOldestAlignmentAsAncestral();
			
			// automatically create a mask from the data 
			if(appModel.getAlignmentLength()>0 && appModel.hasAncestralSequenceAlignmentBeenSet()){
				RateEstimationBehaviour behaviour = RateEstimationBehaviour.NEUTRAL_RATE_AGGREGATED;
				int start = 0;
				int end = appModel.getAlignmentLength()-1;
				int length = appModel.getAlignmentLength();
				TeaspoonMask newMaskTrack;
				newMaskTrack = TeaspoonMaskFactory.initialiseMask(behaviour, start, end, length);
				appModel.addMaskRow(newMaskTrack);
				appView.repaint();
				appView.showGenericDialog("A default mask with aggregated neutral ratio estaimation behaviour covering the whole alignment length has been created.");
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

	/**
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 11 May 2018
	 * @version 0.1
	 * 
	 * Remove the selected analysis masking track.
	 */
	private class TeaspoonCustomGUIcombineMaskTrackListener implements ActionListener{

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			System.out.println("Action Event: Remove a mask_mid track");
			/*
			 * 1. check exactly one row is selected
			 * 2. toggle this row's ancestral property
			 * 3. add ancestral to model data 
			 */
			try{
				if(appView.getMasksTable().getSelectedRows().length == 2){
					// toggle state
					appModel.combineMasksWithSelectedRows(appView.getMasksTable().getSelectedRows());
					appView.repaint();
				}else{
					// invalid
					throw new Exception();
				}
			}catch (Exception ex){
				System.err.println("exactly two masks must be selected");
				ex.printStackTrace();
			}
		}		
	}

	/**
	 * <b>TEASPOON:<b>
	 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
	 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
	 * University of Oxford, 2010-2018.
	 * 
	 * Shows the site-freq spectrum, calculating if needed.
	 * 
	 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
	 * @since 24 Jul 2018
	 * @version 0.1
	 */
	private class TeaspoonCustomGUIshowSpectrumListener implements ActionListener, PropertyChangeListener{
		SiteFreqPlottingTask task;
		JLabel taskLabel;
		JProgressBar taskBar;
		public final String completeText = "Done ";
		public final int numberOfsiteFreqBins = 100;

		/**
		 * Invoked when task's progress property changes.
		 */
		public void propertyChange(PropertyChangeEvent evt) {
			if ("progress" == evt.getPropertyName()) {
				int progress = task.getProgress();
				String message = "Calculating site-frequency spectrum ("+progress+";%)...";
				taskLabel.setText(message);
				taskBar.setValue(progress);
			} 
		}

		/* (non-Javadoc)
		 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
		 */
		@Override
		public void actionPerformed(ActionEvent e) {
			System.out.println("Action Event: show site-freq spectrum");

			taskLabel = appView.getTaskLabel();
			taskBar = appView.getTaskBar();
			
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
				String runID = "rapid";
				//double ratio = 0.7186788;
		//		File maskFile = new File("./HCV_data/sub_053/"+runID+".rapid.mask");
				//File ancestral = new File("./HCV_data/sub_053/FP7_05301_0.fasta");
		//		File output = new File("./HCV_data/sub_053/"+runID+".out");
				// try programmatic set
		//		ArrayList<TeaspoonMask> masks = appModel.getMaskTracks();
		//		try {
		//			TeaspoonMaskFactory.writeMaskFile(maskFile, masks);
		//		} catch (IOException ex) {
		//			// TODO Auto-generated catch block
		//			ex.printStackTrace();
		//		}
				BhattAdaptationParameters parameters = appModel.getParametersSnapshot();
		//		try {
					// hardcoded here
		//			parameters.setMaskFile(maskFile);
		//			parameters.setOutputFile(output);
					parameters.setRunID(runID);
					/*
					 * From GUI model directly:
					 * parameters.setBootstrapReplicates(numBootstraps);
					 * parameters.setAncestralFile(ancestral);
					 * parameters.setInputFileList(inputList);
					 * parameters.setNeutralRate(ratio);
					 */
		//		} catch (IOException ioEx) {
					// TODO Auto-generated catch block
		//			ioEx.printStackTrace();
		//		}
				task = new SiteFreqPlottingTask(taskLabel, taskBar, parameters, numberOfsiteFreqBins);
				task.addPropertyChangeListener((PropertyChangeListener) this);
				task.execute();
				taskLabel.setText(completeText);
				taskBar.setValue(numberOfsiteFreqBins);
			}		
		}
	}

	public class MainBhattAdaptationAnalysisTask extends SwingWorker<Void, Void> {
		JLabel taskLabel;
		JProgressBar taskBar;
		File forceOpen = null;
		public final String completeText = "Done ";
		public  int numberOfsiteFreqBins = 100;
		int numBins;
		BhattAdaptationParameters parameters;
	
		/**
		 * Set an optional tasklabel/bar
		 * @param label
		 * @param bar
		 * @param numBins 
		 * @param parametersBhattAdaptationParameters 
		 */
		public MainBhattAdaptationAnalysisTask(JLabel label, JProgressBar bar, BhattAdaptationParameters parametersBhattAdaptationParameters, int bins){
			taskLabel = label;
			taskBar = bar;
			parameters = parametersBhattAdaptationParameters;
			numBins = bins;
			numberOfsiteFreqBins = bins;
		}
	
		/**
		 * No-arg constructor - risky since taskLabel and taskBar will not be instantiated.
		 */
		public MainBhattAdaptationAnalysisTask() {
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
	//		for(int alignmentFile:files){
				// Attempt to check we have a valid filename.
	//			filesTried++;
				runAnalysis(parameters, this);
	//			progress = Math.round(((float)filesTried / (float)totalFiles)*100f);
				progress = 50;
				String message = "Adaptation analysis ("+numBins+"; "+progress+"%)...";
				taskLabel.setText(message);
				taskBar.setValue(progress);
				setProgress(Math.min(progress, numberOfsiteFreqBins));
				System.out.println("Adaptation analysis ("+numBins+"; "+progress+"%)...");
	//		}
			return null;
		}
	
		/*
		 * Executed in event dispatching thread
		 */
		@Override
		public void done() {
			taskLabel.setText(completeText);
			taskBar.setValue(numberOfsiteFreqBins);
		}

		/**
		 * @param progress
		 */
		public void updateProgress(int progress) {
			String message = "Estimating neutral ratio ("+numBins+"; "+progress+"%)...";
			taskLabel.setText(message);
			taskBar.setValue(progress);		
		}
		
		/**
		 * @param progress
		 */
		public void incrementCountsProgress(int progress) {
			String message = "Estimating non-neutral substitutions ("+numBins+"; "+progress+"%)...";
			taskLabel.setText(message);
			taskBar.setValue(progress);			
		}
		
		/**
		 * @param progress
		 */
		public void incrementCountsProgress(float progress) {
			int roundedProgress = Math.round(progress);
			String message = "Estimating non-neutral substitutions ("+numBins+"; "+roundedProgress+"%)...";
			taskLabel.setText(message);
			taskBar.setValue(roundedProgress);			
		}
		
	}

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
	 * Show the application's About message.
	 */
	class AboutMenuListener implements ActionListener{
		@Override
		public void actionPerformed(ActionEvent arg0) {
			appView.aboutFrame.setVisible(true);
		}
	}

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
	 * Opens a URL in response to user action. Follows pattern example at {@linkplain http://stackoverflow.com/questions/10967451/open-a-link-in-browser-with-java-button}
	 */
	class OpenURLListener implements ActionListener{
		String URL;
		
		/**
		 * No-arg constructor. Points to default project URL https://github.com/jnarag/teaspoon/
		 */
		public OpenURLListener(){
			URL = "https://github.com/jnarag/teaspoon";
		}
	
		/**
		 * String arg constructor. Points to URLtoOpen URL.
		 * @param URLtoOpen - url to point to.
		 */
		public OpenURLListener(String URLtoOpen){
			URL = URLtoOpen;
		}
	
		/**
		 * Opens the URL using the default desktop environment browser.
		 */
		@Override
		public void actionPerformed(ActionEvent arg0) {
			openWebpage(URL);
		}
	
		void openWebpage(URI uri) {
			Desktop desktop = Desktop.isDesktopSupported() ? Desktop.getDesktop() : null;
			if (desktop != null && desktop.isSupported(Desktop.Action.BROWSE)) {
				try {
					desktop.browse(uri);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
	
		void openWebpage(String url) {
			try {
				java.net.URL formedURL = new java.net.URL(url);
				openWebpage(formedURL.toURI());
			} catch (URISyntaxException e) {
				e.printStackTrace();
			} catch (MalformedURLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}	
	}

	public class SiteFreqPlottingTask extends SwingWorker<Void, Void> {
		JLabel taskLabel;
		JProgressBar taskBar;
		File forceOpen = null;
		public final String completeText = "Done ";
		public  int numberOfsiteFreqBins = 100;
		int numBins;
		BhattAdaptationParameters parameters;
	
		/**
		 * Set an optional tasklabel/bar
		 * @param label
		 * @param bar
		 * @param numBins 
		 * @param parametersBhattAdaptationParameters 
		 */
		public SiteFreqPlottingTask(JLabel label, JProgressBar bar, BhattAdaptationParameters parametersBhattAdaptationParameters, int bins){
			taskLabel = label;
			taskBar = bar;
			parameters = parametersBhattAdaptationParameters;
			numBins = bins;
			numberOfsiteFreqBins = bins;
		}
	
		/**
		 * No-arg constructor - risky since taskLabel and taskBar will not be instantiated.
		 */
		public SiteFreqPlottingTask() {
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
	//		for(int alignmentFile:files){
				// Attempt to check we have a valid filename.
	//			filesTried++;
				runFastSiteFreqAnalysis(parameters, numberOfsiteFreqBins, this);
	//			progress = Math.round(((float)filesTried / (float)totalFiles)*100f);
				progress = 50;
				String message = "Counting bins ("+numBins+"; "+progress+"%)...";
				taskLabel.setText(message);
				taskBar.setValue(progress);
				setProgress(Math.min(progress, numberOfsiteFreqBins));
				System.out.println("Counting bins ("+numBins+"; "+progress+"%)...");
	//		}
			return null;
		}
	
		/*
		 * Executed in event dispatching thread
		 */
		@Override
		public void done() {
			taskLabel.setText(completeText);
			taskBar.setValue(numberOfsiteFreqBins);
			toggleSiteFreqHistogramVisible();
		}
	
		/**
		 * @param whichBin
		 */
		public void updateProgress(int whichBin) {
			String message = "Counting bins ("+numBins+"; "+whichBin+"%)...";
			taskLabel.setText(message);
			taskBar.setValue(whichBin);
			
		}
	}
}