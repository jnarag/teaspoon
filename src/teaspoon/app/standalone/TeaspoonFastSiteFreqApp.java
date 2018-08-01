/**
 * 
 */
package teaspoon.app.standalone;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import teaspoon.adaptation.parameters.AbstractBhattParameter;
import teaspoon.adaptation.parameters.BhattMainInputFileParameter;
import teaspoon.app.BhattAdaptationAnalysis;
import teaspoon.app.TeaspoonBootstrap;
import teaspoon.app.TeaspoonMask;
import teaspoon.app.GUI.controllers.TeaspoonController.SiteFreqPlottingTask;
import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.BhattAdaptationResults;
import teaspoon.app.utils.MainAlignmentParser;
import teaspoon.app.utils.NullNeutralRatioException;
import teaspoon.app.utils.TeaspoonMethods;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 29 Jul 2018
 * @version 0.1
 * 
 * FastSiteFreqApp is a simple application to produce an approximate site-
 * frequency distribution, using the BhattAdaptationAnalysis.runFastSiteFreq() 
 * method, which in turn calls BhattMethod.inferExplicitSiteFreqHistogram().
 * It will return approximate counts of the number of substitutions in each 
 * frequency bin in the alignment (from [0,1]) discretized into <i>N</i> bins 
 * of equal width. 
 * 
 * @see GUIAnalysis
 * @see BhattMethod
 */
public class TeaspoonFastSiteFreqApp {

	static File debugAncestralFile = new File("./H7N9_flu/H7_1stWave.fasta");
	static File debugMainFile = new File("./H7N9_flu/PRD_waves_year_W2.fasta");	
	private HashMap<File, float[][]> resultsHash;

	/**
	 * Constructor actually called by main() to run a full analysis.
	 * Key parameter is the bin width which will be 1/N.
	 * 
	 * A SiteFreqPlottingTask is passed so the UI can be updated as 
	 * this method can be slowish with large datasets/bin counts.
	 * 
	 * @param siteFreqPlottingTask 
	 * @param parameters
	 * @param numbins - number of bins in the histogram from 0->1 (equal-sized)
	 * @throws IOException 
	 */
	public TeaspoonFastSiteFreqApp(BhattAdaptationParameters analysisMasterParameters, int numBins, SiteFreqPlottingTask siteFreqPlottingTask) throws IOException {
		// the ancestral file, output file, and masking file
		File ancestralFile = analysisMasterParameters.getAncestralFile();//, outputFile = analysisMasterParameters.getOutputFile(), maskFile = analysisMasterParameters.getMaskFile();

		// the list of main files
		File[] mainFiles = analysisMasterParameters.getInputFileList();

		// the ancestral alignment and combined main alignments
		BhattAdaptationFullSiteMatrix ancestralAlignment, combinedMainAlignment = null;

		// output hash
		resultsHash = new HashMap<File,float[][]>();

		// Read in ancestral alignment and open output file
		ancestralAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(ancestralFile).readFASTA());


		// Walk through main files to get alignments and concatenate them (vertically)
		for(File mainFile:mainFiles){
			BhattAdaptationFullSiteMatrix mainAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(mainFile).readFASTA());
			if(combinedMainAlignment == null){
				combinedMainAlignment = mainAlignment;
			}else{
				combinedMainAlignment = combinedMainAlignment.appendTaxa(mainAlignment);
			}
		}

		// reset parameters for the combined alignment
		analysisMasterParameters.setInputFullSiteMatrix(combinedMainAlignment);
		analysisMasterParameters.setAncestralFullSiteMatrix(ancestralAlignment);
		
		// get the empirical site-freq for this ancestral alignment, combined alignment, and bin count
		float[][] empirical = new BhattAdaptationAnalysis(analysisMasterParameters).runFastSiteFreq(numBins, siteFreqPlottingTask);

		// add to hash
		resultsHash.put(new File("fast"),empirical);

		// write out
		System.out.println("fast"+ " processed");
	}

	/**
	 * No-arg constructor is deprecated.
	 */
	@Deprecated
	public TeaspoonFastSiteFreqApp(){}

	/**
	 * Provides a basic debug capability using the static input/ancestral files, relatively specified:
	 * 	static File debugAncestralFile = new File("./H7N9_flu/H7_1stWave.fasta");
	 *	static File debugMainFile = new File("./H7N9_flu/PRD_waves_year_W2.fasta");	
	 * @param basicDebug (ignored as long as one is passed to call this overloaded constructor.)
	 */
	public TeaspoonFastSiteFreqApp(boolean basicDebug) {
		/*
		 * Example workflow:
		 * 
		 * 1. Create a parameter list
		 * 2. Fill with parameters
		 * 3. Call one of the .runWith..() methods, returning the results
		 * 4. print the results out
		 * 
		 * (optional) Add nre parameters to list, change values
		 * (optional) re-run on the same matrix object
		 */
		BhattAdaptationParameters parameters = new BhattAdaptationParameters();
		// now populate the list
		try {
			parameters.setAncestralFile(debugAncestralFile);
			parameters.setInputFile(debugMainFile);
			parameters.setNeutralRate(0.7186788);
			double[][] customBins = {
					{0.0, 0.15, 0.75},
					{0.15, 0.75, 1.0}
			};
			parameters.setCustomBinSettings(customBins);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// run the thing
		BhattAdaptationResults results = new BhattAdaptationAnalysis(parameters).runWithFixedNR();
		// print results
		results.printToFile(new File("out"));
	}

	/**
	 * @param args
	 * @throws FileNotFoundException, RuntimeException 
	 */
	public static void main(String[] args) throws FileNotFoundException, RuntimeException {
		/*
		 * The real implementation (v0.1.3 and higher)
		 * 
		 * Required arguments:
		 * - Ancestral file (ancestral sequence alignment)
		 * - List of main alignment files 
		 * - Mask file (list of masks and corresponding rate estimation behaviour for each)
		 * - Output file to write to
		 * 
		 * Optional arguments
		 * - Debug flag (true==verbose)
		 * - Number of bootstrap replicates (default: 0)
		 * - neutral ratio
		 * 
		 * Usage:
		 * <mask_mid> <ancestral> <output> <mainfiles>
		 * <mask_mid> <ancestral> <output> <mainfiles> <ratio>
		 * <mask_mid> <ancestral> <output> <mainfiles> <bootstraps>
		 * <mask_mid> <ancestral> <output> <mainfiles> <bootstraps> <ratio>
		 * <debug=true><mask_mid> <ancestral> <output> <mainfiles>
		 * <debug=true><mask_mid> <ancestral> <output> <mainfiles> <ratio>
		 * <debug=true><mask_mid> <ancestral> <output> <mainfiles> <bootstraps>
		 * <debug=true><mask_mid> <ancestral> <output> <mainfiles> <bootstraps> <ratio>
		 * 
		 * 
		 * TODO implement version
		 * TODO implement help
		 */

		/* Now build a parameter set then pass to new runnable instance */
		BhattAdaptationParameters parameters = new BhattAdaptationParameters();

		// first check for debug flag
		if(args[0].equals("true")){
			// parse with debug on
			parameters.setDebugFlag(true);

			// now populate the list.
			// first 4 args are: <mask_mid> <ancestral> <output> <mainfiles>
			File maskFile = null, ancestralFile = null, outputFile = null;
			File[] mainAlignments = null;

			switch(args.length){
			case 5:{
				/* file with the sequence mask_mid */
				maskFile = new File(args[1]);
				if(!maskFile.canRead()){
					throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+maskFile.getAbsolutePath());
				}
				/* file with the sequence alignment for ancestral / outgroup */
				ancestralFile = new File(args[2]);
				if(!ancestralFile.canRead()){
					throw new FileNotFoundException("Could not parse or read ancestral alignment file "+ancestralFile.getAbsolutePath());
				}	
				/* file for output */
				outputFile = new File(args[3]);

				/* 1 or more multiple sequence alignments for the main/focal groups */
				String[] mainfilesList = args[4].split(",");
				mainAlignments = new File[mainfilesList.length];
				for(int mainfile=0;mainfile<mainfilesList.length;mainfile++){
					mainAlignments[mainfile] = new File(mainfilesList[mainfile]);
					if(!mainAlignments[mainfile].canRead()){
						throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+mainAlignments[mainfile].getAbsolutePath());
					}
				}			
				break;
			}
			case 6:{
				// args are <mask_mid> <ancestral> <output> <mainfiles> <ratio|num_bootstraps>; need to work out which

				// first parse files

				/* file with the sequence mask_mid */
				maskFile = new File(args[1]);
				if(!maskFile.canRead()){
					throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+maskFile.getAbsolutePath());
				}
				/* file with the sequence alignment for ancestral / outgroup */
				ancestralFile = new File(args[2]);
				if(!ancestralFile.canRead()){
					throw new FileNotFoundException("Could not parse or read ancestral alignment file "+ancestralFile.getAbsolutePath());
				}	
				/* file for output */
				outputFile = new File(args[3]);

				/* 1 or more multiple sequence alignments for the main/focal groups */
				String[] mainfilesList = args[4].split(",");
				mainAlignments = new File[mainfilesList.length];
				for(int mainfile=0;mainfile<mainfilesList.length;mainfile++){
					mainAlignments[mainfile] = new File(mainfilesList[mainfile]);
					if(!mainAlignments[mainfile].canRead()){
						throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+mainAlignments[mainfile].getAbsolutePath());
					}
				}			

				// see if sixth arg is int or double
				try {
					// probably bootstrap
					int bootstraps = Integer.parseInt(args[4]);
					parameters.setBootstrapReplicates(bootstraps);
				} catch (NumberFormatException e) {
					// probably a rate
					parameters.setNeutralRate(Double.parseDouble(args[4]));
				}
				break;
			}
			case 7:{
				// <mask_mid> <ancestral> <output> <mainfiles> <bootstraps> <ratio>
				// first parse files

				/* file with the sequence mask_mid */
				maskFile = new File(args[1]);
				if(!maskFile.canRead()){
					throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+maskFile.getAbsolutePath());
				}
				/* file with the sequence alignment for ancestral / outgroup */
				ancestralFile = new File(args[2]);
				if(!ancestralFile.canRead()){
					throw new FileNotFoundException("Could not parse or read ancestral alignment file "+ancestralFile.getAbsolutePath());
				}	
				/* file for output */
				outputFile = new File(args[3]);

				/* 1 or more multiple sequence alignments for the main/focal groups */
				String[] mainfilesList = args[4].split(",");
				mainAlignments = new File[mainfilesList.length];
				for(int mainfile=0;mainfile<mainfilesList.length;mainfile++){
					mainAlignments[mainfile] = new File(mainfilesList[mainfile]);
					if(!mainAlignments[mainfile].canRead()){
						throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+mainAlignments[mainfile].getAbsolutePath());
					}
				}	

				// parse bootstraps and ratio, assumed to be in that order
				parameters.setBootstrapReplicates(Integer.parseInt(args[5]));
				parameters.setNeutralRate(Double.parseDouble(args[6]));
				break;
			}
			}

			if(!maskFile.canRead()){
				throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+maskFile.getAbsolutePath());
			}
			if(!ancestralFile.canRead()){
				throw new FileNotFoundException("Could not parse or read ancestral alignment file "+ancestralFile.getAbsolutePath());
			}


			try {
				parameters.setMaskFile(maskFile);
				parameters.setAncestralFile(ancestralFile);
				parameters.setOutputFile(outputFile);
				parameters.setInputFileList(mainAlignments);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}


			/* run it... */

			try {
				new TeaspoonFastSiteFreqApp(parameters, 30,null);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}else{
			// no debug required

			// now populate the list.
			// first 4 args are: <mask_mid> <ancestral> <output> <mainfiles>
			File maskFile = null, ancestralFile = null, outputFile = null;
			File[] mainAlignments = null;

			switch(args.length){
			case 4:{
				/* file with the sequence mask_mid */
				maskFile = new File(args[0]);
				if(!maskFile.canRead()){
					throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+maskFile.getAbsolutePath());
				}
				/* file with the sequence alignment for ancestral / outgroup */
				ancestralFile = new File(args[1]);
				if(!ancestralFile.canRead()){
					throw new FileNotFoundException("Could not parse or read ancestral alignment file "+ancestralFile.getAbsolutePath());
				}	
				/* file for output */
				outputFile = new File(args[2]);

				/* 1 or more multiple sequence alignments for the main/focal groups */
				String[] mainfilesList = args[3].split(",");
				mainAlignments = new File[mainfilesList.length];
				for(int mainfile=0;mainfile<mainfilesList.length;mainfile++){
					mainAlignments[mainfile] = new File(mainfilesList[mainfile]);
					if(!mainAlignments[mainfile].canRead()){
						throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+mainAlignments[mainfile].getAbsolutePath());
					}
				}			
				break;
			}
			case 5:{
				// args are <mask_mid> <ancestral> <output> <mainfiles> <ratio|num_bootstraps>; need to work out which

				// first parse files

				/* file with the sequence mask_mid */
				maskFile = new File(args[0]);
				if(!maskFile.canRead()){
					throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+maskFile.getAbsolutePath());
				}
				/* file with the sequence alignment for ancestral / outgroup */
				ancestralFile = new File(args[1]);
				if(!ancestralFile.canRead()){
					throw new FileNotFoundException("Could not parse or read ancestral alignment file "+ancestralFile.getAbsolutePath());
				}	
				/* file for output */
				outputFile = new File(args[2]);

				/* 1 or more multiple sequence alignments for the main/focal groups */
				String[] mainfilesList = args[3].split(",");
				mainAlignments = new File[mainfilesList.length];
				for(int mainfile=0;mainfile<mainfilesList.length;mainfile++){
					mainAlignments[mainfile] = new File(mainfilesList[mainfile]);
					if(!mainAlignments[mainfile].canRead()){
						throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+mainAlignments[mainfile].getAbsolutePath());
					}
				}			

				// see if fifth arg is int or double
				try {
					// probably bootstrap
					int bootstraps = Integer.parseInt(args[4]);
					parameters.setBootstrapReplicates(bootstraps);
				} catch (NumberFormatException e) {
					// probably a rate
					parameters.setNeutralRate(Double.parseDouble(args[4]));
				}
				break;
			}
			case 6:{
				// <mask_mid> <ancestral> <output> <mainfiles> <bootstraps> <ratio>
				// first parse files

				/* file with the sequence mask_mid */
				maskFile = new File(args[0]);
				if(!maskFile.canRead()){
					throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+maskFile.getAbsolutePath());
				}
				/* file with the sequence alignment for ancestral / outgroup */
				ancestralFile = new File(args[1]);
				if(!ancestralFile.canRead()){
					throw new FileNotFoundException("Could not parse or read ancestral alignment file "+ancestralFile.getAbsolutePath());
				}	
				/* file for output */
				outputFile = new File(args[2]);

				/* 1 or more multiple sequence alignments for the main/focal groups */
				String[] mainfilesList = args[3].split(",");
				mainAlignments = new File[mainfilesList.length];
				for(int mainfile=0;mainfile<mainfilesList.length;mainfile++){
					mainAlignments[mainfile] = new File(mainfilesList[mainfile]);
					if(!mainAlignments[mainfile].canRead()){
						throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+mainAlignments[mainfile].getAbsolutePath());
					}
				}	

				// parse bootstraps and ratio, assumed to be in that order
				parameters.setBootstrapReplicates(Integer.parseInt(args[4]));
				parameters.setNeutralRate(Double.parseDouble(args[5]));
				break;
			}
			}

			if(!maskFile.canRead()){
				throw new FileNotFoundException("Could not parse or read sequence mask_mid file "+maskFile.getAbsolutePath());
			}
			if(!ancestralFile.canRead()){
				throw new FileNotFoundException("Could not parse or read ancestral alignment file "+ancestralFile.getAbsolutePath());
			}


			try {
				parameters.setMaskFile(maskFile);
				parameters.setAncestralFile(ancestralFile);
				parameters.setOutputFile(outputFile);
				parameters.setInputFileList(mainAlignments);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}


			/* run it... */

			try {
				new TeaspoonFastSiteFreqApp(parameters, 30,null);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	/**
	 * Gets the results of this site-frequency analysis. Returns float[][] 
	 * (wrapped in a hash) comprising (N*2) binned site-frequency estimates, 
	 * such that:
	 *  returnArray[N][0] - right-hand bound for Nth bin (bin sequence ranges [0,1])
	 *  returnArray[N][1] - estimated number of replacement substitutions in this bin range.
	 *  
	 * Note that many bins will contain NaNs.
	 * 
	 * @return - float[N][2] return array of bin postions and replacement site frequencies.
	 * @throws NullPointerException
	 */
	public HashMap<File, float[][]> getResults() throws NullPointerException{
		if(this.resultsHash != null){
			return this.resultsHash;
		}else{
			throw new NullPointerException("Results not available or null!");
		}
	}
}
