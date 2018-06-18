/**
 * 
 */
package teaspoon.app.standalone;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

import teaspoon.adaptation.parameters.AbstractBhattParameter;
import teaspoon.adaptation.parameters.BhattMainInputFileParameter;
import teaspoon.app.BhattAdaptationAnalysis;
import teaspoon.app.TeaspoonBootstrap;
import teaspoon.app.TeaspoonMask;
import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.BhattAdaptationResults;
import teaspoon.app.utils.MainAlignmentParser;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 5 Dec 2017
 * @version 0.0.1
 * 
 * This is the runner main-class for the command-line Teaspoon app
 * It will be the entrypoint for all command-line analyses.
 * It will be called by GUI analysis.
 * 
 * @see GUIAnalysis
 */
public class TeaspoonCommandLineApp {

	static File debugAncestralFile = new File("./H7N9_flu/H7_1stWave.fasta");
	static File debugMainFile = new File("./H7N9_flu/PRD_waves_year_W2.fasta");	

	/**
	 * Constructor actually called by main() to run a full analysis.
	 * <pre>
	 * [pseudocode]
	 * 
	 * for f in input main alignments:
	 * 	combined = HashMap<id,FullSitesMatrix>=FullSitesMatrix(f)
	 * 
	 * for m in masks:
	 * 	BhattAdaptationResults estimated = estimateNR(combined.subslice(m))
	 * 	bootstraps = BootstrapFactory(m, replicates)
	 * 
	 *  for t in input main alignemnts:	 // (t==f, above)
	 *  	BhattAdaptationResults empirical = fixedNR(combined.subslice(m,t),empirical.getEstimatedNR())
	 *  
	 *  	for b in bootstraps:
	 *  		BhattAdaptationResults expected[b] = fixedNR(combined.subslice(m,t,b),empirical.getEstimatedNR())
	 *  
	 *  	return DescriptiveStats(empirical,expected[])
	 *  
	 *  [/pseudocode]
	 * </pre>
	 * 
	 * <p>Implementation notes<br>
	 * <b>Note</b>: 
	 * 	- 'mask' and 'sliding window' used interchangeably. this implies windows 
	 * will take ages as all steps recalculated - we'll live with this for now
	 *  - doesn't really matter whether we have a separate FSM for each timepoint
	 *  or combine one with a really stable subslice() method, or both. the
	 *  point though is we may well want aggregate counts across all timepoints.
	 * </p>
  	 * @param parameters
	 * @throws IOException 
	 */
	public TeaspoonCommandLineApp(BhattAdaptationParameters analysisMasterParameters) throws IOException {
		// the ancestral file and masking file
		File ancestralFile = analysisMasterParameters.getAncestralFile(), maskFile = analysisMasterParameters.getMaskFile();
		// the list of main files
		File[] mainFiles = analysisMasterParameters.getInputFileList();
		// the mask positions
		TeaspoonMask[] masks;
		// the bootstrap positions
		int bootstrapReplicates = analysisMasterParameters.getBootstrapReplicates();
		TeaspoonBootstrap[] bootstraps;
		// linked alignments and filenames for main alignments
		HashMap<File,BhattAdaptationFullSiteMatrix> mainAlignments = new HashMap<File,BhattAdaptationFullSiteMatrix>();
		// the ancestral alignment and combined main alignments
		BhattAdaptationFullSiteMatrix ancestralAlignment, combinedMainAlignment = null;

		// WORKFLOW:
		
		// [1] Read in ancestral alignment
		
		// [2] Walk through main files to get alignments
		
		// [3] Walk through masks
		
		// [4] Execute rate estimation behaviour to get neutral rate if needed
		
		// [5] Walk through main files
		
		// [6] Get empirical estimate
		
		// [7] Walk through bootstraps
		
		// [8] Get bootstrap estimate
		
		// [9] Combine bootstraps to get uncertainty for this empirical estimate and output
		
		
		
		// EXECUTION
		
		// [1] Read in ancestral alignment
        ancestralAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(ancestralFile).readFASTA());
		
		// [2] Walk through main files to get alignments
		for(File mainFile:mainFiles){
			BhattAdaptationFullSiteMatrix mainAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(mainFile).readFASTA());
			if(combinedMainAlignment == null){
				combinedMainAlignment = mainAlignment;
			}else{
				combinedMainAlignment.appendTaxa(mainAlignment);
			}
			mainAlignments.put(mainFile, mainAlignment);
		}
		
		// [3] Walk through masks
		masks = TeaspoonMaskFactory.parseFile(maskFile);
		
		for(TeaspoonMask mask:masks){

			// [4] Execute rate estimation behaviour to get neutral rate if needed
			
			switch(mask.estimationBehaviour){
			case NEUTRAL_RATE_AGGREGATED:{
				double estimatedRate = 0;
				// aggregate over all main alignments e.g. use combinedAlignment
				// set the neutral rate based on this
				mask.setNeutralRatio(estimatedRate);
				break;
			}
			case NEUTRAL_RATE_AVERAGED:{
				double estimatedRate = 0;
				int count = 0;
				// average over all partitions 
				Iterator<File> itr = mainAlignments.keySet().iterator();
				while(itr.hasNext()){
					BhattAdaptationFullSiteMatrix alignment = mainAlignments.get(itr.next());
					double alignmentRate = 0;
					// estimate for this alignment
					estimatedRate += alignmentRate;
					count ++;
				}
				// average them
				estimatedRate = estimatedRate / (double) count;
				// set the neutral rate based on this
				mask.setNeutralRatio(estimatedRate);
				break;
			}
			case NEUTRAL_RATE_FIXED:
				// we can just use the neutral rate in this mask
				break;
			default:
				break;		
			}

			// [5] Walk through main files
			for(File mainFile:mainFiles){

				// [6] Get empirical estimate
				BhattAdaptationFullSiteMatrix ancestralPartition = ancestralAlignment.subsampleByMask(mask);
				BhattAdaptationFullSiteMatrix mainPartition = mainAlignments.get(mainFile).subsampleByMask(mask);
				
				// BhattAdaptationResults empirical = new BhattAdaptationAnalysis(ancestralPartition,mainPartition,estimatedRate);
				
				if(bootstrapReplicates > 0){
					// [7] Walk through bootstraps
					BhattAdaptationResults[] bootstrappedResults = new BhattAdaptationResults[bootstrapReplicates];
					int bootstrapCounter = 0;
					int seed = 0;
					bootstraps = TeaspoonBootstrapFactory.generate(mask, mainPartition,bootstrapReplicates,seed);
					
					for(TeaspoonBootstrap bootstrap:bootstraps){
						// [8] Get bootstrap estimate
						BhattAdaptationFullSiteMatrix bootstrappedPartitionAncestral = ancestralPartition.subsampleByBootstrap(bootstrap);
						BhattAdaptationFullSiteMatrix bootstrappedParitionMain = mainPartition.subsampleByBootstrap(bootstrap);

						// BhattAdaptationResults[bootstrappCounter] = new BhattAdaptationAnalysis(bootstrappedPartitionAncestral,bootstrappedParitionMain,estimatedRate);
						
						bootstrapCounter++;
					}
					
					// [9] Combine bootstraps to get uncertainty for this empirical estimate and output	
					// new DescriptiveStats(empirical,bootstrappedResults)
				}else{
					// [9] Output	
				}
				
			}
		}
	}


	/**
	 * 
	 */
	public TeaspoonCommandLineApp() {
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
		 * Basic refactor test harness (v0.1.2 and lower)
		 */
		// TODO Auto-generated method stub
		// FIXME prints out all args then halts for now.
		for(String CLIargument: args){
			System.out.println("\t"+CLIargument);
		}
		new TeaspoonCommandLineApp();
		
		/*
		 * The real implementation (v0.1.3 and higher)
		 * 
		 * Needed arguments:
		 * - Ancestral file (ancestral sequence alignment)
		 * - List of main alignment files 
		 * - Mask file (list of masks and corresponding rate estimation behaviour for each)
		 * - Number of bootstrap replicates
		 */
		if(args.length >= 4){
			// Assume we have at least one main alignment, parse arguments and check as we go
			
			/* how many bootstrap replicates (0 is legal) */
			int bootstrapReplicates = (int)Integer.parseInt(args[0]);
			if(bootstrapReplicates < 0){
				throw new RuntimeException("Could not parse bootstrap replicates argument sensibly. Bootstrap replicates must be 0 or greater.");
			}
			
			/* file with the sequence alignment for ancestral / outgroup */
			File ancestralFile = new File(args[1]);
			if(!ancestralFile.canRead()){
				throw new FileNotFoundException("Could not parse or read ancestral alignment file "+ancestralFile.getAbsolutePath());
			}
			
			/* file with the sequence mask */
			File maskFile = new File(args[2]);
			if(!maskFile.canRead()){
				throw new FileNotFoundException("Could not parse or read sequence mask file "+maskFile.getAbsolutePath());
			}
			
			/* 1 or more multiple sequence alignments for the main/focal groups */
			File [] mainAlignments = new File[args.length-3];
			for(int argIndex=3;argIndex<args.length;argIndex++){
				mainAlignments[argIndex-3] = new File(args[argIndex]);
				if(!mainAlignments[argIndex-3].canRead()){
					throw new FileNotFoundException("Could not parse or read sequence mask file "+mainAlignments[argIndex-3].getAbsolutePath());
				}
			}
			
			/* Now build a parameter set then pass to new runnable instance */
			BhattAdaptationParameters parameters = new BhattAdaptationParameters();
			// now populate the list
			try {
				parameters.setBootstrapReplicates(bootstrapReplicates);
				parameters.setAncestralFile(ancestralFile);
				parameters.setMaskFile(maskFile);
				parameters.setInputFileList(mainAlignments);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			try {
				new TeaspoonCommandLineApp(parameters);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

}
