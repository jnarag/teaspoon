/**
 * 
 */
package teaspoon.app.standalone;

import java.io.File;
import java.io.FileNotFoundException;
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
	 */
	public TeaspoonCommandLineApp(BhattAdaptationParameters parameters) {
		// the ancestral file and masking file
		File ancestralFile = null, maskFile = null;
		// the list of main files
		File[] mainFiles = null;
		// the mask positions
		TeaspoonMask[] masks;
		// the bootstrap positions
		int bootstrapReplicates = 0;
		TeaspoonBootstrap[] bootstraps;
		// linked alignments and filenames for main alignments
		HashMap<File,BhattAdaptationFullSiteMatrix> mainAlignments = null;
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
		masks = TeaspoonMaskFactory(maskFile);
		
		for(TeaspoonMask mask:masks){

			// [4] Execute rate estimation behaviour to get neutral rate if needed
			
			switch(mask.estimationBehaviour){
			case NEUTRAL_RATE_AGGREGATED:{
				double estimatedRate = 0;
				// aggregate over all main alignments e.g. use combinedAlignment
				// set the neutral rate based on this
				mask.setNeutralRate(estimatedRate);
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
				mask.setNeutralRate(estimatedRate);
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
				BhattAdaptationFullSiteMatrix ancestralPartition = ancestralAlignment.maskBy(mask);
				BhattAdaptationFullSiteMatrix mainPartition = mainAlignments.get(mainFile).maskBy(mask);
				
				// BhattAdaptationResults empirical = new BhattAdaptationAnalysis(ancestralPartition,mainPartition,estimatedRate);
				
				if(bootstrapReplicates > 0){
					// [7] Walk through bootstraps
					BhattAdaptationResults[] bootstrappedResults = new BhattAdaptationResults[bootstrapReplicates];
					int bootstrapCounter = 0;
					int seed = 0;
					bootstraps = TeaspoonBootstrapFactory.generate(mainPartition,bootstrapReplicates,seed);
					
					for(TeaspoonBootstrap bootstrap:bootstraps){
						// [8] Get bootstrap estimate
						BhattAdaptationFullSiteMatrix bootstrappedPartitionAncestral = ancestralPartition.obtainBoostrap(bootstrap);
						BhattAdaptationFullSiteMatrix bootstrappedParitionMain = mainPartition.obtainBoostrap(bootstrap);

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
	 * @param maskFile
	 * @return
	 */
	private TeaspoonMask[] TeaspoonMaskFactory(File maskFile) {
		// TODO Auto-generated method stub
		return null;
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
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// FIXME prints out all args then halts for now.
		for(String CLIargument: args){
			System.out.println("\t"+CLIargument);
		}
		new TeaspoonCommandLineApp();
	}

}
