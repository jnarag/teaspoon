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
import teaspoon.app.TEASPOONVersion;
import teaspoon.app.TeaspoonBootstrap;
import teaspoon.app.TeaspoonMask;
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
 * @since 5 Dec 2017
 * @version 0.1.4
 * 
 * This is the runner main-class for the command-line Teaspoon app
 * It will be the entrypoint for all command-line analyses.
 * It will be called by GUI analysis.
 * 
 * @see GUIAnalysis
 */
public class TeaspoonCommandLineApp {

	private final static TEASPOONVersion teaspoonVersion = new TEASPOONVersion();
	static File debugAncestralFile = new File("./H7N9_flu/H7_1stWave.fasta");
	static File debugMainFile = new File("./H7N9_flu/PRD_waves_year_W2.fasta");	
	private HashMap<File,BhattAdaptationResults> resultsHash;

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
	 * 	- 'mask_mid' and 'sliding window' used interchangeably. this implies windows 
	 * will take ages as all steps recalculated - we'll live with this for now
	 *  - doesn't really matter whether we have a separate FSM for each timepoint
	 *  or combine one with a really stable subslice() method, or both. the
	 *  point though is we may well want aggregate counts across all timepoints.
	 * </p>
  	 * @param parameters
	 * @throws IOException 
	 */
	public TeaspoonCommandLineApp(BhattAdaptationParameters analysisMasterParameters) throws IOException {
		// the ancestral file, output file, and masking file
		File ancestralFile = analysisMasterParameters.getAncestralFile(), outputFile = analysisMasterParameters.getOutputFile(), maskFile = analysisMasterParameters.getMaskFile();
		// the list of main files
		File[] mainFiles = analysisMasterParameters.getInputFileList();
		// the mask_mid positions
		TeaspoonMask[] masks;
		// the bootstrap positions
		int bootstrapReplicates = 0;
		try{
			bootstrapReplicates = analysisMasterParameters.getBootstrapReplicates();
		}catch (NullPointerException ex){
			bootstrapReplicates = 0;
		}
		TeaspoonBootstrap[] bootstraps;
		// linked alignments and filenames for main alignments
		HashMap<File,BhattAdaptationFullSiteMatrix> mainAlignments = new HashMap<File,BhattAdaptationFullSiteMatrix>();
		// the ancestral alignment and combined main alignments
		BhattAdaptationFullSiteMatrix ancestralAlignment, combinedMainAlignment = null;
		// output hash
		resultsHash = new HashMap<File,BhattAdaptationResults>();
		
		// the output file writer
		FileWriter writer = new FileWriter(outputFile);
		writer.write(
				"mask_mid\t"+
				"file\t"+
				"N_taxa\t"+
				"N_sites\t"+
				"method\t"+
				"neutral_ratio\t"+
				"rho_high\t"+
				"N_adaptations\t"+
				"bootstrap_information\n"
		);
		
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
		
		// [1] Read in ancestral alignment and open output file
        ancestralAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(ancestralFile).readFASTA());
        
		
		// [2] Walk through main files to get alignments
		for(File mainFile:mainFiles){
			BhattAdaptationFullSiteMatrix mainAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(mainFile).readFASTA());
			if(combinedMainAlignment == null){
				combinedMainAlignment = mainAlignment;
			}else{
				combinedMainAlignment = combinedMainAlignment.appendTaxa(mainAlignment);
			}
			mainAlignments.put(mainFile, mainAlignment);
		}
		
		// [3] Walk through masks
		masks = TeaspoonMaskFactory.parseFile(maskFile);
		
		for(TeaspoonMask mask:masks){

			// [4] Execute rate estimation behaviour to get neutral rate if needed
			
			BhattAdaptationParameters masterMaskParameters = analysisMasterParameters; // set params for this mask_mid run
			
			System.err.println("Deciding netural ratio behaviour.");
			switch(mask.estimationBehaviour){
			case NEUTRAL_RATE_AGGREGATED:{
				System.err.println("Inferring neutral ratio by aggregation.");
				BhattAdaptationParameters trainingParameters = masterMaskParameters;
				trainingParameters.setAncestralFullSiteMatrix(ancestralAlignment.subsampleByMask(mask));
				trainingParameters.setInputFullSiteMatrix(combinedMainAlignment.subsampleByMask(mask));
				
				/*
				 * Actual rate estimation
				 * 
				 * !IMPORTANT!
				 * 
				 * The rate returned may be positive double, zero, infinity, or NaN.
				 * May need to raise a custom RateEstimationError exception
				 */
				double estimatedRate = 0, rhoHigh = 0, adaptationsHigh = 0;
				BhattAdaptationResults trainingResults = new BhattAdaptationAnalysis(trainingParameters).runWithEstimatedNR();
				estimatedRate = trainingResults.getBhattSiteCounter().getNeutralRatio();
				rhoHigh = trainingResults.getBhattSiteCounter().getReplacementSubstitutionsCountArray()[2];
				adaptationsHigh = trainingResults.getBhattSiteCounter().getNonNeutralSubstitutions()[2];
				writer.write(
						mask.toString()+"\t"+
						"<aggregated>\t"+
						trainingResults.getBhattSiteCounter().integerMatrix.length+"\t"+
						trainingResults.getBhattSiteCounter().integerMatrix[0].length+"\t"+
						"A\t"+
						estimatedRate+"\t"+
						rhoHigh+"\t"+
						adaptationsHigh+"\t"+
						"NA\n"
				);
				estimatedRate = (new BhattAdaptationAnalysis(trainingParameters)).runWithEstimatedNR().getBhattSiteCounter().getNeutralRatio();
				try {
					if((estimatedRate >= 0)&&(estimatedRate != Double.POSITIVE_INFINITY)){
					}else{
						throw new NullNeutralRatioException("Infinite or NaN ratio estimated from dataset. This will be ignored from the averaging procedure.");
					}
				} catch (NullNeutralRatioException e) {
					// TODO Auto-generated catch block
					System.err.println("Infinite or NaN ratio estimated from dataset. This will be ignored from the averaging procedure.");
					if(trainingParameters.getDoDebugFlag()){
						e.printStackTrace();							
					}
				}

				// set the neutral rate based on this
				mask.setNeutralRatio(estimatedRate);
				try {
					masterMaskParameters.setNeutralRate(estimatedRate); // add NR to params for this mask_mid
					System.err.println("Inferred ratio was "+estimatedRate+". This will be used for the analysis.");
				} catch (IllegalArgumentException e) {
					// TODO Auto-generated catch block
					System.err.println("Could not infer a valid ratio (was: ["+estimatedRate+"]). Will use 0.0 for analysis.");
					estimatedRate = 0.0;
					masterMaskParameters.setNeutralRate(estimatedRate); // add NR to params for this mask_mid
					e.printStackTrace();
				}
				break;
			}
			case NEUTRAL_RATE_AVERAGED:{
				double estimatedRate = 0;
				int count = 0;
				// average over all partitions 
				System.err.println("Inferring neutral ratio by averaging.");
				Iterator<File> itr = mainAlignments.keySet().iterator();
				while(itr.hasNext()){
					File thisInput = itr.next();
					BhattAdaptationFullSiteMatrix alignment = mainAlignments.get(thisInput);
					BhattAdaptationParameters trainingParameters = masterMaskParameters;
					trainingParameters.setAncestralFullSiteMatrix(ancestralAlignment.subsampleByMask(mask));
					trainingParameters.setInputFullSiteMatrix(alignment.subsampleByMask(mask));

					/*
					 * Actual rate estimation
					 * 
					 * !IMPORTANT!
					 * 
					 * The rate returned may be positive double, zero, infinity, or NaN.
					 * May need to raise a custom RateEstimationError exception
					 */
					double alignmentRate = 0, rhoHigh = 0, adaptationsHigh = 0;
					
					BhattAdaptationResults trainingResults = new BhattAdaptationAnalysis(trainingParameters).runWithEstimatedNR();
					alignmentRate = trainingResults.getBhattSiteCounter().getNeutralRatio();
					rhoHigh = trainingResults.getBhattSiteCounter().getReplacementSubstitutionsCountArray()[2];
					adaptationsHigh = trainingResults.getBhattSiteCounter().getNonNeutralSubstitutions()[2];
					writer.write(
							mask.toString()+"\t"+
							thisInput.getName()+"\t"+
							trainingResults.getBhattSiteCounter().integerMatrix.length+"\t"+
							trainingResults.getBhattSiteCounter().integerMatrix[0].length+"\t"+
							"M\t"+
							alignmentRate+"\t"+
							rhoHigh+"\t"+
							adaptationsHigh+"\t"+
							"NA\n"
					);
					// estimate for this alignment
					try {
						if((alignmentRate >= 0)&&(alignmentRate != Double.POSITIVE_INFINITY)){
							estimatedRate += alignmentRate;
							count ++;
						}else{
							throw new NullNeutralRatioException("Infinite or NaN ratio estimated from dataset. This will be ignored from the averaging procedure.");
						}
					} catch (NullNeutralRatioException e) {
						// TODO Auto-generated catch block
						System.err.println("Infinite or NaN ratio estimated from dataset. This will be ignored from the averaging procedure.");
						if(trainingParameters.getDoDebugFlag()){
							e.printStackTrace();							
						}
					}
				}
				// average them
				estimatedRate = estimatedRate / (double) count;
				// set the neutral rate based on this - use a try-catch in case it is <0 or NaN
				System.err.println("Inferred ratio was "+estimatedRate+". This will be used for the analysis.");
				try {
					masterMaskParameters.setNeutralRate(estimatedRate); // add NR to params for this mask_mid
					System.err.println("Inferred ratio was "+estimatedRate+". This will be used for the analysis.");
				} catch (IllegalArgumentException e) {
					// TODO Auto-generated catch block
					System.err.println("Could not infer a valid ratio (was: ["+estimatedRate+"]). Will use 0.0 for analysis.");
					estimatedRate = 0.0;
					masterMaskParameters.setNeutralRate(estimatedRate); // add NR to params for this mask_mid
					e.printStackTrace();
				}
				writer.write(
						mask.toString()+"\t"+
						"<overall>\t"+
						"N_taxa_sum_of_above\t"+
						"N_sites_as_above\t"+
						"M\t"+
						estimatedRate+"\t"+
						"NA\t"+
						"NA\t"+
						"NA\n"
				);
				break;
			}
			case NEUTRAL_RATE_FIXED:
				// we can just use the neutral rate in this mask_mid
				System.err.println("Existing neutral ratio "+mask.getNeutralRatio()+" found in maskfile. This will be used.");
				masterMaskParameters.setNeutralRate(mask.getNeutralRatio()); // add NR to params for this mask_mid
				// write output
				writer.write(
						mask.toString()+"\t"+
						"<overall>\t"+
						"N_taxa\t"+
						"N_sites\t"+
						"F\t"+
						mask.getNeutralRatio()+"\t"+
						"NA\t"+
						"NA\t"+
						"NA\n"
				);
				break;
			default:
				break;		
			}

			// [5] Walk through main files
			for(File mainFile:mainFiles){

				BhattAdaptationParameters empiricalMaskedFileParameters = masterMaskParameters;
				
				// [6] Get empirical estimate
				BhattAdaptationFullSiteMatrix ancestralPartition = ancestralAlignment.subsampleByMask(mask);
				empiricalMaskedFileParameters.setAncestralFullSiteMatrix(ancestralPartition);

				BhattAdaptationFullSiteMatrix mainPartition = mainAlignments.get(mainFile).subsampleByMask(mask);
				empiricalMaskedFileParameters.setInputFullSiteMatrix(mainPartition);
				
				BhattAdaptationResults empirical = new BhattAdaptationAnalysis(empiricalMaskedFileParameters).runWithFixedNR();
				
				if(bootstrapReplicates > 0){
					// [7] Walk through bootstraps
					BhattAdaptationResults[] bootstrappedResults = new BhattAdaptationResults[bootstrapReplicates];
					int bootstrapCounter = 0;
					int seed = 0;
					bootstraps = TeaspoonBootstrapFactory.generate(mask, mainPartition,bootstrapReplicates,seed);

                    DescriptiveStatistics r_m = new DescriptiveStatistics();	// not sure why Jayna
                    DescriptiveStatistics s_m = new DescriptiveStatistics();	// uses these any more..?
                    DescriptiveStatistics bootstrappedAdaptations = new DescriptiveStatistics(); // holder for adaptations

					for(TeaspoonBootstrap bootstrap:bootstraps){
						BhattAdaptationParameters bootstrapMaskedFileParameters = masterMaskParameters;
						
						// [8] Get bootstrap estimate
						BhattAdaptationFullSiteMatrix bootstrappedPartitionAncestral = ancestralPartition.subsampleByBootstrap(bootstrap);
						bootstrapMaskedFileParameters.setAncestralFullSiteMatrix(bootstrappedPartitionAncestral);
						BhattAdaptationFullSiteMatrix bootstrappedParitionMain = mainPartition.subsampleByBootstrap(bootstrap);
						bootstrapMaskedFileParameters.setInputFullSiteMatrix(bootstrappedParitionMain);

						bootstrappedResults[bootstrapCounter] = new BhattAdaptationAnalysis(bootstrapMaskedFileParameters).runWithFixedNR();

                        // ? why Jayna ?
						if (!Double.isNaN(bootstrappedResults[bootstrapCounter].getBhattSiteCounter().getReplacementSubstitutionsCountArray()[1])) {
                            r_m.addValue(bootstrappedResults[bootstrapCounter].getBhattSiteCounter().getReplacementSubstitutionsCountArray()[1]);
                        }
                        
                        // ? why Jayna ?
                        if (!Double.isNaN(bootstrappedResults[bootstrapCounter].getBhattSiteCounter().getSilentSubstitutionsCountArray()[1])) {
                            s_m.addValue(bootstrappedResults[bootstrapCounter].getBhattSiteCounter().getSilentSubstitutionsCountArray()[1]);
                        }

                        // add the number of adapations
                        if(!Double.isNaN(bootstrappedResults[bootstrapCounter].getBhattSiteCounter().getNonNeutralSubstitutions()[2])){
                        	bootstrappedAdaptations.addValue(bootstrappedResults[bootstrapCounter].getBhattSiteCounter().getNonNeutralSubstitutions()[2]);
                        }
                        bootstrapCounter++;
					}
					
					// [9] Output if bootstraps
					// Combine bootstraps to get uncertainty for this empirical estimate and output	
					// new DescriptiveStats(empirical,bootstrappedResults)

					// add to hash
					resultsHash.put(mainFile,empirical);
					
					System.out.println(mainFile.toString()+ " processed.");
					System.out.println("empirical nonNeutral subs:");
					System.out.println(
							"\t"+empirical.getBhattSiteCounter().getNonNeutralSubstitutions()[0]+
							"\t"+empirical.getBhattSiteCounter().getNonNeutralSubstitutions()[1]+
							"\t"+empirical.getBhattSiteCounter().getNonNeutralSubstitutions()[2]
							);

					// BEGIN jayna
                    if(masterMaskParameters.getDoDebugFlag()){
    					System.out.println("boostrap results:");
    					for(int bs = 0; bs<bootstrappedResults.length;bs++){
    						System.out.println(bs+"\t"+bootstrappedResults[bs].getBhattSiteCounter().getNonNeutralSubstitutions()[2]);
    					}
                        System.out.println(">" + mainFile.getName() + ": r_m = " + r_m.getSum() + ", s_m = " + s_m.getSum() + " average_nr = " + r_m.getSum() / s_m.getSum());

                        StringBuffer sb1 = new StringBuffer();
                    	StringBuffer mid = new StringBuffer();
                    	StringBuffer low = new StringBuffer();
                    	StringBuffer high = new StringBuffer();

                    	sb1.append("patient,timepoint,range,total_sites,no_silent_sites,no_replacement_sites,Replacement/Silent Ratio,No_of_NonNeutral_changes\n");
                    	mid.append("gene,time,total_sites_mid,no_silent_sites_mid,no_replacement_sites_mid,r_mid/s_mid,no_of_adaptations\n");
                    	low.append("gene,time,total_sites_low,no_silent_sites_low,no_replacement_sites_low,r_mid/s_mid,no_of_noneutral_sites\n");
                    	high.append("gene,time,total_sites_high,no_silent_sites_high,no_replacement_sites_high,r_mid/s_mid,no_of_adaptations\n");

                    	TeaspoonMethods.record(low, mainFile.getName(), new double[]{0, 0.0, 0}, empirical.getBhattSiteCounter());
                    	TeaspoonMethods.record(mid, mainFile.getName(), new double[]{0, 0.0, 1}, empirical.getBhattSiteCounter());
                    	TeaspoonMethods.record(high, mainFile.getName(), new double[]{0, 0.0, 2}, empirical.getBhattSiteCounter());

                    	System.err.println("low:\n"+low.toString());
                    	System.err.println("mid:\n"+mid.toString());
                    	System.err.println("hi:\n"+high.toString());
                    }
		            // END jayna

					double rhoHigh = 0, adaptationsHigh = 0;
					rhoHigh = empirical.getBhattSiteCounter().getReplacementSubstitutionsCountArray()[2];
					adaptationsHigh = empirical.getBhattSiteCounter().getNonNeutralSubstitutions()[2];
					writer.write(
							mask.toString()+"\t"+
							mainFile.getName()+"\t"+
							empirical.getBhattSiteCounter().integerMatrix.length+"\t"+
							empirical.getBhattSiteCounter().integerMatrix[0].length+"\t"+
							"F\t"+
							mask.getNeutralRatio()+"\t"+
							rhoHigh+"\t"+
							adaptationsHigh+"\t"+
							"Bootstrap results (no.adaptations),"+
							"N:,"+bootstrapReplicates+","+
							"min:,"+bootstrappedAdaptations.getMin()+","+
							"25%:,"+bootstrappedAdaptations.getPercentile(25)+","+
							"mean:,"+bootstrappedAdaptations.getMean()+","+
							"median:,"+bootstrappedAdaptations.getPercentile(50)+","+
							"75%:,"+bootstrappedAdaptations.getPercentile(75)+","+
							"max:,"+bootstrappedAdaptations.getMax()+","+
							"\n"
					);

				}else{
					// [9] Output if single empirical no bootstraps

					// BEGIN jayna
					if(masterMaskParameters.getDoDebugFlag()){
						StringBuffer sb1 = new StringBuffer();
						StringBuffer mid = new StringBuffer();
						StringBuffer low = new StringBuffer();
						StringBuffer high = new StringBuffer();

						sb1.append("patient,timepoint,range,total_sites,no_silent_sites,no_replacement_sites,Replacement/Silent Ratio,No_of_NonNeutral_changes\n");
						mid.append("gene,time,total_sites_mid,no_silent_sites_mid,no_replacement_sites_mid,r_mid/s_mid,no_of_adaptations\n");
						low.append("gene,time,total_sites_low,no_silent_sites_low,no_replacement_sites_low,r_mid/s_mid,no_of_noneutral_sites\n");
						high.append("gene,time,total_sites_high,no_silent_sites_high,no_replacement_sites_high,r_mid/s_mid,no_of_adaptations\n");

						TeaspoonMethods.record(low, mainFile.getName(), new double[]{0, 0.0, 0}, empirical.getBhattSiteCounter());
						TeaspoonMethods.record(mid, mainFile.getName(), new double[]{0, 0.0, 1}, empirical.getBhattSiteCounter());
						TeaspoonMethods.record(high, mainFile.getName(), new double[]{0, 0.0, 2}, empirical.getBhattSiteCounter());

						System.err.println("low:\n"+low.toString());
						System.err.println("mid:\n"+mid.toString());
						System.err.println("hi:\n"+high.toString());
					}
		            // END jayna

					double rhoHigh = 0, adaptationsHigh = 0;
					rhoHigh = empirical.getBhattSiteCounter().getReplacementSubstitutionsCountArray()[2];
					adaptationsHigh = empirical.getBhattSiteCounter().getNonNeutralSubstitutions()[2];
					writer.write(
							mask.toString()+"\t"+
							mainFile.getName()+"\t"+
							empirical.getBhattSiteCounter().integerMatrix.length+"\t"+
							empirical.getBhattSiteCounter().integerMatrix[0].length+"\t"+
							"F\t"+
							mask.getNeutralRatio()+"\t"+
							rhoHigh+"\t"+
							adaptationsHigh+"\t"+
							"NA\n"
					);
					
					// add to hash
					resultsHash.put(mainFile,empirical);
					
					// write out
					System.out.println(mainFile.toString()+ "processed");
					System.out.println("empirical nonNeutral subs"+ empirical.getBhattSiteCounter().getNonNeutralSubstitutions()[2]);
				}
				
			}
		}
		// Lastly close the filewriter
		writer.close();
	}

	/**
	 * No-arg constructor is deprecated.
	 */
	@Deprecated
	public TeaspoonCommandLineApp(){}
	
	/**
	 * Provides a basic debug capability using the static input/ancestral files, relatively specified:
	 * 	static File debugAncestralFile = new File("./H7N9_flu/H7_1stWave.fasta");
	 *	static File debugMainFile = new File("./H7N9_flu/PRD_waves_year_W2.fasta");	
	 * @param basicDebug (ignored as long as one is passed to call this overloaded constructor.)
	 */
	public TeaspoonCommandLineApp(boolean basicDebug) {
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

		// first check for help, version, or debug flags
		if(args[0].equals("-h")){
			// print version string
			System.out.println(teaspoonVersion.getHTMLCredits());
		}else if(args[0].equals("-v")){
			// print version string
			System.out.println(teaspoonVersion.getVersion());
		}else if(args[0].equals("true")){
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
				new TeaspoonCommandLineApp(parameters);
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
				new TeaspoonCommandLineApp(parameters);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	public HashMap<File,BhattAdaptationResults> getResults() throws NullPointerException{
		if(this.resultsHash != null){
			return this.resultsHash;
		}else{
			throw new NullPointerException("Results not available or null!");
		}
	}
}
