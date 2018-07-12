/**
 * 
 */
package teaspoon.app;

import teaspoon.adaptation.BhattMethod;
import teaspoon.app.utils.AncestralAlignmentParser;
import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.BhattAdaptationResults;
import teaspoon.app.utils.MainAlignmentParser;
import teaspoon.app.utils.NullNeutralRatioException;


/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
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
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattAdaptationAnalysis {

	private BhattAdaptationParameters analysisParameters;
	private BhattAdaptationFullSiteMatrix mainAlignment, ancestralAlignment;

	// Bins which define the frequency ranges. defaults set here but see BhattSiteFreqBinsParameter
	private double[][] bins = {
			{0.0, 0.15, 0.75},
			{0.15, 0.75, 1.0}
	};
	
	// Which bin(s) to use to calculate a neutral ration when estimating
	private final boolean[] Nvec = {false,true,false};
	
	// Priors on sites
	private final double[] prior = {1.0,1.0,1.0,1.0};

	/**
	 * No-arg constructor is deprecated.
	 */
	@Deprecated
	public BhattAdaptationAnalysis() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param parameterList
	 */
	public BhattAdaptationAnalysis(BhattAdaptationParameters parameters) {
		analysisParameters = parameters;
		// TODO Auto-generated constructor stub
	}

	/**
	 * Runs the Bhatt/Raghwani analysis with fixed neutral ratio.
	 * @return - a results object 
	 */
	public BhattAdaptationResults runWithFixedNR() {
		// TODO Auto-generated method stub
		/*
		 * Workflow:
		 * 
		 * 0. check NR exists and is non-negative. die if not
		 * 1. read in the input file, create a master matrix for anc and main
		 * 2. create a consensus. clean both
		 * 3. for each mask_mid, run submatrix bhatt counts
		 * 4. populate and return results
		 */
		
		// see if bins exist
		if(analysisParameters.hasCustomBinSettings()){
			this.bins = analysisParameters.getCustomBinSettings();
		}
			
        // load main - check to see if the main FSM exists already
		if(analysisParameters.hasFullSiteMatrixMain()){
			mainAlignment = analysisParameters.getFullSiteMatrixMain();
		}else{
			mainAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(analysisParameters.getInputFile()).readFASTA());
		}
		
		// load ancestral
		if(analysisParameters.hasFullSiteMatrixMain()){
			ancestralAlignment = analysisParameters.getFullSiteMatrixAncestral();
		}else{
			ancestralAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(analysisParameters.getAncestralFile()).readFASTA());
		}
        
		// assume cleaning occurs somewhere
		
		// count via the BhattMethod and get results; check for debug flag though
		BhattMethod siteCounter;
		if(analysisParameters.getDoDebugFlag()){
			// do verbose debug
			siteCounter = new BhattMethod(mainAlignment.getSiteMatrix(), ancestralAlignment.deriveConsensus(),true);
		}else{
			// no debug needed
			siteCounter = new BhattMethod(mainAlignment.getSiteMatrix(), ancestralAlignment.deriveConsensus());
		}
		
		// run counts
        try {
        	// bins, prior, Nvec all hardcoded for now.
			siteCounter.inferCountsFixedNR(bins, prior, true, Nvec, analysisParameters.getNeutralRate());
		} catch (NullNeutralRatioException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        // get results, read into a new BhattAdaptationResults object and return it
        siteCounter.equals(null);

        if(analysisParameters.getDoDebugFlag()){
        	System.out.println( 
    				siteCounter.getSilentSubstitutionsCountArray()[(int) 0] 		+ "," + 
    						siteCounter.getReplacementSubstitutionsCountArray()[(int) 0] + "," + 
    						siteCounter.getReplacementToSilentRatesRatio()[(int) 0] 		+ "," + 
    						siteCounter.getNonNeutralSubstitutions()[(int) 0]
    						);
        }

        return new BhattAdaptationResults(siteCounter,analysisParameters);
	}

	/**
	 * Runs the Bhatt/Raghwani analysis with no neutral ratio which will be estimated
	 * @return - a results object which should include an estimated neutral rate
	 */
	public BhattAdaptationResults runWithEstimatedNR() {
		// TODO Auto-generated method stub
		/*
		 * Workflow:
		 * 
		 * 1. read in the input file, create a master matrix for anc and main
		 * 2. create a consensus. clean both
		 * 3. for each mask_mid, run submatrix bhatt counts
		 * 4. populate and return results
		 * (5) (should we auto-update NR..?)
		 */

		// see if bins exist
		if(analysisParameters.hasCustomBinSettings()){
			this.bins = analysisParameters.getCustomBinSettings();
		}
        // load main - check to see if the main FSM exists already
		if(analysisParameters.hasFullSiteMatrixMain()){
			mainAlignment = analysisParameters.getFullSiteMatrixMain();
		}else{
			mainAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(analysisParameters.getInputFile()).readFASTA());
		}
		
		// load ancestral - check to see if the ancestral FSM exists alreadt
		if(analysisParameters.hasFullSiteMatrixAncestral()){
			ancestralAlignment = analysisParameters.getFullSiteMatrixAncestral();
		}else{
			ancestralAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(analysisParameters.getAncestralFile()).readFASTA());
		}
        
		// assume cleaning occurs somewhere
		
		// count via the BhattMethod and get results; check for debug flag though
		BhattMethod siteCounter;
		if(analysisParameters.getDoDebugFlag()){
			// do verbose debug
			siteCounter = new BhattMethod(mainAlignment.getSiteMatrix(), ancestralAlignment.deriveConsensus(),true);
		}else{
			// no debug needed
			siteCounter = new BhattMethod(mainAlignment.getSiteMatrix(), ancestralAlignment.deriveConsensus());
		}
		
		// run counts
		// bins, prior, Nvec all hardcoded for now.
		siteCounter.inferCountsEstimatedNR(bins, prior, true, Nvec);
        
        // get results, read into a new BhattAdaptationResults object and return it
        siteCounter.equals(null);

        if(analysisParameters.getDoDebugFlag()){
        	System.out.println( 
    			siteCounter.getSilentSubstitutionsCountArray()[(int) 0] 		+ "," + 
    			siteCounter.getReplacementSubstitutionsCountArray()[(int) 0] 	+ "," + 
    			siteCounter.getReplacementToSilentRatesRatio()[(int) 0] 		+ "," + 
    			siteCounter.getNonNeutralSubstitutions()[(int) 0]				+ "," +
    			siteCounter.getNeutralRatio() 								 
    		);
        }

        return new BhattAdaptationResults(siteCounter,analysisParameters);
	}

	/**
	 * Runs the Bhatt/Raghwani analysis with fixed neutral ratio.
	 * Bootstraps results.
	 * @return - a results object 
	 */
	public BhattAdaptationResults runWithFixedNRandBootstrap() {
		// TODO Auto-generated method stub
		/*
		 * Workflow:
		 * 
		 * 0. check NR exists and is non-negative. die if not
		 * 1. read in the input file, create a master matrix for anc and main
		 * 2. create a consensus. clean both
		 * 3. for each mask_mid, run submatrix bhatt counts
		 * 4. populate and return results
		 */
		return null;
	}

	/**
	 * Runs the Bhatt/Raghwani analysis with no neutral ratio which will be estimated
	 * Bootstraps results.
	 * @return - a results object which should include an estimated neutral rate
	 */
	public BhattAdaptationResults runWithEstimatedNRandBootstrap() {
		// TODO Auto-generated method stub
		/*
		 * Workflow:
		 * 
		 * 1. read in the input file, create a master matrix for anc and main
		 * 2. create a consensus. clean both
		 * 3. for each mask_mid, run submatrix bhatt counts
		 * 4. populate and return results
		 * (5) (should we auto-update NR..?)
		 */
		return null;
	}

	/**
	 * @return the analysisParameters
	 */
	public BhattAdaptationParameters getAnalysisParameters() {
		return analysisParameters;
	}

	/**
	 * @param analysisParameters the analysisParameters to set
	 */
	public void setAnalysisParameters(BhattAdaptationParameters analysisParameters) {
		this.analysisParameters = analysisParameters;
	}

}
