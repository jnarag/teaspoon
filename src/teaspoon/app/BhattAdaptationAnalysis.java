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
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattAdaptationAnalysis {

	private BhattAdaptationParameters analysisParameters;
	private BhattAdaptationFullSiteMatrix mainAlignment, ancestralAlignment;

	private final double[][] bins = {
			{0.0, 0.15, 0.75},
			{0.15, 0.75, 1.0}
	};
	private final boolean[] Nvec = {false,true,false};
	private final double[] prior = {1.0,1.0,1.0,1.0};

	/**
	 * 
	 */
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
		 * 3. for each mask, run submatrix bhatt counts
		 * 4. populate and return results
		 */
		
        // load main
		mainAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(analysisParameters.getInputFile()).readFASTA());
        
		// load ancestral
		ancestralAlignment = new BhattAdaptationFullSiteMatrix(new MainAlignmentParser(analysisParameters.getAncestralFile()).readFASTA());
        
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

		System.out.println( 
				siteCounter.getSilentSubstitutionsCountArray()[(int) 0] 		+ "," + 
						siteCounter.getReplacementSubstitutionsCountArray()[(int) 0] + "," + 
						siteCounter.getReplacementToSilentRatesRatio()[(int) 0] 		+ "," + 
						siteCounter.getNonNeutralSubstitutions()[(int) 0]
						);

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
		 * 3. for each mask, run submatrix bhatt counts
		 * 4. populate and return results
		 * (5) (should we auto-update NR..?)
		 */
		return null;
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
		 * 3. for each mask, run submatrix bhatt counts
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
		 * 3. for each mask, run submatrix bhatt counts
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
