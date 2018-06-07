/**
 * 
 */
package teaspoon.app;

import teaspoon.app.utils.BhattAdaptaionFullSiteMatrix;
import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.BhattAdaptationResults;


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
		BhattAdaptationFullSiteMatrix siteData = new BhattAdaptationFullSiteMatrix();
		siteData.loadAlignmentFile(analysisParameters.getInputFile()); //possibly remove no-arg FSM constructor and do this in init
		
		return null;
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
