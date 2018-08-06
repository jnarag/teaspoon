/**
 * 
 */
package teaspoon.app.utils;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.Iterator;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import teaspoon.adaptation.BhattMethod;
import teaspoon.adaptation.parameters.BhattParameterType;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * Class holding the results data for a BhattAdaptationAnalysis
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattAdaptationResults {

	private BhattMethod siteCounter;
	private BhattAdaptationParameters parameters;
	private DescriptiveStatistics bootstrapAdaptationEstimates;
	public boolean hasBootstraps = false;
	
	/**
	 * no-arg constructor is deprecated
	 */
	@Deprecated
	public BhattAdaptationResults() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * 
	 * @param inferredSiteCounter
	 * @param inputParameters
	 */
	public BhattAdaptationResults(BhattMethod inferredSiteCounter, BhattAdaptationParameters inputParameters) {
		// TODO Auto-generated constructor stub
		this.siteCounter = inferredSiteCounter;
		this.parameters = inputParameters;
	}

	/**
	 * @param file
	 */
	public void printToFile(File file) {
		// TODO Auto-generated method stub	
	}
	
	/**
	 * Prints the results data to text
	 */
	public void printToText(){
		System.out.println("input parameters: ");
		Iterator<BhattParameterType> itr = parameters.getParameterIterator();
		while(itr.hasNext()){
			BhattParameterType param = itr.next();
			System.out.println("\t" + param + "\t" + parameters.getByKey(param));
		}
		// Check debug flag
		if(parameters.getDoDebugFlag()){
			// Print verbose output
			double[][] observedData = this.siteCounter.getObservedSiteDebugData();
			for(int siteRow=0; siteRow<observedData.length; siteRow++){
				StringBuilder output = new StringBuilder(siteRow + "\t");
				for(int siteCol=0; siteCol<observedData[0].length; siteCol++){
					output.append(observedData[siteRow][siteCol]+"\t");
				}
				System.out.println(output);
			}
			
		}
		System.out.println( "output: " +
			this.siteCounter.getSilentSubstitutionsCountArray()[(int) 0] 		+ "," + 
			this.siteCounter.getReplacementSubstitutionsCountArray()[(int) 0] 	+ "," + 
			this.siteCounter.getReplacementToSilentRatesRatio()[(int) 0] 		+ "," + 
			this.siteCounter.getNonNeutralSubstitutions()[(int) 0]
		);
	}

	/**
	 * @return the BhattMethod site counter
	 */
	public BhattMethod getBhattSiteCounter() {
		return siteCounter;
	}

	/**
	 * Add a bootstrap replicate series to this result.
	 * @param bootstrappedAdaptations
	 */
	public void setBootstrappedResults(DescriptiveStatistics bootstrappedAdaptations) {
		this.setBootstrapAdaptationEstimates(bootstrappedAdaptations);
		this.hasBootstraps = true;
	}

	/**
	 * @return the bootstrapAdaptationEstimates
	 */
	public DescriptiveStatistics getBootstrapAdaptationEstimates() {
		return bootstrapAdaptationEstimates;
	}

	/**
	 * @param bootstrapAdaptationEstimates the bootstrapAdaptationEstimates to set
	 */
	public void setBootstrapAdaptationEstimates(
			DescriptiveStatistics bootstrapAdaptationEstimates) {
		this.bootstrapAdaptationEstimates = bootstrapAdaptationEstimates;
	}

}
