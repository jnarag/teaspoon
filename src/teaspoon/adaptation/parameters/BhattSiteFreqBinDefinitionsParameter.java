/**
 * 
 */
package teaspoon.adaptation.parameters;


/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * Site frequency bin intervals parameter for BhattAdaptationAnalysis.
 * This is specified equivalently to:
 * 
	// Bins which define the frequency ranges
	private final double[][] bins = {
			{0.0, 0.15, 0.75},
			{0.15, 0.75, 1.0}
	};

 * ... Meaning also by implication that the central bin is
 * specified as the 'neutral' one for purposes of setting the
 * rates ratio, e.g. 
 * 
	// Which bin(s) to use to calculate a neutral ration when estimating
	private final boolean[] Nvec = {false,true,false};
 * 
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattSiteFreqBinDefinitionsParameter extends AbstractBhattParameter {

	/**
	 * Explicit constructor with a double[][] input type.
	 * @param double[][] - low and high intervals for bins
	 */
	public BhattSiteFreqBinDefinitionsParameter(double[][] bins){
		this.paramType = BhattParameterType.SITE_FREQUENCY_BINS;
		this.paramValue = bins;
	}
	
	/**
	 * Returns the double[][] - low and high intervals for bins.
	 * @return the double[][] - low and high intervals for bins
	 */
	public double[][] getValue(){
		return (double[][]) paramValue;
	}
}
