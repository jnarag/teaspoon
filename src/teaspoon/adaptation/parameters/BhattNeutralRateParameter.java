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
 * Neutral rate parameter for BhattAdaptationAnalysis
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattNeutralRateParameter extends AbstractBhattParameter {

	/**
	 * Explicit constructor with a double input type.
	 * @param double - substitution rate
	 */
	public BhattNeutralRateParameter(double rate ){
		this.paramType = BhattParameterType.INPUT_FILE_MAIN;
		this.paramValue = rate;
	}
	
	/**
	 * Returns the neutral rate for the middle class.
	 * @return the substitution rate
	 */
	public double getValue(){
		return (double) paramValue;
	}
}
