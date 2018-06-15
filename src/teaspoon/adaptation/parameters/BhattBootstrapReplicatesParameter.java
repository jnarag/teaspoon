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
 * Number of bootstrap replicates for BhattAdaptationAnalysis
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattBootstrapReplicatesParameter extends AbstractBhattParameter {

	/**
	 * Explicit constructor with a double input type.
	 * @param double - substitution rate
	 */
	public BhattBootstrapReplicatesParameter(int replicates){
		this.paramType = BhattParameterType.BOOTSTRAP_REPLICATES;
		this.paramValue = replicates;
	}
	
	/**
	 * Returns the neutral rate for the middle class.
	 * @return the substitution rate
	 */
	public int getValue(){
		return (int) paramValue;
	}
}
