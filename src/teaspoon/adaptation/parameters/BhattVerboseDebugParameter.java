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
public class BhattVerboseDebugParameter extends AbstractBhattParameter {

	/**
	 * Explicit constructor with a boolean input type.
	 * @param boolean - toggles do verbose debug behaviour
	 */
	public BhattVerboseDebugParameter(boolean doDebug ){
		this.paramType = BhattParameterType.INPUT_FILE_MAIN;
		this.paramValue = doDebug;
	}
	
	/**
	 * Returns the state of the debug flag.
	 * @return the debug flag
	 */
	public boolean getValue(){
		return (boolean) paramValue;
	}

}
