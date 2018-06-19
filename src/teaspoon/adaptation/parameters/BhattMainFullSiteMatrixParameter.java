/**
 * 
 */
package teaspoon.adaptation.parameters;

import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 19 Jun 2018
 * @version 0.1
 */
public class BhattMainFullSiteMatrixParameter extends AbstractBhattParameter {

	/**
	 * Explicit constructor with a double input type.
	 * @param double - substitution rate
	 */
	public BhattMainFullSiteMatrixParameter(BhattAdaptationFullSiteMatrix matrix){
		this.paramType = BhattParameterType.INT_MATRIX_MAIN;
		this.paramValue = matrix;
	}
	
	/**
	 * Returns the neutral rate for the middle class.
	 * @return the substitution rate
	 */
	public BhattAdaptationFullSiteMatrix getValue(){
		return (BhattAdaptationFullSiteMatrix) paramValue;
	}
}
