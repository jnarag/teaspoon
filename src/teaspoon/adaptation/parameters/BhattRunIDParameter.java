/**
 * 
 */
package teaspoon.adaptation.parameters;

import java.io.File;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * File parameter for BhattAdaptationAnalysis
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattRunIDParameter extends AbstractBhattParameter {

	/**
	 * Explicit constructor with a java.io.File input type.
	 * @param inputFile
	 */
	public BhattRunIDParameter(String runid){
		this.paramType = BhattParameterType.RUN_ID;
		this.paramValue = runid;
	}
	
	/**
	 * Returns the input file
	 * @return the input file
	 */
	public String getValue(){
		return (String) paramValue;
	}
}
