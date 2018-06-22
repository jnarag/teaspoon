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
public class BhattOuputFileParameter extends AbstractBhattParameter {

	/**
	 * Explicit constructor with a java.io.File input type.
	 * @param outputFile
	 */
	public BhattOuputFileParameter(File outputFile){
		this.paramType = BhattParameterType.OUTPUT_FILE;
		this.paramValue = outputFile;
	}
	
	/**
	 * Returns the output file
	 * @return the output file
	 */
	public File getValue(){
		return (File) paramValue;
	}
}
