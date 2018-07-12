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
 * Mask file parameter for BhattAdaptationAnalysis
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattMaskInputFileParameter extends AbstractBhattParameter {

	/**
	 * Explicit constructor with a java.io.File input type.
	 * @param inputFile - sequence mask_mid
	 */
	public BhattMaskInputFileParameter(File inputFile){
		this.paramType = BhattParameterType.INPUT_FILE_MASK;
		this.paramValue = inputFile;
	}
	
	/**
	 * Returns the mask_mid file
	 * @return the mask_mid file
	 */
	public File getValue(){
		return (File) paramValue;
	}
}
