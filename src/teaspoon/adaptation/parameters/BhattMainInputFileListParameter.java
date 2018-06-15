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
 * Main file list parameter for BhattAdaptationAnalysis
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattMainInputFileListParameter extends AbstractBhattParameter {

	/**
	 * Explicit constructor with a java.io.File input type.
	 * @param inputFile - list of main files
	 */
	public BhattMainInputFileListParameter(File[] inputFileList){
		this.paramType = BhattParameterType.INPUT_FILE_LIST_MAIN;
		this.paramValue = inputFileList;
	}
	
	/**
	 * Returns the main file list
	 * @return the main file list
	 */
	public File[] getValue(){
		return (File[]) paramValue;
	}
}
