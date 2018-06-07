/**
 * 
 */
package teaspoon.app.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;

import teaspoon.adaptation.parameters.AbstractBhattParameter;
import teaspoon.adaptation.parameters.BhattInputFileParameter;
import teaspoon.adaptation.parameters.BhattParameterType;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * Holds all the parameters for an BhattAdaptationAnalysis
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattAdaptationParameters {

	private HashMap<BhattParameterType,AbstractBhattParameter> parameters;

	/**
	 * no-arg constructor
	 */
	public BhattAdaptationParameters() {
		// TODO Auto-generated constructor stub
		parameters = new HashMap<BhattParameterType,AbstractBhattParameter>();
	}

	/**
	 * Set or update the input file
	 * @param input
	 * @return
	 * @throws FileNotFoundException
	 */
	public boolean setInputFile(File input) throws FileNotFoundException{
		if(!input.exists()){
			throw new FileNotFoundException();
		}else{
			parameters.put(BhattParameterType.INPUT_FILE, new BhattInputFileParameter(input));
			return true;
		}
	}
	
	/**
	 * Return the input file value
	 * @return - an input file
	 */
	public File getInputFile(){
		// Walk through the parameters
		return (File) parameters.get(BhattParameterType.INPUT_FILE).getParamValue();
	}
}
