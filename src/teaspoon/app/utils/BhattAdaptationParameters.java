/**
 * 
 */
package teaspoon.app.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Iterator;

import teaspoon.adaptation.parameters.AbstractBhattParameter;
import teaspoon.adaptation.parameters.BhattBootstrapReplicatesParameter;
import teaspoon.adaptation.parameters.BhattMainInputFileListParameter;
import teaspoon.adaptation.parameters.BhattMainInputFileParameter;
import teaspoon.adaptation.parameters.BhattMaskInputFileParameter;
import teaspoon.adaptation.parameters.BhattNeutralRateParameter;
import teaspoon.adaptation.parameters.BhattParameterType;
import teaspoon.adaptation.parameters.BhattVerboseDebugParameter;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * Holds all the parameters for an BhattAdaptationAnalysis.
 * Values are stored in a HashMap with BhattParameterType as 
 * key - this <b>STRONGLY</b> implies each parameter type binds
 * to one and only one value.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattAdaptationParameters {

	private HashMap<BhattParameterType,AbstractBhattParameter> parameters;

	/**
	 * no-arg constructor. debug flag defaults to false
	 */
	public BhattAdaptationParameters() {
		// TODO Auto-generated constructor stub
		parameters = new HashMap<BhattParameterType,AbstractBhattParameter>();
		parameters.put(BhattParameterType.DO_VERBOSE_DEBUG, new BhattVerboseDebugParameter(false));
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
			parameters.put(BhattParameterType.INPUT_FILE_MAIN, new BhattMainInputFileParameter(input));
			return true;
		}
	}
	
	/**
	 * Return the input file value
	 * @return - an input file
	 */
	public File getInputFile(){
		// Walk through the parameters
		return (File) parameters.get(BhattParameterType.INPUT_FILE_MAIN).getParamValue();
	}

	/**
	 * Set or update the input file
	 * @param input
	 * @return
	 * @throws FileNotFoundException
	 */
	public boolean setInputFileList(File[] inputList) throws FileNotFoundException{
		for(File input:inputList){
			if(!input.exists()){
				throw new FileNotFoundException();
			}
		}
		parameters.put(BhattParameterType.INPUT_FILE_LIST_MAIN, new BhattMainInputFileListParameter(inputList));
		return true;
	}
	
	/**
	 * Return the input file value
	 * @return - an input file
	 */
	public File[] getInputFileList(){
		// Walk through the parameters
		return (File[]) parameters.get(BhattParameterType.INPUT_FILE_LIST_MAIN).getParamValue();
	}

	/**
	 * Set or update the mask file
	 * @param input
	 * @return
	 * @throws FileNotFoundException
	 */
	public boolean setMaskFile(File input) throws FileNotFoundException{
		if(!input.exists()){
			throw new FileNotFoundException();
		}else{
			parameters.put(BhattParameterType.INPUT_FILE_MASK, new BhattMaskInputFileParameter(input));
			return true;
		}
	}
	
	/**
	 * Return the mask file value
	 * @return - an input file
	 */
	public File getMaskFile(){
		// Walk through the parameters
		return (File) parameters.get(BhattParameterType.INPUT_FILE_MASK).getParamValue();
	}

	/**
	 * Set or update the ancestral input file
	 * @param input
	 * @return
	 * @throws FileNotFoundException
	 */
	public boolean setAncestralFile(File input) throws FileNotFoundException{
		if(!input.exists()){
			throw new FileNotFoundException("File not found: " + input.getAbsolutePath());
		}else{
			parameters.put(BhattParameterType.INPUT_FILE_ANCESTRAL, new BhattMainInputFileParameter(input));
			return true;
		}
	}
	
	/**
	 * Return the ancestral input file value
	 * @return - an input file
	 */
	public File getAncestralFile(){
		// Walk through the parameters
		return (File) parameters.get(BhattParameterType.INPUT_FILE_ANCESTRAL).getParamValue();
	}

	/**
	 * @return
	 */
	public double getNeutralRate() {
		// TODO Auto-generated method stub
		return (double) parameters.get(BhattParameterType.FIXED_NEUTRAL_RATE).getParamValue();
	}
	
	/**
	 * sets the neutral rate of substitutions
	 * @param rate
	 */
	public void setNeutralRate(double rate) throws IllegalArgumentException{
		if( (!Double.isNaN(rate)) && rate >= 0){
			parameters.put(BhattParameterType.FIXED_NEUTRAL_RATE, new BhattNeutralRateParameter(rate));
		}else{
			throw new IllegalArgumentException("Substitution rate must be nonnegative: " + rate);
		}
	}
	
	/**
	 * Set boolean to toggle verbose output for debugging.
	 * @param toggleDoVerboseDebug
	 */
	public void setDebugFlag(boolean toggleDoVerboseDebug){
		parameters.put(BhattParameterType.DO_VERBOSE_DEBUG, new BhattVerboseDebugParameter(toggleDoVerboseDebug));
	}
	
	/**
	 * Sets how many bootstrap replicates to use
	 * @param bootstrapReplicates
	 */
	public void setBootstrapReplicates(int bootstrapReplicates) {
		parameters.put(BhattParameterType.BOOTSTRAP_REPLICATES, new BhattBootstrapReplicatesParameter(bootstrapReplicates));	
	}
	
	public int getBootstrapReplicates(){
		return (int) parameters.get(BhattParameterType.BOOTSTRAP_REPLICATES).getParamValue();
	}

	/**
	 * Check for debug behaviour.
	 * @return boolean - true if debug verbose requested
	 */
	public boolean getDoDebugFlag(){
		return (boolean) parameters.get(BhattParameterType.DO_VERBOSE_DEBUG).getParamValue();
	}

	/**
	 * @return
	 */
	public Iterator<BhattParameterType> getParameterIterator() {
		// TODO Auto-generated method stub
		return parameters.keySet().iterator();
	}

	/**
	 * Gets a specific parameter value by type.
	 * @param param - BhattParameterType
	 * @return -  Object, which could be anything.. so assume it's toString() does something sensible
	 */
	public Object getByKey(BhattParameterType param) {
		// TODO Auto-generated method stub
		return parameters.get(param);
	}

}
