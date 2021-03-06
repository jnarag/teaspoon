/**
 * 
 */
package teaspoon.app.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

import teaspoon.adaptation.parameters.AbstractBhattParameter;
import teaspoon.adaptation.parameters.BhattAncestralFullSiteMatrixParameter;
import teaspoon.adaptation.parameters.BhattBootstrapReplicatesParameter;
import teaspoon.adaptation.parameters.BhattMainInputFileListParameter;
import teaspoon.adaptation.parameters.BhattMainInputFileParameter;
import teaspoon.adaptation.parameters.BhattMaskInputFileParameter;
import teaspoon.adaptation.parameters.BhattNeutralRateParameter;
import teaspoon.adaptation.parameters.BhattParameterType;
import teaspoon.adaptation.parameters.BhattRunIDParameter;
import teaspoon.adaptation.parameters.BhattSiteFreqBinDefinitionsParameter;
import teaspoon.adaptation.parameters.BhattVerboseDebugParameter;
import teaspoon.adaptation.parameters.BhattAncestralFullSiteMatrixParameter;
import teaspoon.adaptation.parameters.BhattMainFullSiteMatrixParameter;

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
	 * @param output
	 * @return
	 * @throws FileNotFoundException
	 */
	public void setOutputFile(File output){
		parameters.put(BhattParameterType.OUTPUT_FILE, new BhattMainInputFileParameter(output));
	}
	
	/**
	 * Return the input file value
	 * @return - an input file
	 */
	public File getOutputFile(){
		// Walk through the parameters
		return (File) parameters.get(BhattParameterType.OUTPUT_FILE).getParamValue();
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
	 * Set or update the mask_mid file
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
	 * Return the mask_mid file value
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
			parameters.put(BhattParameterType.FIXED_NEUTRAL_RATE, new BhattNeutralRateParameter(0.0));
			throw new IllegalArgumentException("Substitution rate must be nonnegative: " + rate);
		}
	}

	/**
	 * @param subsampleByMask
	 */
	public void setAncestralFullSiteMatrix(BhattAdaptationFullSiteMatrix ancestralMatrix) {
		parameters.put(BhattParameterType.INT_MATRIX_ANCESTRAL, new BhattAncestralFullSiteMatrixParameter(ancestralMatrix));
	}


	/**
	 * @param subsampleByMask
	 */
	public void setInputFullSiteMatrix(BhattAdaptationFullSiteMatrix mainMatrix) {
		parameters.put(BhattParameterType.INT_MATRIX_MAIN, new BhattMainFullSiteMatrixParameter(mainMatrix));
		
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

	/**
	 * @return true if the FSM has non-zero length
	 */
	public boolean hasFullSiteMatrixMain() {
		try{
			if( ((BhattAdaptationFullSiteMatrix) parameters.get(BhattParameterType.INT_MATRIX_MAIN).getParamValue()).alignmentLength() > 0){
				return true;
			}else{
				return false;
			}		
		}catch (NullPointerException ex){
			ex.printStackTrace();
			return false;
		}
	}

	/**
	 * @return
	 */
	public boolean hasFullSiteMatrixAncestral() {
		try {
			if( ((BhattAdaptationFullSiteMatrix) parameters.get(BhattParameterType.INT_MATRIX_ANCESTRAL).getParamValue()).alignmentLength() > 0){
				return true;
			}else{
				return false;
			}
		}catch (NullPointerException ex){
			// TODO Auto-generated catch block
			ex.printStackTrace();
			return false;
		}
	}

	/**
	 * @return
	 */
	public BhattAdaptationFullSiteMatrix getFullSiteMatrixMain() {
		// TODO Auto-generated method stub
		return (BhattAdaptationFullSiteMatrix) parameters.get(BhattParameterType.INT_MATRIX_MAIN).getParamValue();
	}

	/**
	 * @return
	 */
	public BhattAdaptationFullSiteMatrix getFullSiteMatrixAncestral() {
		// TODO Auto-generated method stub
		return (BhattAdaptationFullSiteMatrix) parameters.get(BhattParameterType.INT_MATRIX_ANCESTRAL).getParamValue();
	}

	/**
	 * Tests whether there is a double[2][3] matrix describing the bins
	 * but doesn't check for NaNs or negative values in the matrix itself.
	 * 
	 * @return true if there seems to be custom bins boundaries here
	 */
	public boolean hasCustomBinSettings() {
		if( 
			((double[][]) parameters.get(BhattParameterType.SITE_FREQUENCY_BINS).getParamValue())[0].length == 3 &&
			((double[][]) parameters.get(BhattParameterType.SITE_FREQUENCY_BINS).getParamValue()).length == 2
		){
			return true;
		}else{
			return false;			
		}
	}
	
	/**
	 * Returns the custom bin settings, should be an double[3][2] matrix
	 * @return
	 * @see BhattSiteFreqBinDefinitionsParameter
	 */
	public double[][] getCustomBinSettings(){
		return ((double[][]) parameters.get(BhattParameterType.SITE_FREQUENCY_BINS).getParamValue());
	}

	/**
	 * Set the custom bin settings, should be an double[3][2] matrix
	 * @see BhattSiteFreqBinDefinitionsParameter
	 */
	public void setCustomBinSettings(double[][] binIntervals){
		this.parameters.put(BhattParameterType.SITE_FREQUENCY_BINS, new BhattSiteFreqBinDefinitionsParameter(binIntervals));
	}

	/**
	 * @return the run ID as String
	 */
	public String getRunID() {
		return ((String) parameters.get(BhattParameterType.RUN_ID).getParamValue());
	}

	/**
	 * Set the custom run ID settings, should be a String
	 * @see BhattRunIDParameter
	 */
	public void setRunID(String run_ID) {
		this.parameters.put(BhattParameterType.RUN_ID, new BhattRunIDParameter(run_ID));
	}
}
