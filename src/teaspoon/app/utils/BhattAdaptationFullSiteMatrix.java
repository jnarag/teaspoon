/**
 * 
 */
package teaspoon.app.utils;

import java.io.File;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * The 'full' site matrix contains all positions in an input main- or ancestral-file
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public class BhattAdaptationFullSiteMatrix {

	// holds all the data for a file //
	int[][] siteMatrix;
	
	/**
	 * constructor takes an 
	 * @param readFASTA
	 */
	public BhattAdaptationFullSiteMatrix(int[][] fasta) {
		// TODO Auto-generated constructor stub
		siteMatrix = fasta;
	}

	/**
	 * @return the siteMatrix
	 */
	public int[][] getSiteMatrix() {
		return siteMatrix;
	}

	/**
	 * @param siteMatrix the siteMatrix to set
	 */
	public void setSiteMatrix(int[][] siteMatrix) {
		this.siteMatrix = siteMatrix;
	}

	/**
	 * Loads input to int matrix
	 */
	public void loadAlignmentFile(File input){}
	
	/**
	 * Does all the clean-up nonsense
	 */
	public void cleanUpSites(){}
	
	/**
	 * Given a list of mask positions, returns submatrix
	 * @param mask
	 * @return
	 */
	public int[][] subsampleByMask(boolean[] mask){
		//TODO implement subslice of positions
		return siteMatrix;
	}
	
	/**
	 * Get a consensus
	 */
	public int[] deriveConsensus(){
		//TODO implement as per BhattMethod
		return MainAlignmentParser.consensusArray(siteMatrix);
	}
	
	/**
	 * Return a single bootstrap replicate.
	 * Overload with a mask variant too
	 * @return
	 */
	public int[][] obtainBoostrap(){
		return null;
	}
}
