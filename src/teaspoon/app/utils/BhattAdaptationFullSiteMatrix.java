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
		return null;
	}
	
	/**
	 * Get a consensus
	 */
	public void deriveConsensus(){}
	
	/**
	 * Return a single bootstrap replicate.
	 * Overload with a mask variant too
	 * @return
	 */
	public int[][] obtainBoostrap(){
		return null;
	}
}
