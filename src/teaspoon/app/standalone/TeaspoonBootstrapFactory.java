/**
 * 
 */
package teaspoon.app.standalone;

import teaspoon.app.TeaspoonBootstrap;
import teaspoon.app.TeaspoonMask;
import teaspoon.app.utils.BhattAdaptationFullSiteMatrix;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 13 Jun 2018
 * @version 0.1
 */
public class TeaspoonBootstrapFactory {

	/**
	 * Contains static utility method to generate bootstrap replicates from an empirical sequence alignment.
	 */
	public TeaspoonBootstrapFactory() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * Generates the specified of bootstraps with from the seed given.
	 * For each <i>n</i>th bootstrap produced, the seed is equal to 
	 * (<i>seed + n</i>), e.g.
	 * <pre>
	 * Replicate	Seed
	 * 1	12345
	 * 2	12346
	 * 3	12347
	 * ..	...
	 * s	s+(n-1)
	 * </pre>
	 * 
	 * @param mask - a mask of alignment positions to include/exclude 
	 * @param mainPartition - the sequence alignment as int[][]
	 * @param replicates - number of replicates for each bootstram
	 * @param seed - seed to use for each of the replicates
	 * @return - An array of TeaspoonBootstrap objects, each one a wrapper for an int[] specifying how often each position should be sampled
	 */
	public static TeaspoonBootstrap[] generate(TeaspoonMask mask, BhattAdaptationFullSiteMatrix mainPartition, int numberOfReplicates, int seed) {
		TeaspoonBootstrap[] bootstraps = new TeaspoonBootstrap[numberOfReplicates];
		for(int replicate=0;replicate<bootstraps.length;replicate++){
			// increment the seed deterministically for each replicate
			int replicateSeed = seed + replicate;
			// sample a bootstrap. passes the alignment length explicitly in case mask is shorter than alignment
			bootstraps[replicate] = new TeaspoonBootstrap(mainPartition.alignmentLength(),mask,replicateSeed);
		}
		return bootstraps;
	}

}
