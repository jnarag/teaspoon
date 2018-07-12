/**
 * 
 */
package teaspoon.app;

import org.apache.commons.math3.random.MersenneTwister;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * A bootstrap resamples positions in an alignment according to a mask_mid. 
 * The mask_mid <i>may</i> be shorter than the alignment (though not recommended)
 * but not longer or an ArrayIndexOutOfBounds exception will occur.
 * 
 * The returned vector is as long as the original alignment matrix but contains 
 * zeroes for sites that aren't sampled at all. Use TeaspoonBootstrap.numberOfValidSites()
 * to see how many sites are sampled.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 13 Jun 2018
 * @version 0.1
 */
public class TeaspoonBootstrap {

	private final int[] bootstraps;
	private final int seed;
	
	/**
	 * No-arg constructor is deprecated
	 */
	@Deprecated
	public TeaspoonBootstrap() {
		// TODO Auto-generated constructor stub
		bootstraps = null;
		seed = -1;
	}
	
	/**
	 * Constructs a new resampling across the whole alignemnt length.
	 * Uses System.currentTimeMillis() as seed.
	 * @param length
	 */
	public TeaspoonBootstrap(int length){
		// set seed
		seed = (int)System.currentTimeMillis();
		bootstraps = new int[length];
		// initialise RNG
		MersenneTwister randomNumberGenerator = new MersenneTwister(seed);
		// draw with replacement from sites
		for(int draw=0;draw<length;draw++){
			int nextInt = randomNumberGenerator.nextInt(length);
			bootstraps[nextInt]++;
		}
	}
	
	/**
	 * Constructs a new resampling across the whole alignemnt length.
	 * Uses the supplied seed.
	 * @param length
	 * @param specifiedSeed
	 */
	public TeaspoonBootstrap(int length, int specifiedSeed){
		seed = specifiedSeed;
		bootstraps = new int[length];
		// initialise RNG
		MersenneTwister randomNumberGenerator = new MersenneTwister(seed);
		// draw with replacement from sites
		for(int draw=0;draw<length;draw++){
			int nextInt = randomNumberGenerator.nextInt(length);
			bootstraps[nextInt]++;
		}
	}

	/**
	 * Constructs a new resampling across the range of sites defined by the mask_mid.
	 * Uses System.currentTimeMillis() as seed.
	 * @param length
	 * @param mask_mid
	 */
	public TeaspoonBootstrap(int length, TeaspoonMask mask){
		seed = (int)System.currentTimeMillis();
		bootstraps = new int[length];
		boolean[] maskPositions = mask.getPositions();
		// initialise RNG
		MersenneTwister randomNumberGenerator = new MersenneTwister(seed);
		// draw with replacement from sites
		for(int draw=0;draw<length;draw++){
			// only accept draws within mask_mid range
			boolean validDraw = false;
			while(!validDraw){
				// draw a new int
				int nextInt = randomNumberGenerator.nextInt(length);
				// check if int is in range
				if(maskPositions[nextInt]){
					// accept this draw
					bootstraps[nextInt]++;
					validDraw = true;
				}
			}
		}
	}
	
	/**
	 * Constructs a new resampling across the range of sites defined by the mask_mid.
	 * Uses the supplied integer seed.
	 * @param length
	 * @param mask_mid
	 * @param specifiedSeed
	 */
	public TeaspoonBootstrap(int length, TeaspoonMask mask, int specifiedSeed){
		seed = specifiedSeed;
		bootstraps = new int[length];
		boolean[] maskPositions = mask.getPositions();
		// initialise RNG
		MersenneTwister randomNumberGenerator = new MersenneTwister(seed);
		// draw with replacement from sites
		for(int draw=0;draw<length;draw++){
			// only accept draws within mask_mid range
			boolean validDraw = false;
			while(!validDraw){
				// draw a new int
				int nextInt = randomNumberGenerator.nextInt(length);
				// check if int is in range. in case mask_mid is shorter than length, test draw range first
				if(nextInt >= maskPositions.length){break;}else{
					if(maskPositions[nextInt]){
						// accept this draw
						bootstraps[nextInt]++;
						validDraw = true;
					}
				}
			}
		}
	}

	/**
	 * @return the bootstraps
	 */
	public int[] getBootstraps() {
		return bootstraps;
	}

	/**
	 * @return the seed
	 */
	public int getSeed() {
		return seed;
	}

	/**
	 * total number of alignment positions this bootstrap will sample from, e.g.
	 * <pre>
	 * {0,2,1} 	= 3
	 * {3,8,1} 	= 12
	 * {0,0,0} 	= 0
	 * {1,1,1} 	= 3
	 * </pre>	 
	 * @return total number of alignment positions this bootstrap will sample from.
	 */
	public int getNumberOfValidPositions() {
		// TODO Auto-generated method stub
		int numberOfValidPositions = 0;
		for(int site:this.bootstraps){
			numberOfValidPositions = numberOfValidPositions + site;
		}
		return numberOfValidPositions;
	}
}
