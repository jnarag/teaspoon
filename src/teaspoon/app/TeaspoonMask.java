/**
 * 
 */
package teaspoon.app;

import teaspoon.app.utils.RateEstimationBehaviour;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * Key parameters:
 * 	- estimationBehaviour enum RateEstimationBehaviour one of 'fixed', 'aggregated', 'averaged'. (default: aggregated)
 *  - neutral ratio, double
 *  - maskValues, boolean[] representing the alignment positions. 'true' to use in the analysis.
 *  
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 13 Jun 2018
 * @version 0.1
 */
public class TeaspoonMask {

	public RateEstimationBehaviour estimationBehaviour;
	private double neutralRatio;
	private final boolean[] maskValues;
	
	/**
	 * No-arg constructor is deprecated
	 */
	@Deprecated
	public TeaspoonMask() {
		// TODO Auto-generated constructor stub
		maskValues = null;
		neutralRatio = 0;
		estimationBehaviour = null;
	}
	
	public TeaspoonMask(RateEstimationBehaviour behaviour, boolean[] mask){
		maskValues = mask;
		estimationBehaviour = behaviour;
	}
	
	
	public TeaspoonMask(RateEstimationBehaviour behaviour, boolean[] mask, double rate){
		maskValues = mask;
		estimationBehaviour = behaviour;
		neutralRatio = rate;
	}

	/**
	 * Get the neutral ratio if available. Warning! If using ratio estimation, there may be a null or zero ratio.
	 * @return
	 */
	public double getNeutralRatio(){
		return this.neutralRatio;
	}
	
	/**
	 * @param estimatedRatio
	 */
	public void setNeutralRatio(double estimatedRatio) {
		// TODO Auto-generated method stub
		this.neutralRatio = estimatedRatio;
	}
	
	/**
	 * Returns the total length of the sequence
	 * @return
	 */
	public int getLength(){
		return this.maskValues.length;
	}

	/**
	 * @return the boolean[] of mask positions (true=include)
	 */
	public boolean[] getPositions() {
		// TODO Auto-generated method stub
		return this.maskValues;
	}

	/**
	 * total number of alignment positions this maks will sample from, e.g.
	 * <pre>
	 * {false,true,true} 	= 2
	 * {true,true,true} 	= 3
	 * {false,false,false} 	= 0
	 * </pre>	 
	 * @return total number of alignment positions this maks will sample from.
	 */
	public int getNumberOfValidPositions() {
		// TODO Auto-generated method stub
		int numberOfValidPositions = 0;
		for(boolean site:this.maskValues){
			if(site){numberOfValidPositions++;}
		}
		return numberOfValidPositions;
	}

	/**
	 * Gets the first alignment position included in the mask (indexed to zero).
	 * Where mask regions are discontinuous this will be the 5-prime start.
	 * @return first position included in the mask
	 */
	public int getFirstStart() {
		int firstStart = -1;
		for(int i=0;i<this.maskValues.length;i++){
			if(maskValues[i] && firstStart < 0){
				firstStart = i;
			}
		}
		return firstStart;
	}

	/**
	 * Gets the last alignment position included in the mask (indexed to zero).
	 * Where mask regions are discontinuous this will be the 3-prime end.
	 * @return last position included in the mask
	 */
	public int getLastEnd() {
		int lastEnd = maskValues.length+1;
		for(int i=this.maskValues.length-1;i>-1;i--){
			if(maskValues[i] && lastEnd>maskValues.length){
				lastEnd = i;
			}
		}
		return lastEnd;
	}

}
