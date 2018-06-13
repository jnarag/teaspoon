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
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 13 Jun 2018
 * @version 0.1
 */
public class TeaspoonMask {

	public RateEstimationBehaviour estimationBehaviour;
	private double neutralRate;
	
	/**
	 * 
	 */
	public TeaspoonMask() {
		// TODO Auto-generated constructor stub
	}
	/**
	 * @param estimatedRate
	 */
	public void setNeutralRate(double estimatedRate) {
		// TODO Auto-generated method stub
		this.neutralRate = estimatedRate;
	}

}
