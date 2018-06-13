/**
 * 
 */
package teaspoon.app.utils;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * Defines the neutral rate estimation behaviour for a mask (partition):
 * <pre>
  	NEUTRAL_RATE_FIXED		Rate is fixed (positive double)
	NEUTRAL_RATE_AGGREGATED	Rate estimated by aggregating site frequencies over all main alignments
	NEUTRAL_RATE_AVERAGED	Rates estimated separately for all main alignments, then averaged.
 * </pre>
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 13 Jun 2018
 * @version 0.1
 */
public enum RateEstimationBehaviour {
	NEUTRAL_RATE_FIXED,
	NEUTRAL_RATE_AGGREGATED,
	NEUTRAL_RATE_AVERAGED
}
