/**
 * 
 */
package teaspoon.adaptation.parameters;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * Bhatt Parameter Type defines the type of parameters used in a Bhatt adaptayion analysis.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public enum BhattParameterType {
	INPUT_FILE_ANCESTRAL,
	INPUT_FILE_MAIN,
	INPUT_FILE_LIST_MAIN,
	INPUT_FILE_MASK,
	OUTPUT_FILE,
	FIXED_NEUTRAL_RATE,
	DO_VERBOSE_DEBUG,
	BOOTSTRAP_REPLICATES, 
	INT_MATRIX_ANCESTRAL,
	INT_MATRIX_MAIN,
	SITE_FREQUENCY_BINS, 
	RUN_ID
}
