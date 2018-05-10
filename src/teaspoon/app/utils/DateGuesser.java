/**
 * 
 */
package teaspoon.app.utils;

import java.io.File;
import java.util.Date;
import java.util.HashMap;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 5 Dec 2017
 * @version 0.0.1
 * 
 * Guesses dates from sequence file IDs
 * Example usage:
 * <pre>
 * try{
 *  // Try and get dates from the alignment
 * 	HashMap<String,Date> sequenceDates = new DateGuesser(input_alignment).getDates()
 * }catch (Exception ex){
 *  // throws a not-very-useful Exception if anything goes wrong.
 * 	ex.printStackTrace();
 * }
 * </pre>
 */
public class DateGuesser {

	private HashMap<String, Date> dates;
	
	/**
	 * @return the dates
	 */
	public HashMap<String, Date> getDates() {
		return dates;
	}

	/**
	 * Default no-arg constructor. Deprecated
	 */
	@Deprecated
	public DateGuesser() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * Constructor
	 * Guesses dates from an alignmentfile and populates the dates object
	 * @param alignmentFile
	 */
	public DateGuesser(File alignmentFile) throws Exception{
		// @TODO
	}
}
