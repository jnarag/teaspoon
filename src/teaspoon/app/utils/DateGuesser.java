/**
 * 
 */
package teaspoon.app.utils;

import java.io.File;
import java.util.Date;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
	
	/**
	 * Parses a date from a string, looking for currentTimeMillis
	 * @param textString
	 * @return Date
	 */
	public static Date extractDate(String textString){
		//TODO just uses current time
		//FIXME parse a long
		Date date = new Date(System.currentTimeMillis());
		return date;
	}
	
	/**
	 * Parses a string into a decimal. The last-occuring float/decimal present in the string.
	 * @param textString
	 * @return
	 */
	public static float extractFloatDate(String textString){
		//Pattern pattern = Pattern.compile("([0-9]{1,}[\\.]+[])");
		Pattern pattern = Pattern.compile("([0-9]+(\\.{1}[0-9]+)*)");
		Matcher match = pattern.matcher(textString);
		match.matches();
		float inferredDate = Float.NaN;
		int matchCount=0;
		while(match.find()){
			String substring = textString.substring(match.start(), match.end());
			System.out.println(matchCount+" match, span "+match.start()+" "+match.end()+": "+substring);
			inferredDate = Float.parseFloat(substring);
			matchCount++;
		}
		return inferredDate;
	}
	
}
