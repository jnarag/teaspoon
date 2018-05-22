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
 * Custom Exception thrown if you try to run a counting method based on fixed neutral ratio, but pass a null value.
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 22 May 2018
 * @version 0.1.2
 */
public class NullNeutralRatioException extends Exception {

	/**
	 * 
	 */
	public NullNeutralRatioException() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param message
	 */
	public NullNeutralRatioException(String message) {
		super(message);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param cause
	 */
	public NullNeutralRatioException(Throwable cause) {
		super(cause);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param message
	 * @param cause
	 */
	public NullNeutralRatioException(String message, Throwable cause) {
		super(message, cause);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param message
	 * @param cause
	 * @param enableSuppression
	 * @param writableStackTrace
	 */
	public NullNeutralRatioException(String message, Throwable cause,
			boolean enableSuppression, boolean writableStackTrace) {
		super(message, cause, enableSuppression, writableStackTrace);
		// TODO Auto-generated constructor stub
	}

}
