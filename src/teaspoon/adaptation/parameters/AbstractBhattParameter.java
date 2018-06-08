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
 * Abstract Bhatt Parameter defines values/types of parameter.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 7 Jun 2018
 * @version 0.1
 */
public abstract class AbstractBhattParameter {

	protected BhattParameterType paramType;
	protected Object paramValue;
	
	/**
	 * No-arg constructor
	 */
	public AbstractBhattParameter(){
		paramType = null;
		paramValue = null;
	}
	
	/**
	 * Value-only parameter for use by concrete subclasses tied to a type e.g. BhattMainInputFileParameter etc
	 * @param value
	 */
	public AbstractBhattParameter(Object value) {
		paramType = null;
		paramValue = value;
	}
	
	/**
	 * Totally explicit constructor for general purpose concrete subclass, relies on a sensible param and value.
	 * @param type
	 * @param value
	 */
	public AbstractBhattParameter(BhattParameterType type, Object value) {
		paramType = type;
		paramValue = value;
	}

	/**
	 * @return the paramType
	 */
	public BhattParameterType getParamType() {
		return paramType;
	}

	/**
	 * @return the paramValue
	 */
	public Object getParamValue() {
		return paramValue;
	}

}
