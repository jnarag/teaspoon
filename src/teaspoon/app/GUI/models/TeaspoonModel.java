/**
 * 
 */
package teaspoon.app.GUI.models;

/**
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 5 Dec 2017
 * @version 0.0.1
 */
public class TeaspoonModel {

	/*
	 * Data model
	 * @TODO checks and other sensible things
	 */
	private Object[][] data; 
	
	/**
	 * Default no-arg constructor
	 */
	public TeaspoonModel(){
		this.setData(new Object[1][1]);
	}

	/**
	 * @return the data
	 */
	public Object[][] getData() {
		return data;
	}

	/**
	 * @param data the data to set
	 */
	private void setData(Object[][] data) {
		this.data = data;
	}
}
