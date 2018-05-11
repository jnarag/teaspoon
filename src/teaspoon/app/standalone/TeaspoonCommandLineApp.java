/**
 * 
 */
package teaspoon.app.standalone;

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
 * This is the runner main-class for the command-line Teaspoon app
 * It will be the entrypoint for all command-line analyses.
 * It will be called by GUI analysis.
 * 
 * @see GUIAnalysis
 */
public class TeaspoonCommandLineApp {

	
	/**
	 * 
	 */
	public TeaspoonCommandLineApp() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// FIXME prints out all args then halts for now.
		for(String CLIargument: args){
			System.out.println("\t"+CLIargument);
		}
	}

}
