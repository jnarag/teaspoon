/**
 * 
 */
package teaspoon.app.standalone;

import java.io.File;
import java.io.FileNotFoundException;

import teaspoon.adaptation.parameters.AbstractBhattParameter;
import teaspoon.adaptation.parameters.BhattInputFileParameter;
import teaspoon.app.BhattAdaptationAnalysis;
import teaspoon.app.utils.BhattAdaptationParameters;
import teaspoon.app.utils.BhattAdaptationResults;

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
		/*
		 * Example workflow:
		 * 
		 * 1. Create a parameter list
		 * 2. Fill with parameters
		 * 3. Call one of the .runWith..() methods, returning the results
		 * 4. print the results out
		 * 
		 * (optional) Add nre parameters to list, change values
		 * (optional) re-run on the same matrix object
		 */
		BhattAdaptationParameters parameters = new BhattAdaptationParameters();
		// now populate the list
		try {
			parameters.setInputFile(null);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// run the thing
		BhattAdaptationResults results = new BhattAdaptationAnalysis(parameters).runWithFixedNR();
		// print results
		results.printToFile(new File(""));
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
