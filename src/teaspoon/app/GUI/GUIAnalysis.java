/**
 * 
 */
package teaspoon.app.GUI;

import teaspoon.app.TEASPOONVersion;
import teaspoon.app.GUI.controllers.TeaspoonController;
import teaspoon.app.GUI.models.TeaspoonModel;
import teaspoon.app.GUI.views.TeaspoonView;

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
 * This is the main application entrypoint for the Teaspoon GUI app.
 * 
 * @see TeaspoonCommandLineApp
 */
public class GUIAnalysis implements Runnable{

	private final static TEASPOONVersion version = new TEASPOONVersion();
	TeaspoonView appView;
	TeaspoonModel appModel;
	TeaspoonController appController;
	
	/**
	 * Default no-arg constructor
	 */
	public GUIAnalysis(){
		appView = new TeaspoonView(version);
		appModel = new TeaspoonModel(version);
		appController = new TeaspoonController(appView,appModel, version);
	}
	
	/**
	 * Main runner method
	 */
	public void run(){}

	/**
	 * Main app method entrypoint (GUI)
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		javax.swing.SwingUtilities.invokeLater(new GUIAnalysis());
	}
}
