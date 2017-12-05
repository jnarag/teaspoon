/**
 * 
 */
package teaspoon.app.GUI;

import teaspoon.app.GUI.controllers.TeaspoonController;
import teaspoon.app.GUI.models.TeaspoonModel;
import teaspoon.app.GUI.views.TeaspoonView;

/**
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 5 Dec 2017
 * @version 0.0.1
 * 
 * This is the main application entrypoint for the Teaspoon GUI app.
 * 
 * @see TeaspoonCommandLineApp
 */
public class GUIAnalysis {

	TeaspoonView appView;
	TeaspoonModel appModel;
	TeaspoonController appController;
	
	/**
	 * Default no-arg constructor
	 */
	public GUIAnalysis(){
		appView = new TeaspoonView();
		appModel = new TeaspoonModel();
		appController = new TeaspoonController(appView,appModel);
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
		new GUIAnalysis().run();
	}
}
