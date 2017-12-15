package teaspoon.app.GUI.controllers;

import teaspoon.app.GUI.models.TeaspoonModel;
import teaspoon.app.GUI.views.TeaspoonView;

/**
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 5 Dec 2017
 * @version 0.0.1
 */
public class TeaspoonController {

	TeaspoonView appView;
	TeaspoonModel appModel;
	
	/**
	 * No-arg constructor. Deprecated
	 */
	@Deprecated
	public TeaspoonController(){}
	
	/**
	 * Main contstructor, requires a view and model
	 * @param appView
	 * @param appModel
	 */
	public TeaspoonController(TeaspoonView globalAppView, TeaspoonModel globalAppModel) {
		this.appView = globalAppView;
		this.appModel = globalAppModel;
	}

}
