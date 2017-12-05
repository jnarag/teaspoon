/**
 * 
 */
package teaspoon.app.GUI.views;

import java.awt.GraphicsConfiguration;
import java.awt.HeadlessException;

import javax.swing.JFrame;

/**
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 5 Dec 2017
 * @version 0.0.1
 */
public class TeaspoonView extends JFrame {

	/**
	 * Default no-arg constructor
	 * @throws HeadlessException
	 */
	public TeaspoonView() throws HeadlessException {
		super();
		this.setTitle("Teaspoon");
		this.setDefaultCloseOperation(EXIT_ON_CLOSE);
		this.setSize(720, 320);
		this.setVisible(true);
	}

	/**
	 * @param gc
	 */
	public TeaspoonView(GraphicsConfiguration gc) {
		super(gc);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param title
	 * @throws HeadlessException
	 */
	public TeaspoonView(String title) throws HeadlessException {
		super(title);
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param title
	 * @param gc
	 */
	public TeaspoonView(String title, GraphicsConfiguration gc) {
		super(title, gc);
		// TODO Auto-generated constructor stub
	}

}
