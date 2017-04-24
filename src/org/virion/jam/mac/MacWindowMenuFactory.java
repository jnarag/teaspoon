package org.virion.jam.mac;

import java.awt.event.KeyEvent;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

import org.virion.jam.framework.AbstractFrame;
import org.virion.jam.framework.Application;
import org.virion.jam.framework.MenuBarFactory;
import org.virion.jam.framework.MenuFactory;

/**
 * @author rambaut
 *         Date: Dec 26, 2004
 *         Time: 11:03:39 AM
 */
public class MacWindowMenuFactory implements MenuFactory {
	public String getMenuName() {
		return "Window";
	}

	public void populateMenu(JMenu menu, AbstractFrame frame) {

		Application application = Application.getApplication();

		JMenuItem item;

		item = new JMenuItem(frame.getZoomWindowAction());
		menu.add(item);

		item = new JMenuItem(frame.getMinimizeWindowAction());
		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_M, MenuBarFactory.MENU_MASK));
		menu.add(item);

	}

	public int getPreferredAlignment() {
		return RIGHT;
	}
}
