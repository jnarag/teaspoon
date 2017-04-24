package org.virion.jam.demo;

import org.virion.jam.demo.menus.DemoMenuFactory;
import org.virion.jam.framework.DefaultEditMenuFactory;
import org.virion.jam.framework.DefaultFileMenuFactory;
import org.virion.jam.framework.DefaultHelpMenuFactory;
import org.virion.jam.framework.DefaultMenuBarFactory;
import org.virion.jam.mac.MacFileMenuFactory;
import org.virion.jam.mac.MacHelpMenuFactory;
import org.virion.jam.mac.MacWindowMenuFactory;
import org.virion.jam.mac.Utils;


public class DemoMenuBarFactory extends DefaultMenuBarFactory {

	public DemoMenuBarFactory() {
		if (Utils.isMacOSX()) {
			registerMenuFactory(new MacFileMenuFactory(true));
			registerMenuFactory(new DefaultEditMenuFactory());
			registerMenuFactory(new DemoMenuFactory());
			registerMenuFactory(new MacWindowMenuFactory());
			registerMenuFactory(new MacHelpMenuFactory());
		} else {
			registerMenuFactory(new DefaultFileMenuFactory(true));
			registerMenuFactory(new DefaultEditMenuFactory());
			registerMenuFactory(new DemoMenuFactory());
			registerMenuFactory(new DefaultHelpMenuFactory());
		}
	}

}