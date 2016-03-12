package virion.jam.demo;

import virion.jam.framework.*;
import virion.jam.mac.*;
import virion.jam.demo.menus.DemoMenuFactory;


public class DemoMenuBarFactory extends DefaultMenuBarFactory {

	public DemoMenuBarFactory() {
		if (virion.jam.mac.Utils.isMacOSX()) {
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