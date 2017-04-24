/*
 * Copyright (c) 2005 Biomatters LTD. All Rights Reserved.
 */

package virion.jam.framework;

public class SingleDocMenuBarFactory extends DefaultMenuBarFactory {

	public SingleDocMenuBarFactory() {
		if (virion.jam.mac.Utils.isMacOSX()) {
			registerMenuFactory(new MacFileMenuFactory(false));
			registerMenuFactory(new DefaultEditMenuFactory());
			registerMenuFactory(new MacWindowMenuFactory());
			registerMenuFactory(new MacHelpMenuFactory());
		} else {
			registerMenuFactory(new DefaultFileMenuFactory(false));
			registerMenuFactory(new DefaultEditMenuFactory());
			registerMenuFactory(new DefaultHelpMenuFactory());
		}
	}

}