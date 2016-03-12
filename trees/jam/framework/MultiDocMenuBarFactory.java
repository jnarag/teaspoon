/*
 * Copyright (c) 2005 Biomatters LTD. All Rights Reserved.
 */

package virion.jam.framework;

import virion.jam.mac.MacFileMenuFactory;
import virion.jam.mac.MacHelpMenuFactory;
import virion.jam.mac.MacWindowMenuFactory;


public class MultiDocMenuBarFactory extends DefaultMenuBarFactory {


	public MultiDocMenuBarFactory() {
		if (virion.jam.mac.Utils.isMacOSX()) {
			registerMenuFactory(new MacFileMenuFactory(true));
			registerMenuFactory(new DefaultEditMenuFactory());
			registerMenuFactory(new MacHelpMenuFactory());
			registerMenuFactory(new MacWindowMenuFactory());
		} else {
			registerMenuFactory(new DefaultFileMenuFactory(true));
			registerMenuFactory(new DefaultEditMenuFactory());
			registerMenuFactory(new DefaultHelpMenuFactory());
		}
	}
}