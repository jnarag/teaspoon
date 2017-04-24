/**
* ConsoleMenuBarFactory.java
*/

package virion.jam.console;

public class ConsoleMenuBarFactory extends DefaultMenuBarFactory {

	public ConsoleMenuBarFactory() {
        // org.virion stuff shouldn't be called from here - it's a separate project!

		// no its not. This class is part of JAM.
        if (virion.jam.mac.Utils.isMacOSX()) {
        //if (System.getProperty("mrj.version") != null) {
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