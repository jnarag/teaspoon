/**
 * MenuBarFactory.java
 */

package org.virion.jam.framework;

import java.awt.Toolkit;

import javax.swing.JMenuBar;

public interface MenuBarFactory {

    public final static int MENU_MASK = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();

    void populateMenuBar(JMenuBar menuBar, AbstractFrame frame);
    void deregisterMenuFactories();
    void registerPermanentMenuFactory(MenuFactory menuFactory);
    void registerMenuFactory(MenuFactory menuFactory);
}