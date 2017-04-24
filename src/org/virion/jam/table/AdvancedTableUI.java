package org.virion.jam.table;

import java.awt.Font;
import java.awt.Point;
import java.awt.event.MouseEvent;

import javax.swing.JComponent;
import javax.swing.event.MouseInputListener;
import javax.swing.plaf.basic.BasicTableUI;

/**
 * @author rambaut
 *         Date: Oct 20, 2004
 *         Time: 10:16:52 PM
 */
public class AdvancedTableUI extends BasicTableUI {

    public void installUI(JComponent c) {
        super.installUI(c);

        if (org.virion.jam.mac.Utils.isMacOSX()) {
            c.setFont(new Font("Lucida Grande", Font.PLAIN, 9));
        }
    }

	protected MouseInputListener createMouseInputListener() {
		return new AdvancedMouseInputHandler();
	}

	class AdvancedMouseInputHandler extends MouseInputHandler {
		public void mousePressed(MouseEvent e) {
			Point origin = e.getPoint();
			int row = table.rowAtPoint(origin);
			int column = table.columnAtPoint(origin);
			if (row != -1 && column != -1) {
				if (table.isCellSelected(row, column)) {
					e.consume();
				}
			}

			super.mousePressed(e);
		}
	}
}
