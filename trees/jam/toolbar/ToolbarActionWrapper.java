package virion.jam.toolbar;

import java.awt.event.ActionEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.Icon;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class ToolbarActionWrapper extends ToolbarAction {
    public ToolbarActionWrapper(AbstractAction action, String toolTipText, Icon icon) {
        super((String)action.getValue(Action.NAME), toolTipText, icon);
        this.action = action;

        action.addPropertyChangeListener(new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent event) {
                if (event.getPropertyName().equals("enabled")) {
                    setEnabled((Boolean)event.getNewValue());
                }
            }
        });
    }

    public void actionPerformed(ActionEvent ae) {
        action.actionPerformed(ae);
    }


    private final AbstractAction action;
}
