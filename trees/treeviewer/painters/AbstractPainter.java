package jebl.gui.trees.treeviewer.painters;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id: AbstractPainter.java 181 2006-01-23 17:31:10Z rambaut $
 */
public abstract class AbstractPainter<T> implements Painter<T> {
    public void addPainterListener(jebl.gui.trees.treeviewer.painters.PainterListener listener) {
        listeners.add(listener);
    }

    public void removePainterListener(jebl.gui.trees.treeviewer.painters.PainterListener listener) {
        listeners.remove(listener);
    }

    public void firePainterChanged() {
        for (jebl.gui.trees.treeviewer.painters.PainterListener listener : listeners) {
            listener.painterChanged();
        }
    }
    private final List<jebl.gui.trees.treeviewer.painters.PainterListener> listeners = new ArrayList<jebl.gui.trees.treeviewer.painters.PainterListener>();
}
