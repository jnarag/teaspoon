package jebl.gui.trees.treeviewer;

import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;

import jebl.evolution.trees.Tree;
import jebl.evolution.trees.Utils;
import jebl.gui.trees.treeviewer.painters.Painter;

/**
 * @author Joseph Heled
 * @version $Id$
 *          <p/>
 *          Created by IntelliJ IDEA.
 *          User: joseph
 *          Date: 19/12/2006
 *          Time: 12:38:58
 */
public class TreeDrawableElementNodeLabel extends TreeDrawableElementLabel {

    private Rectangle2D defaultBounds;
    private Rectangle2D minBounds = null;
    // when size for node lable is to be taken from another node, known to be larger
    private jebl.evolution.graphs.Node nodeSizeReference;
    private jebl.gui.trees.treeviewer.painters.BasicLabelPainter painter;
    private int defaultSize;
    private int curSize;
    private Painter.Justification taxonLabelJustification;
    //private AffineTransform save;
    // debug
    private Tree tree;
    String dtype;

    TreeDrawableElementNodeLabel(Tree tree, jebl.evolution.graphs.Node node, Painter.Justification taxonLabelJustification,
                                  Rectangle2D labelBounds, AffineTransform transform, int priority,
                                  jebl.evolution.graphs.Node nodeSizeReference, jebl.gui.trees.treeviewer.painters.BasicLabelPainter painter,
                                  String dtype) {
        super(node, labelBounds, transform, priority);

        defaultBounds = labelBounds;
        this.nodeSizeReference = nodeSizeReference != null ? nodeSizeReference : node;
        //this.node = node;
        this.taxonLabelJustification = taxonLabelJustification;
        this.painter = painter;
        curSize = defaultSize = (int)painter.getFontSize();

        //save = new AffineTransform(transform);

        this.tree = tree;
        this.dtype = dtype;
    }

    public void setSize(int size, Graphics2D g2) {
        if( curSize != size ) {
            // System.out.println("Set size of " + getDebugName() + " to " + size);
            Rectangle2D newBounds;
            if( size == getMaxSize() ) {
                newBounds = defaultBounds;
            } else if( size == getMinSize() && minBounds != null ) {
                newBounds = minBounds;
            } else {
                // set it up
                // do it more effciently, share between all
                float s = painter.getFontSize();
                painter.setFontSize(size, false);
                painter.calibrate(g2);
                newBounds = new Rectangle2D.Double(0.0, 0.0, painter.getWidth(g2, nodeSizeReference), painter.getPreferredHeight());
                if( size == getMinSize() ) {
                    minBounds = newBounds;
                }
                painter.setFontSize(s, false);
                painter.calibrate(g2);
            }        

            if(prints>1) System.out.println("before " + getDebugName() + " at " + curSize + " to " + size  + " transorm " + transform);

            final double dx = newBounds.getWidth() - bounds.getWidth();
            final double dy = newBounds.getHeight() - bounds.getHeight();

            bounds = newBounds;
            if( taxonLabelJustification != Painter.Justification.CENTER ) {
                transform.translate(taxonLabelJustification == Painter.Justification.RIGHT ? -dx : 0 , -dy/2);
            } else {
                transform.translate(-dx/2, -dy);
            }
            if(prints>1) System.out.println("change size of " + getDebugName() + " to " + size  + " transorm " + transform);
            curSize = size;
        }
    }

    public String getDebugName() {
        String name;
        if( tree.isExternal(node) ) {
            name = tree.getTaxon(node).getName();

        } else {
            name = Utils.DEBUGsubTreeRep(Utils.rootTheTree(tree), node);
        }
        return name + (dtype != null ? "(" + dtype + ")" : "");
    }

    public int getCurrentSize() {
        return curSize;
    }

    protected void drawIt(Graphics2D g2) {
        AffineTransform oldTransform = g2.getTransform();

        g2.transform(transform);
        float s = painter.getFontSize();
        if( painter.setFontSize(curSize, false) ) {
            painter.calibrate(g2);
        }
        painter.setForeground(foreground);
        painter.paint(g2, node, taxonLabelJustification, bounds);

        if( painter.setFontSize(s, false) ) {
            painter.calibrate(g2);
        }
        g2.setTransform(oldTransform);
    }

    public int getMinSize() {
        return Math.min((int)painter.getFontMinSize(), defaultSize);
    }

    public int getMaxSize() {
        return defaultSize;
    }
}