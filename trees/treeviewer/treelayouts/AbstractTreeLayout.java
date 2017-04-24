package jebl.gui.trees.treeviewer.treelayouts;

import java.awt.Shape;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import jebl.evolution.trees.Tree;

/**
 * @author Andrew Rambaut
 * @version $Id: AbstractTreeLayout.java 181 2006-01-23 17:31:10Z rambaut $
 */
public abstract class AbstractTreeLayout implements TreeLayout {
    public void setTree(Tree tree) {
        this.tree = (jebl.evolution.trees.RootedTree)tree;
	    invalidate();
    }

    public void invalidate() {
        invalid = true;
        fireTreeLayoutChanged();
    }

    public Point2D getNodePoint(jebl.evolution.graphs.Node node) {
        checkValidation();
        return nodePoints.get(node);
    }

    public Shape getBranchPath(jebl.evolution.graphs.Node node) {
        checkValidation();
        return branchPaths.get(node);
    }

    public Line2D getTaxonLabelPath(jebl.evolution.graphs.Node node) {
        checkValidation();
        return taxonLabelPaths.get(node);
    }

    public Line2D getBranchLabelPath(jebl.evolution.graphs.Node node) {
        checkValidation();
        return branchLabelPaths.get(node);
    }

    public Line2D getNodeLabelPath(jebl.evolution.graphs.Node node) {
        checkValidation();
        return nodeLabelPaths.get(node);
    }

	public Shape getCalloutPath(jebl.evolution.graphs.Node node) {
	    checkValidation();
	    return calloutPaths.get(node);
	}

    private void checkValidation() {
        if (invalid) {
            validate();
            invalid = false;
        }
    }

    public void addTreeLayoutListener(TreeLayoutListener listener) {
        listeners.add(listener);
    }

    public void removeTreeLayoutListener(TreeLayoutListener listener) {
        listeners.remove(listener);
    }

    protected void fireTreeLayoutChanged() {
        for (TreeLayoutListener listener : listeners) {
            listener.treeLayoutChanged();
        }
    }

    protected abstract void validate();

    private boolean invalid = true;
    protected jebl.evolution.trees.RootedTree tree = null;
    protected Map<jebl.evolution.graphs.Node, Point2D> nodePoints = new HashMap<jebl.evolution.graphs.Node, Point2D>();
    protected Map<jebl.evolution.graphs.Node, Shape> branchPaths = new HashMap<jebl.evolution.graphs.Node, Shape>();
    protected Map<jebl.evolution.graphs.Node, Line2D> taxonLabelPaths = new HashMap<jebl.evolution.graphs.Node, Line2D>();
    protected Map<jebl.evolution.graphs.Node, Line2D> branchLabelPaths = new HashMap<jebl.evolution.graphs.Node, Line2D>();
    protected Map<jebl.evolution.graphs.Node, Line2D> nodeLabelPaths = new HashMap<jebl.evolution.graphs.Node, Line2D>();
	protected Map<jebl.evolution.graphs.Node, Shape> calloutPaths = new HashMap<jebl.evolution.graphs.Node, Shape>();

    private Set<TreeLayoutListener> listeners = new HashSet<TreeLayoutListener>();
}
