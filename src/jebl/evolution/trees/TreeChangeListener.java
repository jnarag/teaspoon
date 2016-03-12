package jebl.evolution.trees;

/**
 * A listener for notifying about changes to a tree or to the selected nodes in a tree.
 *
 * @author Matt Kearse
 * @version $Id: TreeChangeListener.java 913 2008-05-22 04:53:39Z matt_kearse $
 */

public abstract class TreeChangeListener {
    /**
     * The tree has changed. The tree contained in the TreeChangeEvent must not
     * be the original tree. Instead it must be a new instance of a tree,
     * first cloned using {@link Utils#copyTree(RootedTree)} and {@link Utils#rootTheTree(Tree)}  or {@link Utils#rootTreeAtCenter(Tree)} if necessary
     * before changes are made.
     *
     * @param treeChangeEvent the changed tree.
     */
    public abstract void treeChanged(TreeChangeEvent treeChangeEvent);

    /**
     * The selected nodes in the tree have changed.
     * @param treeChangeEvent the new set of selected nodes.
     */
    public abstract void selectionChanged(TreeSelectionChangeEvent treeChangeEvent);
}
