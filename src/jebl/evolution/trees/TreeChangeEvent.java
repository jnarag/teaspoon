package jebl.evolution.trees;

/**
 * Represents a change to a tree. As at 2008-5-22, this class just encapsulates the tree
 * but in future it may contain more details on the type of change made to the tree.
 * @author Matt Kearse
 * @version $Id: TreeChangeEvent.java 913 2008-05-22 04:53:39Z matt_kearse $
 */

public final class TreeChangeEvent {
    private Tree tree;

    public TreeChangeEvent(Tree tree){
        this.tree = tree;
    }

    public Tree getTree(){
        return tree;
    }
}
