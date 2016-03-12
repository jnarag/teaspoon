package jebl.evolution.trees;

import jebl.evolution.graphs.Node;

import java.util.Set;

/**
 * Represents a change to the selected nodes in a tree. As at 2008-5-22, this class just encapsulates the set
 * of selected nodes but in future it may contain more details on the type of selection change.
 * @author Matt Kearse
 * @version $Id: TreeSelectionChangeEvent.java 913 2008-05-22 04:53:39Z matt_kearse $
 */

public final class TreeSelectionChangeEvent {
    private Set<Node> selectedNodes;

    public TreeSelectionChangeEvent(Set<Node> selectedNodes){
        this.selectedNodes =   selectedNodes;
    }

    public Set<Node> getSelectedNodes(){
        return selectedNodes;
    }

}
