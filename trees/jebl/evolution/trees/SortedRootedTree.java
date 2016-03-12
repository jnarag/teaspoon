package jebl.evolution.trees;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: SortedRootedTree.java 627 2007-01-15 03:50:40Z pepster $
 */
public class SortedRootedTree extends FilteredRootedTree {

    public enum BranchOrdering {
		INCREASING_NODE_DENSITY("increasing"),
		DECREASING_NODE_DENSITY("decreasing");

		BranchOrdering(String name) {
			this.name = name;
		}

		public String toString() { return name; }

		private String name;
	}

    public SortedRootedTree(final RootedTree source, BranchOrdering branchOrdering) {
        super(source);
	    switch (branchOrdering) {
		    case INCREASING_NODE_DENSITY:
			    this.comparator = new Comparator<jebl.evolution.graphs.Node>() {
			        public int compare(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) {
			            return jebl.evolution.trees.Utils.getExternalNodeCount(source, node1) -
					            jebl.evolution.trees.Utils.getExternalNodeCount(source, node2);
			        }

			        public boolean equals(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) {
			            return compare(node1, node2) == 0;
			        }
			    };
			break;
		    case DECREASING_NODE_DENSITY:
			    this.comparator = new Comparator<jebl.evolution.graphs.Node>() {
			        public int compare(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) {
			            return jebl.evolution.trees.Utils.getExternalNodeCount(source, node2) -
					            jebl.evolution.trees.Utils.getExternalNodeCount(source, node1);
			        }

			        public boolean equals(jebl.evolution.graphs.Node node1, jebl.evolution.graphs.Node node2) {
			            return compare(node1, node2) == 0;
			        }
			    };
			break;
		    default:
			    throw new IllegalArgumentException("Unknown enum value");
	    }
    }

    public SortedRootedTree(RootedTree source, Comparator<jebl.evolution.graphs.Node> comparator) {
	    super(source);
        this.comparator = comparator;
    }

    public List<jebl.evolution.graphs.Node> getChildren(jebl.evolution.graphs.Node node) {
        List<jebl.evolution.graphs.Node> sourceList = source.getChildren(node);
        Collections.sort(sourceList, comparator);
        return sourceList;
    }

	// PRIVATE members

    private final Comparator<jebl.evolution.graphs.Node> comparator;
}