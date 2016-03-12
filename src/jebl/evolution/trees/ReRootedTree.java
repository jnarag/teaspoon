package jebl.evolution.trees;

import jebl.evolution.graphs.Node;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ReRootedTree.java 776 2007-09-05 11:17:12Z rambaut $
 */
public class ReRootedTree extends FilteredRootedTree {

    public enum RootingType {
		MID_POINT("midpoint"),
		LEAST_SQUARES("least squares");

		RootingType(String name) {
			this.name = name;
		}

		public String toString() { return name; }

		private String name;
	}

    public ReRootedTree(final RootedTree source, RootingType rootingType) {
        super(source);
	    switch (rootingType) {
		    case MID_POINT:
			break;
		    case LEAST_SQUARES:
			break;
		    default:
			    throw new IllegalArgumentException("Unknown enum value");
	    }
    }

	// PRIVATE members
//    private final Node rootChild1;
//    private final Node rootChild2;
//    private final double rootLength1;
//    private final double rootLength2;
}