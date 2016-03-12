package jebl.evolution.trees;

import java.util.*;

/**
 * Static utility functions for rooted trees.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: RootedTreeUtils.java 889 2008-02-27 01:13:21Z matt_kearse $
 */

public class RootedTreeUtils {
    private RootedTreeUtils() { }  // make class uninstantiable

    /**
	 * Return the number of leaves under this node.
	 * @param tree
	 * @param node
	 * @return the number of leaves under this node.
	 */
	public static final int getTipCount(RootedTree tree, jebl.evolution.graphs.Node node) {
		int tipCount = 0;
		for (jebl.evolution.graphs.Node child : tree.getChildren(node)) {
			tipCount += getTipCount(tree, child);
		}

		// is external
		if (tipCount == 0) return 1;

		return tipCount;
	}

    public static double getMinTipHeight(RootedTree tree, jebl.evolution.graphs.Node node) {
        if (tree.isExternal(node)) {
            return tree.getHeight(node);
        }

        double minTipHeight = Double.MAX_VALUE;
        for (jebl.evolution.graphs.Node child : tree.getChildren(node)) {
           double h = getMinTipHeight(tree, child);
            if (h < minTipHeight) {
                minTipHeight = h;
            }
        }

        return minTipHeight;
    }

    public static double getMaxTipHeight(RootedTree tree, jebl.evolution.graphs.Node node) {
        if (tree.isExternal(node)) {
            return tree.getHeight(node);
        }

        double maxTipHeight = -Double.MAX_VALUE;
        for (jebl.evolution.graphs.Node child : tree.getChildren(node)) {
           double h = getMaxTipHeight(tree, child);
            if (h > maxTipHeight) {
                maxTipHeight = h;
            }
        }

        return maxTipHeight;
    }

    /**
	 * @return true only if all tips have height 0.0
	 */
	public static boolean isUltrametric(RootedTree tree, double tolerance) {
		for (jebl.evolution.graphs.Node node : tree.getExternalNodes()) {
			if (Math.abs(tree.getHeight(node)) > tolerance) return false;
		}
		return true;
	}

	/**
	 * @return true only if internal nodes have 2 children
	 */
	public static boolean isBinary(RootedTree tree) {
		for (jebl.evolution.graphs.Node node : tree.getInternalNodes()) {
			if (tree.getChildren(node).size() > 2) return false;
		}
		return true;
	}

	/**
	 * Gets a set of external nodes that correspond to the given taxa.
	 */
	public static Set<jebl.evolution.graphs.Node> getTipsForTaxa(RootedTree tree, Collection<jebl.evolution.taxa.Taxon> taxa) throws jebl.evolution.taxa.MissingTaxonException {

		Set<jebl.evolution.graphs.Node> tipNodes = new LinkedHashSet<jebl.evolution.graphs.Node>();

		for (jebl.evolution.taxa.Taxon taxon : taxa) {

			jebl.evolution.graphs.Node node = tree.getNode(taxon);

			if (node == null) {
				throw new jebl.evolution.taxa.MissingTaxonException(taxon);
			}

			tipNodes.add(node);
		}

		return tipNodes;
	}

	/**
	 * Gets a set of tip nodes descended from the given node.
	 */
	public static Set<jebl.evolution.graphs.Node> getDescendantTips(RootedTree tree, jebl.evolution.graphs.Node node) {

		Set<jebl.evolution.graphs.Node> tipNodes = new LinkedHashSet<jebl.evolution.graphs.Node>();
		getDescendantTips(tree, node, tipNodes);
		return tipNodes;
	}

	/**
	 * Private recursive function used by getDescendantTips.
	 */
	private static void getDescendantTips(RootedTree tree, jebl.evolution.graphs.Node node, Set<jebl.evolution.graphs.Node> tipNodes) {

		for (jebl.evolution.graphs.Node child : tree.getChildren(node)) {
			if (tree.isExternal(child)) {
				tipNodes.add(child);
			} else {
				getDescendantTips(tree, child, tipNodes);
			}
		}
	}

	/**
	 * Gets the most recent common ancestor (MRCA) node of a set of tip nodes.
	 * @param tree the Tree
	 * @param tipNodes a set of tip nodes
	 * @return the Node of the MRCA
	 */
	public static jebl.evolution.graphs.Node getCommonAncestorNode(RootedTree tree, Set<jebl.evolution.graphs.Node> tipNodes) {

		if (tipNodes.size() == 0) {
			throw new IllegalArgumentException("No leaf nodes selected");
		}

		if (tipNodes.size() == 1) return tipNodes.iterator().next();

		jebl.evolution.graphs.Node[] mrca = new jebl.evolution.graphs.Node[] { null };
		getCommonAncestorNode(tree, tree.getRootNode(), tipNodes, mrca);

		return mrca[0];
	}

	/**
	 * Private recursive function used by getCommonAncestorNode.
	 */
	private static int getCommonAncestorNode(RootedTree tree, jebl.evolution.graphs.Node node,
	                                         Set<jebl.evolution.graphs.Node> tipNodes,
	                                         jebl.evolution.graphs.Node[] mrca) {


		int matches = 0;

		for (jebl.evolution.graphs.Node child : tree.getChildren(node)) {
			if (tree.isExternal(child)) {
				if (tipNodes.contains(child)) {
					matches ++;
				}
			} else {

				matches += getCommonAncestorNode(tree, child, tipNodes, mrca);

				if (mrca[0] != null) {
					return matches;
				}
			}
		}

		// If we haven't already found the MRCA, test this node
		if (matches == tipNodes.size()) {
			mrca[0] = node;
		}

		return matches;
	}

	/**
	 * Performs the a monophyly test on a set of tip nodes. The nodes are monophyletic
	 * if there is a node in the tree which subtends all the tips in the set (and
	 * only those tips).
	 * @param tree a tree object to perform test on
	 * @param tipNodes a set containing the tip node.
	 * @return boolean is monophyletic?
	 */
	public static boolean isMonophyletic(RootedTree tree, Set<jebl.evolution.graphs.Node> tipNodes) {

		if (tipNodes.size() == 0) {
			throw new IllegalArgumentException("No tip nodes selected");
		}

		if (tipNodes.size() == 1) {
			// A single selected leaf is always monophyletic
			return true;
		}

		if (tipNodes.size() == tree.getExternalNodes().size()) {
			// All leaf nodes are selected
			return true;
		}

		int[] matchCount = new int[] { 0 };
		int[] tipCount = new int[] { 0 };

		Boolean result = isMonophyletic(tree, tree.getRootNode(), tipNodes, matchCount, tipCount);

		if (result != null) return result;

		return false;
	}

	/**
	 * Private recursive function used by isMonophyletic.
	 */
	private static Boolean isMonophyletic(RootedTree tree, jebl.evolution.graphs.Node node,
	                                      Set<jebl.evolution.graphs.Node> tipNodes,
	                                      int[] matchCount, int[] tipCount) {

		int mc = 0;
		int tc = 0;

		for (jebl.evolution.graphs.Node child : tree.getChildren(node)) {
			if (tree.isExternal(child)) {
				if (tipNodes.contains(child)) {
					mc ++;
				}
				tc ++;
			} else {

				Boolean result = isMonophyletic(tree, child, tipNodes, matchCount, tipCount);

				if (result != null) {
					return result;
				}

				mc += matchCount[0];
				tc += tipCount[0];
			}
		}


		matchCount[0] = mc;
		tipCount[0] = tc;

		// If we haven't already found the MRCA, test this node
		if (mc == tc && tc == tipNodes.size()) {
			// monophyletic
			return Boolean.TRUE;
		}

		if (mc != 0 && mc != tc) {
			// not monophyletic
			return Boolean.FALSE;
		}

		// no result yet
		return null;
	}

	/**
	 * Recursive function for constructing a newick tree representation in the given buffer.
	 */
	public static String uniqueNewick(RootedTree tree, jebl.evolution.graphs.Node node) {
		if (tree.isExternal(node)) {
			return tree.getTaxon(node).getName();
		} else {
			StringBuffer buffer = new StringBuffer("(");

			List<String> subtrees = new ArrayList<String>();

			for (jebl.evolution.graphs.Node child : tree.getChildren(node)) {
				subtrees.add(uniqueNewick(tree, child));
			}
			Collections.sort(subtrees);

			for (int i = 0; i < subtrees.size(); i++) {
				buffer.append(subtrees.get(i));
				if (i < subtrees.size() - 1) {
					buffer.append(",");
				}
			}
			buffer.append(")");

			return buffer.toString();
		}
	}

	/**
	 * Compares 2 trees and returns true if they have the same topology.
	 */
	public static boolean equal(RootedTree tree1, RootedTree tree2) {

		return uniqueNewick(tree1, tree1.getRootNode()).equals(uniqueNewick(tree2, tree2.getRootNode()));
	}


}

