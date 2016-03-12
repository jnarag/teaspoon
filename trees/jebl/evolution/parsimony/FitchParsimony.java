/*
 * FitchParsimony.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package jebl.evolution.parsimony;

import jebl.evolution.alignments.Patterns;
import jebl.evolution.sequences.SequenceType;
import jebl.evolution.trees.Tree;
import jebl.evolution.trees.Utils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for reconstructing characters using Fitch parsimony. This is intended to be much faster
 * than the static methods in the utility "Parsimony" class.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: FitchParsimony.java 604 2007-01-04 20:22:42Z msuchard $
 */
public class FitchParsimony implements ParsimonyCriterion {

	private final SequenceType sequenceType;
	private final int stateCount;
	private final boolean gapsAreStates;

	private Map<jebl.evolution.graphs.Node, boolean[][]> stateSets = new HashMap<jebl.evolution.graphs.Node, boolean[][]>();
	private Map<jebl.evolution.graphs.Node, jebl.evolution.sequences.State[]> states = new HashMap<jebl.evolution.graphs.Node, jebl.evolution.sequences.State[]>();

//    private boolean[][] union;         // Must now be local to recursive function
//    private boolean[][] intersection;  // as nodes are not guaranteed to be called in post-order

	private jebl.evolution.trees.RootedTree tree = null;
	private final List<jebl.evolution.alignments.Pattern> patterns;
	private List<jebl.evolution.taxa.Taxon> taxa;

	private boolean hasCalculatedSteps = false;
	private boolean hasRecontructedStates = false;

	private final double[] siteScores;

	public FitchParsimony(List<jebl.evolution.alignments.Pattern> patterns, boolean gapsAreStates) {
		if (patterns == null || patterns.size() == 0) {
			throw new IllegalArgumentException("The patterns cannot be null or empty");
		}

		this.sequenceType = patterns.get(0).getSequenceType();
		this.gapsAreStates = gapsAreStates;
		this.taxa = patterns.get(0).getTaxa();

		if (gapsAreStates) {
			stateCount = sequenceType.getCanonicalStateCount() + 1;
		} else {
			stateCount = sequenceType.getCanonicalStateCount();

		}

		this.patterns = patterns;

		this.siteScores = new double[patterns.size()];
	}

	public FitchParsimony(Patterns patterns, boolean gapsAreStates) {
		this(patterns.getPatterns(), gapsAreStates);
	}

	/**
	 * Calculates the minimum number of siteScores for the parsimony reconstruction of a
	 * a set of character patterns on a tree. This only does the first pass of the
	 * Fitch algorithm so it does not store ancestral state reconstructions.
	 *
	 * @param tree a tree object to reconstruct the characters on
	 * @return number of parsimony siteScores
	 */
	public double[] getSiteScores(Tree tree) {

		if (tree == null) {
			throw new IllegalArgumentException("The tree cannot be null");
		}

		if (!(tree instanceof jebl.evolution.trees.RootedTree)) {
			throw new IllegalArgumentException("The tree must be an instance of rooted tree");
		}

		if (this.tree == null || this.tree != tree) {
			this.tree = (jebl.evolution.trees.RootedTree) tree;

			if (!Utils.isBinary(this.tree)) {
				throw new IllegalArgumentException("The Fitch algorithm can only reconstruct ancestral states on binary trees");
			}

			initialize();
		}

		if (!hasCalculatedSteps) {
			for (int i = 0; i < siteScores.length; i++) {
				siteScores[i] = 0;
			}
			calculateSteps(this.tree); //this.tree.getRootNode());
			hasCalculatedSteps = true;
		}


		return siteScores;
	}

	public double getScore(Tree tree) {

		getSiteScores(tree);

		double score = 0;

		int i = 0;
		for (jebl.evolution.alignments.Pattern pattern : patterns) {
			score += siteScores[i] * pattern.getWeight();
			i++;
		}
		return score;
	}

	/**
	 * Returns the reconstructed character states for a given node in the tree. If this method is repeatedly
	 * called with the same tree and patterns then only the first call will reconstruct the states and each
	 * subsequent call will return the stored states.
	 *
	 * @param tree a tree object to reconstruct the characters on
	 * @param node the node of the tree
	 * @return an array containing the reconstructed states for this node
	 */
	public jebl.evolution.sequences.State[] getStates(Tree tree, jebl.evolution.graphs.Node node) {

		getSiteScores(tree);

		if (!hasRecontructedStates) {
			reconstructStates(this.tree.getRootNode(), null);
			hasRecontructedStates = true;
		}

		return states.get(node);
	}

	private void initialize() {
		hasCalculatedSteps = false;
		hasRecontructedStates = false;

		for (jebl.evolution.graphs.Node node : tree.getNodes()) {
			boolean[][] stateSet = new boolean[patterns.size()][stateCount];
			stateSets.put(node, stateSet);

			jebl.evolution.sequences.State[] stateArray = new jebl.evolution.sequences.State[patterns.size()];
			states.put(node, stateArray);
		}
	}

	/**
	 * This is the first pass of the Fitch algorithm. This calculates the set of states
	 * at each node and counts the total number of siteScores (the score). If that is all that
	 * is required then the second pass is not necessary.
	 */
	private void calculateSteps(jebl.evolution.trees.RootedTree tree) {

		// nodes in pre-order
		final List<jebl.evolution.graphs.Node> nodes = Utils.getNodes(tree, tree.getRootNode());

		// used as locals in the loop below, allocated once
		boolean[][] union = new boolean[patterns.size()][stateCount];
		boolean[][] intersection = new boolean[patterns.size()][stateCount];

		// iterate in reverse - post order. State of child is gurantted to be reasy before parent

		for (int k = nodes.size() - 1; k >= 0; --k) {
			final jebl.evolution.graphs.Node node = nodes.get(k);
			final boolean[][] nodeStateSet = stateSets.get(node);

			if (tree.isExternal(node)) {
				boolean[][] stateSet = stateSets.get(node);
				jebl.evolution.sequences.State[] stateArray = states.get(node);

				for (int i = 0; i < patterns.size(); ++i) {
					jebl.evolution.alignments.Pattern pattern = patterns.get(i);

					jebl.evolution.taxa.Taxon taxon = tree.getTaxon(node);
					int index = taxa.indexOf(taxon);

					if (index == -1)
						throw new IllegalArgumentException("Unknown taxon, " + taxon.getName() + " in tree");

					jebl.evolution.sequences.State state = pattern.getState(index);
					stateArray[i] = state;
					if (gapsAreStates && state.isGap()) {
						stateSet[i][stateCount - 1] = true;
					} else {
						for (jebl.evolution.sequences.State canonicalState : state.getCanonicalStates()) {
							stateSet[i][canonicalState.getIndex()] = true;
						}
					}
				}
			} else {
				boolean first = true;
				for (jebl.evolution.graphs.Node child : tree.getChildren(node)) {
					boolean[][] childStateSet = stateSets.get(child);
					if (first) {
						for (int i = 0; i < patterns.size(); i++) {
							copyOf(childStateSet[i], union[i]);
							copyOf(childStateSet[i], intersection[i]);
						}
						first = false;
					} else {
						for (int i = 0; i < patterns.size(); i++) {
							unionOf(union[i], childStateSet[i], union[i]);
							intersectionOf(intersection[i], childStateSet[i], intersection[i]);
						}
					}
				}

				for (int i = 0; i < patterns.size(); i++) {
					if (sizeOf(intersection[i]) > 0) {
						copyOf(intersection[i], nodeStateSet[i]);
					} else {
						copyOf(union[i], nodeStateSet[i]);
						siteScores[i]++;
					}
				}
			}
		}
	}


	private String printState(boolean[][] stateSet) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0, n = stateSet.length; i < n; i++) {
			sb.append("site " + i);
			for (int j = 0, l = stateSet[i].length; j < l; j++) {
				sb.append(" " + (stateSet[i][j] ? "T" : "F"));
			}
			sb.append("\n");
		}
		return sb.toString();
	}


	private String printState(boolean[] stateSet) {
		StringBuffer sb = new StringBuffer();
//			for(int i=0,n=stateSet.length; i<n; i++) {
//		int i = 0;
//				sb.append("site "+i);
		for (int j = 0, l = stateSet.length; j < l; j++) {
			sb.append(" " + (stateSet[j] ? "T" : "F"));
		}
//				sb.append("\n");
//			}
		return sb.toString();
	}


	/**
	 * The second pass of the Fitch algorithm. This reconstructs the ancestral states at
	 * each node.
	 *
	 * @param node
	 * @param parentStates
	 */
	private void reconstructStates(jebl.evolution.graphs.Node node, jebl.evolution.sequences.State[] parentStates) {

		if (!tree.isExternal(node)) {
			boolean[][] nodeStateSet = stateSets.get(node);
			jebl.evolution.sequences.State[] nodeStates = states.get(node);

			for (int i = 0; i < patterns.size(); i++) {

				if (parentStates != null && nodeStateSet[i][parentStates[i].getIndex()]) {
					nodeStates[i] = parentStates[i];
				} else {
					int first = firstIndexOf(nodeStateSet[i]);
					nodeStates[i] = sequenceType.getState(first);
				}
			}

			for (jebl.evolution.graphs.Node child : tree.getChildren(node)) {
				reconstructStates(child, nodeStates);
			}
		}
	}

	private static void copyOf(boolean[] s, boolean[] d) {

		for (int i = 0; i < d.length; i++) {
			d[i] = s[i];
		}
	}

	private static void unionOf(boolean[] s1, boolean[] s2, boolean[] d) {

		for (int i = 0; i < d.length; i++) {
			d[i] = s1[i] || s2[i];
		}
	}

	private static void intersectionOf(boolean[] s1, boolean[] s2, boolean[] d) {

		for (int i = 0; i < d.length; i++) {
			d[i] = s1[i] && s2[i];
		}
	}

	private static int firstIndexOf(boolean[] s1) {

		for (int i = 0; i < s1.length; i++) {
			if (s1[i]) {
				return i;
			}
		}
		return -1;
	}

	private static int sizeOf(boolean[] s1) {

		int count = 0;
		for (int i = 0; i < s1.length; i++) {
			if (s1[i]) count += 1;
		}
		return count;
	}

}