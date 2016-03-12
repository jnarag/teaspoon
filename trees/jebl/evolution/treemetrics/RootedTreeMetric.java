package jebl.evolution.treemetrics;

/**
 * @author rambaut
 *         Date: Jun 25, 2006
 *         Time: 12:12:34 AM
 */
public interface RootedTreeMetric {
	/**
	 * calculates the metric between two rooted trees
	 * @param tree1 first tree
	 * @param tree2 second tree
	 * @return the tree metric value
	 */
	double getMetric(jebl.evolution.trees.RootedTree tree1, jebl.evolution.trees.RootedTree tree2);
}
