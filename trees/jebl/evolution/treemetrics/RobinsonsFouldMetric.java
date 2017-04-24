package jebl.evolution.treemetrics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class RobinsonsFouldMetric implements jebl.evolution.treemetrics.RootedTreeMetric {

	public RobinsonsFouldMetric() {
		taxonMap = null;
	}

	public RobinsonsFouldMetric(List<jebl.evolution.taxa.Taxon> taxa) {
		taxonMap = new HashMap<jebl.evolution.taxa.Taxon, Integer>();
		for (int i = 0; i < taxa.size(); i++) {
			taxonMap.put(taxa.get(i), i);
		}
	}

	public double getMetric(jebl.evolution.trees.RootedTree tree1, jebl.evolution.trees.RootedTree tree2) {

		Map<jebl.evolution.taxa.Taxon, Integer> tm = taxonMap;

		if (tm == null) {
			List<jebl.evolution.taxa.Taxon> taxa = new ArrayList<jebl.evolution.taxa.Taxon>(tree1.getTaxa());

			tm = new HashMap<jebl.evolution.taxa.Taxon, Integer>();
			for (int i = 0; i < taxa.size(); i++) {
				tm.put(taxa.get(i), i);
			}
		}

		Set<String> clades1 = getClades(tm, tree1);
		Set<String> clades2 = getClades(tm, tree2);

		clades1.removeAll(clades2);

		return clades1.size();
	}

	private Set<String> getClades(Map<jebl.evolution.taxa.Taxon, Integer> taxa, jebl.evolution.trees.RootedTree tree) {

		Set<String> clades = new HashSet<String>();

		getTips(taxa, tree, tree.getRootNode(), clades);

		return clades;
	}

	private Set<Integer> getTips(Map<jebl.evolution.taxa.Taxon, Integer> taxa, jebl.evolution.trees.RootedTree tree, jebl.evolution.graphs.Node node, Set<String> clades) {

		Set<Integer> tips = new TreeSet<Integer>();

		if (tree.isExternal(node)) {
			tips.add(taxa.get(tree.getTaxon(node)));
		} else {
			jebl.evolution.graphs.Node child1 = tree.getChildren(node).get(0);
			Set<Integer> tips1 = getTips(taxa, tree, child1, clades);

			jebl.evolution.graphs.Node child2 = tree.getChildren(node).get(1);
			Set<Integer> tips2 = getTips(taxa, tree, child2, clades);

			tips.addAll(tips1);
			tips.addAll(tips2);

			clades.add(getCladeString(tips));
		}

		return tips;
	}

	private static String getCladeString(Set<Integer> tips) {
		Iterator<Integer> iter = tips.iterator();
		StringBuffer buffer = new StringBuffer();
		buffer.append(iter.next());
		while (iter.hasNext()) {
			buffer.append(",");
			buffer.append(iter.next());
		}
		return buffer.toString();
	}

	private final Map<jebl.evolution.taxa.Taxon, Integer> taxonMap;
}
