package jebl.evolution.treemetrics;

import java.util.ArrayList;
import java.util.List;

/**
 * Billera tree distance - sum of change in branch lengths required to transform one tree to the second
 *
 * Note that this interface is not optimal for a large set where all pairs are required.
 * Creating TreeBiPartitionInfo's as a pre step is better unless memory is an issue.
 * 
 * @author Joseph Heled
 * @version $Id$
 */
public class BilleraMetric implements RootedTreeMetric {
    public double getMetric(jebl.evolution.trees.RootedTree tree1, jebl.evolution.trees.RootedTree tree2) {
        List<jebl.evolution.taxa.Taxon> taxa = new ArrayList<jebl.evolution.taxa.Taxon>(tree1.getTaxa());
        jebl.evolution.trees.TreeBiPartitionInfo p1 = new jebl.evolution.trees.TreeBiPartitionInfo(tree1, taxa);
        jebl.evolution.trees.TreeBiPartitionInfo p2 = new jebl.evolution.trees.TreeBiPartitionInfo(tree2, taxa);
        return jebl.evolution.trees.TreeBiPartitionInfo.distance(p1, p2, jebl.evolution.trees.TreeBiPartitionInfo.DistanceNorm.NORM1);
    }
}
