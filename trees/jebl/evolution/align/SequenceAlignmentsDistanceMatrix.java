package jebl.evolution.align;

import jebl.evolution.distances.BasicDistanceMatrix;
import jebl.evolution.distances.CannotBuildDistanceMatrixException;
import jebl.evolution.taxa.Taxon;
import jebl.util.ProgressListener;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * Builds a distance matrix by performing a series of pairwise alignments between the
 * specified sequences (unlike the methods in jebl.evolution.distances, which
 * extract the pairwise distances from a multiple sequence alignment).
 *
 * @author Joseph Heled
 * @version $Id: SequenceAlignmentsDistanceMatrix.java 842 2007-11-12 22:33:48Z twobeers $
 *
 */
public class SequenceAlignmentsDistanceMatrix extends BasicDistanceMatrix {
    public SequenceAlignmentsDistanceMatrix(List<jebl.evolution.sequences.Sequence> seqs, PairwiseAligner aligner, ProgressListener progress)
            throws CannotBuildDistanceMatrixException
    {
        super(getTaxa(seqs), getDistances(seqs, aligner, progress));
    }

    static List<Taxon> getTaxa(List<jebl.evolution.sequences.Sequence> seqs) {
        List<Taxon> t = new ArrayList<Taxon>();
        for( jebl.evolution.sequences.Sequence s : seqs ) {
            t.add(s.getTaxon());
        }
        return t;
    }



    private static double[][] getDistances(List<jebl.evolution.sequences.Sequence> seqs, PairwiseAligner aligner, final ProgressListener progress) throws CannotBuildDistanceMatrixException {
        final int n = seqs.size();
        double [][] d = new double[n][n];
        boolean isProtein = seqs.get(0).getSequenceType().getCanonicalStateCount()> 4;

        CompoundAlignmentProgressListener compoundProgress = new CompoundAlignmentProgressListener(progress,(n * (n - 1)) / 2);

        for(int i = 0; i < n; ++i) {
            for(int j = i+1; j < n; ++j) {
                PairwiseAligner.Result result = aligner.doAlignment(seqs.get(i), seqs.get(j), compoundProgress.getMinorProgress());
                compoundProgress.incrementSectionsCompleted(1);
                if(compoundProgress.isCanceled()) return d;
                if(isProtein) {
                    d[i][j] = new jebl.evolution.distances.JukesCantorDistanceMatrix(result.alignment, null).getDistances()[0][1];
                } else {
                    d[i][j] = new jebl.evolution.distances.F84DistanceMatrix(result.alignment).getDistances()[0][1];

                }
                d[j][i] = d[i][j];
            }
        }
        return d;
    }
}