package jebl.evolution.align;

import jebl.evolution.alignments.Alignment;
import jebl.util.ProgressListener;

/**
 * @author Joseph Heled
 * @version $Id: PairwiseAligner.java 315 2006-05-03 02:13:54Z alexeidrummond $
 *
 */
public interface PairwiseAligner {

    public class Result {
        final public Alignment alignment;
        final public double score;

        Result(Alignment alignment, double score) {
            this.alignment = alignment;
            this.score = score;
        }
    }

    Result doAlignment(jebl.evolution.sequences.Sequence seq1, jebl.evolution.sequences.Sequence seq2, ProgressListener progress);

    double getScore(jebl.evolution.sequences.Sequence seq1, jebl.evolution.sequences.Sequence seq2);
}
