package jebl.evolution.alignments;

import jebl.evolution.sequences.BasicSequence;
import jebl.evolution.sequences.SequenceType;
import jebl.evolution.taxa.Taxon;

import java.util.List;
import java.util.Set;

/**
 * Date: 17/01/2006
 * Time: 08:08:44
 *
 * @author Joseph Heled
 * @version $Id: ResampledAlignment.java 482 2006-10-25 06:30:57Z twobeers $
 *
 * Provide a re-sampled alignment. This means an alignment constructed by choosing a set of sites from
 * the source alignment and concataneting them. The set may be of any length and may contain duplications
 * (sampling with replacment).
 *
 * Due to Java restrictions on constructors, class is implemented using delegation.
 */

public class ResampledAlignment implements Alignment {
    protected BasicAlignment alignment;

    /**
     *  Setup resampled alignment.
     *
     * @param srcAlignment  sample sites from this alignment
     * @param siteIndices   Use this set to construct the resampled alignment
     */
    public void init(Alignment srcAlignment, int[] siteIndices) {
        final int nNewSites = siteIndices.length;
        final int nSeqs = srcAlignment.getSequences().size();

        // Work directly with states (fastest)
        jebl.evolution.sequences.State[][] newSeqsStates = new jebl.evolution.sequences.State[nSeqs][];

        for(int k = 0; k < nSeqs; ++k) {
            newSeqsStates[k] = new jebl.evolution.sequences.State[nNewSites];
        }

        final List<jebl.evolution.sequences.Sequence> seqs = srcAlignment.getSequenceList();
        for(int n = 0; n < nNewSites; ++n) {
            final int fromSite = siteIndices[n];
            for(int k = 0; k < nSeqs; ++k) {
                newSeqsStates[k][n] = seqs.get(k).getState(fromSite);
            }
        }

        jebl.evolution.sequences.Sequence[] newSeqs = new jebl.evolution.sequences.Sequence[nSeqs];
        for(int k = 0; k < nSeqs; ++k) {
            jebl.evolution.sequences.Sequence src = seqs.get(k);
            newSeqs[k] = new BasicSequence(src.getSequenceType(), src.getTaxon(), newSeqsStates[k]);
        }
        alignment = new BasicAlignment(newSeqs);
    }

    public List<jebl.evolution.sequences.Sequence> getSequenceList() {
        return alignment.getSequenceList();
    }

    public int getPatternCount() {
        return alignment.getPatternCount();
    }

    public int getPatternLength() {
        return alignment.getPatternLength();
    }

    public List<Pattern> getPatterns() {
        return alignment.getPatterns();
    }

    public List<Taxon> getTaxa() {
        return alignment.getTaxa();
    }

    public SequenceType getSequenceType() {
        return alignment.getSequenceType();
    }

    public int getSiteCount() {
        return alignment.getSiteCount();
    }

    public Set<jebl.evolution.sequences.Sequence> getSequences() {
        return alignment.getSequences();
    }

    public jebl.evolution.sequences.Sequence getSequence(Taxon taxon) {
        return alignment.getSequence(taxon);
    }
}
