package jebl.evolution.io;

import jebl.evolution.sequences.BasicSequence;
import jebl.evolution.sequences.SequenceType;
import jebl.evolution.taxa.Taxon;

import java.io.EOFException;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;

/**
 * Class for importing PHYLIP sequential file format
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PhylipSequentialImporter.java 840 2007-11-09 04:52:39Z twobeers $
 */
public class PhylipSequentialImporter implements SequenceImporter {

    /**
     * Constructor
     */
    public PhylipSequentialImporter(Reader reader, SequenceType sequenceType, int maxNameLength) {
        helper = new jebl.evolution.io.ImportHelper(reader);

        this.sequenceType = sequenceType;
        this.maxNameLength = maxNameLength;
    }

    /**
     * importSequences.
     */
    public List<jebl.evolution.sequences.Sequence> importSequences() throws IOException, ImportException {

        List<jebl.evolution.sequences.Sequence> sequences = new ArrayList<jebl.evolution.sequences.Sequence>();

        try {

            int taxonCount = helper.readInteger();
            int siteCount = helper.readInteger();

            String firstSeq = null;

            for (int i = 0; i < taxonCount; i++) {
                StringBuilder name = new StringBuilder();

                char ch = helper.read();
                int n = 0;
                while (!Character.isWhitespace(ch) && (maxNameLength < 1 || n < maxNameLength)) {
                    name.append(ch);
                    ch = helper.read();
                    n++;
                }

                StringBuilder seq = new StringBuilder(siteCount);
                helper.readSequence(seq, sequenceType, "", siteCount, "-", "?", ".", firstSeq);

                if (firstSeq == null) {
                    firstSeq = seq.toString();
                }
                sequences.add(new BasicSequence(sequenceType, Taxon.getTaxon(name.toString()), seq.toString()));
            }
        } catch (EOFException e) {
        }

        return sequences;
    }

    private final jebl.evolution.io.ImportHelper helper;
    private final SequenceType sequenceType;
    private int maxNameLength = 10;
}
