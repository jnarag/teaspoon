package jebl.evolution.io;

import java.io.IOException;
import java.util.Collection;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: SequenceExporter.java 429 2006-08-26 18:17:39Z rambaut $
 */
public interface SequenceExporter {

    /**
     * exportSequences.
     */
    void exportSequences(Collection<? extends jebl.evolution.sequences.Sequence> sequences) throws IOException;
}
