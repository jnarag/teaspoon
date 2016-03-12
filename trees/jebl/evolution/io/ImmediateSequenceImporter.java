package jebl.evolution.io;

import jebl.util.ProgressListener;

import java.io.IOException;

/**
 *
 * A sequence importer sending the sequences back one by one, which makes it
 * possible to import larger documents if handled wisely on the other side.
 * 
 * @author Joseph Heled
 * @version $Id: ImmediateSequenceImporter.java 465 2006-10-04 04:24:20Z twobeers $
 *
 */
public interface ImmediateSequenceImporter {
    public interface Callback {
        void add(jebl.evolution.sequences.Sequence seq);
    }

    void importSequences(Callback callback, ProgressListener progressListener) throws IOException, ImportException;
}
