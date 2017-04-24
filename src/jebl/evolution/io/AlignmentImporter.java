package jebl.evolution.io;

import java.io.IOException;
import java.util.List;

import jebl.evolution.alignments.Alignment;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: AlignmentImporter.java 185 2006-01-23 23:03:18Z rambaut $
 */
public interface AlignmentImporter {

	/**
	 * importAlignment.
	 */
	List<Alignment> importAlignments() throws IOException, ImportException;
}
