package jebl.evolution.io;

import java.io.IOException;
import java.util.List;

import jebl.evolution.distances.DistanceMatrix;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DistanceMatrixImporter.java 185 2006-01-23 23:03:18Z rambaut $
 */
public interface DistanceMatrixImporter {

    enum Triangle { LOWER, UPPER, BOTH };

    /**
     * importDistances.
     */
    List<DistanceMatrix> importDistanceMatrices() throws IOException, ImportException;
}
