package jebl.evolution.distances;

import jebl.evolution.taxa.Taxon;
import jebl.util.ProgressListener;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BasicDistanceMatrix.java 842 2007-11-12 22:33:48Z twobeers $
 */
public class BasicDistanceMatrix implements DistanceMatrix {

    public BasicDistanceMatrix(Collection<Taxon> taxa, double[][] distances) {

        if (distances == null || distances.length == 0) {
            throw new IllegalArgumentException("Source distance matrix is null or empty");
        }

        if (distances.length != distances[0].length) {
            throw new IllegalArgumentException("Source distance matrix is not square");
        }

        if (distances.length != taxa.size()) {
            throw new IllegalArgumentException("Source distance matrix dimensions do not match the number of taxa");
        }

        this.taxa = new ArrayList<Taxon>(taxa);
        this.distances = distances;
    }

    /**
     * Gets the size of the matrix (which is square), i.e., number of rows or columns.
     *
     * @return the size
     */
    public int getSize() {
        return distances.length;
    }

    /**
     * @return the list of taxa that the state values correspond to.
     */
    public List<Taxon> getTaxa() {
        return taxa;
    }

    /**
     * Gets the distance at a particular row and column
     *
     * @param row    the row index
     * @param column the column index
     * @return the distance
     */
    public double getDistance(int row, int column) {
        return distances[row][column];
    }

    /**
     * Gets the distance between 2 taxa
     *
     * @param taxonRow
     * @param taxonColumn
     * @return the distance
     */
    public double getDistance(Taxon taxonRow, Taxon taxonColumn) {
        int row = taxa.indexOf(taxonRow);
        if (row == -1) {
            throw new IllegalArgumentException("The row taxon, " + taxonRow.getName() + " is not found in this matrix");
        }

        int column = taxa.indexOf(taxonColumn);
        if (column == -1) {
            throw new IllegalArgumentException("The column taxon, " + taxonColumn.getName() + " is not found in this matrix");
        }

        return getDistance(row, column);
    }

    /**
     * Gets a sub-matrix for only those taxa in the collection (all
     * of which should be present in this matrix).
     *
     * @param taxonSubset
     * @return the new submatrix
     */
    public DistanceMatrix getSubmatrix(Collection<Taxon> taxonSubset) {
        double[][] newDistances = new double[taxonSubset.size()][taxonSubset.size()];
        int i = 0;
        for (Taxon taxonRow : taxonSubset) {

            int row = taxa.indexOf(taxonRow);
            if (row == -1) {
                throw new IllegalArgumentException("The taxon, " + taxonRow.getName() + " is not found in this matrix");
            }

            int j = 0;
            for (Taxon taxonColumn : taxonSubset) {
                int column = taxa.indexOf(taxonColumn);
                if (column == -1) {
                    throw new IllegalArgumentException("The taxon, " + taxonColumn.getName() + " is not found in this matrix");
                }

                newDistances[i][j] = getDistance(row, column);
            }
            i++;
        }
        return new BasicDistanceMatrix(taxonSubset, newDistances);
    }

    /**
     * Gets a 2-dimensional array containing the distances
     *
     * @return the distances
     */
    public double[][] getDistances() {
        return distances;
    }

    private final List<Taxon> taxa;
    private final double[][] distances;

    protected interface PairwiseDistanceCalculator {
        double calculatePairwiseDistance(int taxon1, int taxon2) throws CannotBuildDistanceMatrixException;
    }

    /**
     * Utility method for building a matrix of distances based on a calculator for each pairwise sequence distance
     * @throws CannotBuildDistanceMatrixException only if useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable is false
     */
    protected static double[][] buildDistancesMatrix(PairwiseDistanceCalculator pairwiseDistanceCalculator, int dimension, boolean useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable, ProgressListener progress)
            throws CannotBuildDistanceMatrixException
    {
        double[][] distances = new double[dimension][dimension];

        float tot = (dimension * (dimension - 1)) / 2;
        int done = 0;
        final double noDistance=-1;
        double maxDistance=-1;
        for(int i = 0; i < dimension; ++i) {
            for(int j = i+1; j < dimension; ++j) {
                try {
                    distances[i][j] = pairwiseDistanceCalculator.calculatePairwiseDistance(i, j);
                    maxDistance=Math.max(distances[i][j],maxDistance);
                } catch (CannotBuildDistanceMatrixException e) {
                    if (!useTwiceMaximumDistanceWhenPairwiseDistanceNotCalculatable) {
                        throw e;
                    }
                   distances[i][j]=noDistance;
                }
                distances[j][i] = distances[i][j];
                if( progress != null ) progress.setProgress( ++done / tot);
            }
        }
        if (maxDistance<0) {
            throw new CannotBuildDistanceMatrixException("It is not possible to compute the Tamura-Nei genetic distance " +
                    "for these sequences because no pair of sequences overlap in the alignment.");
        }
        for(int i = 0; i < dimension; ++i) {
            for(int j = i+1; j < dimension; ++j) {
                if (distances[i][j]==noDistance) {
                    distances[i][j]= distances[j][i] =maxDistance*2;
                }
            }
        }
        return distances;
    }
}