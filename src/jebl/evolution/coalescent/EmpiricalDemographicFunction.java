/*
 * EmpiricalDemographicFunction.java
 *
 * (c) 2005 JEBL Development Team
 *
 * This package is distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package jebl.evolution.coalescent;

import jebl.evolution.coalescent.DemographicFunction;

/**
 * @author Oliver Pybus
 * @author Andrew Rambaut
 * @version $Id$
 */
public class EmpiricalDemographicFunction implements DemographicFunction {

    public EmpiricalDemographicFunction(double[] populationSizes, double[] times) {
        this.populationSizes = populationSizes;
        this.times = times;
    }

    /**
     * Gets the value of the demographic function N(t) at time t.
     */
    public double getDemographic(double t) {

        assert(t >= 0.0);

        // If time is beyond end of timeseries, then return last value of series
        if (t > times[times.length - 1]) {
            return populationSizes[populationSizes.length - 1];
        }

        for (int i = 0; i < times.length; i++) {
            if (times[i] == t) {
                return populationSizes[i];
            } else if (times[i] > t) {
                // Do linear interpolation. I think this works for both t[x]>t[x-1] and t[x]<t[x-1]
                double proportion = (t - times[i - 1]) / (times[i] - times[i - 1]);
                double popSize = populationSizes[i-1] + (proportion*(populationSizes[i] - populationSizes[i - 1]));

                return popSize;
            }
        }

        throw new RuntimeException("Error in jebl.evolution.treesimulation.EmpiricalDemographicFunction.getDemographic: went off the end of the array");
    }


    /**
     * Returns value of demographic intensity function at time t
     * (= integral 1/N(x) dx from 0 to t).
     */
    public double getIntensity(double t) {
        throw new UnsupportedOperationException("getIntensity is not implemented in jebl.evolution.treesimulation.EmpiricalDemographicFunction");
    }

    /**
     * Returns value of inverse demographic intensity function
     * (returns time, needed for simulation of coalescent intervals).
     */
    public double getInverseIntensity(double x) {
        throw new UnsupportedOperationException("getInverseIntensity is not implemented in jebl.evolution.treesimulation.EmpiricalDemographicFunction");
    }

    public boolean hasIntegral() {
        return false;
    }

    public double getIntegral(double start, double finish) {
        return 0;
    }

    /**
     * Returns the number of arguments for this function.
     */
    public int getArgumentCount() {
        return 0;
    }

    /**
     * Returns the name of the nth argument of this function.
     */
    public String getArgumentName(int n) {
        return null;
    }

    /**
     * Returns the value of the nth argument of this function.
     */
    public double getArgument(int n) {
        return 0;
    }

    /**
     * Sets the value of the nth argument of this function.
     */
    public void setArgument(int n, double value) {
    }

    /**
     * Returns the lower bound of the nth argument of this function.
     */
    public double getLowerBound(int n) {
        return 0;
    }

    /**
     * Returns the upper bound of the nth argument of this function.
     */
    public double getUpperBound(int n) {
        return 0;
    }

    private final double[] populationSizes;
    private final double[] times;
}
