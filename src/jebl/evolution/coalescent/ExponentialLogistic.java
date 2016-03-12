// ConstLogistic.java
//
// (c) 2002-2004 BEAST Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package jebl.evolution.coalescent;

/**
 * This class models logistic growth from an initial exponential phase.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: ExponentialLogistic.java 937 2008-08-08 00:11:32Z rambaut $
 *
 */
public class ExponentialLogistic extends LogisticGrowth {

    /**
     * Construct demographic model with default settings
     */
    public ExponentialLogistic() {
        // empty constructor
    }

	/**
     * Construct demographic model with given settings
     * @param N0    Current population size
     * @param r     Logistic-phase growth rate
     * @param c     Logistic shape
     * @param time  Exponential-logistic phase transition time
     * @param r2    Exponential-phase growth rate
     */
	public ExponentialLogistic(double N0, double r, double c, double time, double r2) {

		super(N0, r, c);
        this.r2 = r2;
        this.time = time;
    }

	public double getR2() { return r2; }
	public void setR2(double r2) { this.r2 = r2; }

    public double getTime() { return time; }
    public void setTime(double time) { this.time = time; }

	// Implementation of abstract methods

	public double getDemographic(double t) {

        double transition_time = getTime();

        // size of the population under the logistic at transition_time
        if (t < transition_time) {
            return super.getDemographic(t);
        } else {
            double r2 = getR2();
            double N1 = super.getDemographic(transition_time);
            return N1 * Math.exp(-r2*(t - transition_time));
        }
	}

	public double getIntensity(double t) {
        throw new UnsupportedOperationException();
	}

	public double getInverseIntensity(double x) {
        throw new UnsupportedOperationException();
	}

    public boolean hasIntegral() {
        return false;
    }

	public double getIntegral(double start, double finish) {
        throw new UnsupportedOperationException();
	}

    @Override
    public int getArgumentCount() {
        return 4;
    }

    @Override
    public String getArgumentName(int n) {
        switch (n) {
            case 0: return "N0";
            case 1: return "r";
            case 2: return "c";
            case 3: return "r2";
            case 4: return "t";
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    @Override
    public double getArgument(int n) {
        switch (n) {
            case 0: return getN0();
            case 1: return getGrowthRate();
            case 2: return getShape();
            case 3: return getR2();
            case 4: return getTime();
        }
        throw new IllegalArgumentException("Argument " + n + " does not exist");
    }

    @Override
    public void setArgument(int n, double value) {
        switch (n) {
            case 0: setN0(value); break;
            case 1: setGrowthRate(value); break;
            case 2: setShape(value); break;
            case 3: setR2(value); break;
            case 4: setTime(value); break;
            default: throw new IllegalArgumentException("Argument " + n + " does not exist");

        }
    }

    @Override
    public double getLowerBound(int n) {
        return 0.0;
    }

    @Override
    public double getUpperBound(int n) {
        return Double.POSITIVE_INFINITY;
    }

	//
	// private stuff
	//

    // the exponential phase growth rate
    private double r2 = 0.1;
    
    // the transition time
    private double time = 1.0;
}