package teaspoon.adaptation;

import java.util.Calendar;
import java.util.Random;

public class TeaspoonRandomSampler {
	// Seeded static random number generator
	private static Random seededRandomNumberGenerator = new Random(Calendar.getInstance().getTimeInMillis() + Thread.currentThread().getId());
	// The values we will resample from
	double[] valuesToSampleFrom;

	/**
	 * No-arg constructor is deprecated
	 * @deprecated
	 */
	@Deprecated
	public TeaspoonRandomSampler() {}

	/**
	 * Initialise a new sampler with observations and a seed. All calls to this sampler instance share a seed.
	 * @param observations
	 */
	public TeaspoonRandomSampler(double[] observations){
		this.valuesToSampleFrom=observations;
	}


	public static double sampleGamma(double k, double theta) {
		boolean accept = false;
		if (k < 1) {
			// Weibull algorithm
			double c = (1 / k);
			double d = ((1 - k) * Math.pow(k, (k / (1 - k))));
			double u, v, z, e, x;
			do {
				u = seededRandomNumberGenerator.nextDouble();
				v = seededRandomNumberGenerator.nextDouble();
				z = -Math.log(u);
				e = -Math.log(v);
				x = Math.pow(z, c);
				if ((z + e) >= (d + x)) {
					accept = true;
				}
			} while (!accept);
			return (x * theta);
		} else {
			// Cheng's algorithm
			double b = (k - Math.log(4));
			double c = (k + Math.sqrt(2 * k - 1));
			double lam = Math.sqrt(2 * k - 1);
			double cheng = (1 + Math.log(4.5));
			double u, v, x, y, z, r;
			do {
				u = seededRandomNumberGenerator.nextDouble();
				v = seededRandomNumberGenerator.nextDouble();
				y = ((1 / lam) * Math.log(v / (1 - v)));
				x = (k * Math.exp(y));
				z = (u * v * v);
				r = (b + (c * y) - x);
				if ((r >= ((4.5 * z) - cheng)) ||
						(r >= Math.log(z))) {
					accept = true;
				}
			} while (!accept);
			return (x * theta);
		}
	}

	public double[] sampleDirichlet(){
		double total = 0;
		double[] D = new double[4];
		D[0] = sampleGamma(valuesToSampleFrom[0],1);total+=D[0];
		D[1] = sampleGamma(valuesToSampleFrom[1],1);total+=D[1];
		D[2] = sampleGamma(valuesToSampleFrom[2],1);total+=D[2];
		D[3] = sampleGamma(valuesToSampleFrom[3],1);total+=D[3];

		D[0]=D[0]/total;
		D[1]=D[1]/total;
		D[2]=D[2]/total;
		D[3]=D[3]/total;
		return D;  
	}

	public double sampleBeta() {

		double total = 0;
		double[] D = new double[2];

		D[0] = sampleGamma(valuesToSampleFrom[0],1);total+=D[0];
		D[1] = sampleGamma(valuesToSampleFrom[1],1);total+=D[1];

		D[0] = D[0]/total;
		D[1] = D[1]/total;

		return D[0];
	}
}