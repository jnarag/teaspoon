package teaspoon.adaptation;

public class FuAndLi extends DiversityStats {
	
	public FuAndLi() {
		// default no-arg constructor
		throw new RuntimeException(
		"ERROR: please input the intger matrix and the ancestral matrix or just integer matrix");
	}

	public FuAndLi(int[][] m){
		super(m);
	}
	
	public FuAndLi(int[][] m,int[] a){
		super(m,a);
	}

	public double numberOfSingletons(){
		double numberSingletons=0;
		for(int site=0;site<integer_matrix[0].length;site++){
			double[] C = new double[4];
			if(bad_sites_list[site]==false){
				C[0] = preprocess.num_of_base(integer_matrix, 1, site);
				C[1] = preprocess.num_of_base(integer_matrix, 2, site);
				C[2] = preprocess.num_of_base(integer_matrix, 3, site);
				C[3] = preprocess.num_of_base(integer_matrix, 4, site);
				if(C[0]==1||C[1]==1||C[2]==1||C[3]==1){
					numberSingletons++;
				}
			}
		}
		return numberSingletons;
	}

	// explicit counts
	public double[] counts(){
		double S = fastSS();
		double K = numberOfPairwiseDifferences();
		double n = numberOfSingletons();
		double[] values =  {S,K,n};
		return values;
		
	}
	// watterson estimates
	public double[] wattersonEstimates(){
		double[] wattersons = new double[3];
		double Es = 0.0;
		for(int i=1;i< n;i++){
			Es += 1.0/(i);		// harmonic series estimated from coalescent
		}
		double thetaPi = numberOfPairwiseDifferences();	// unbiased pairwise differences estimator
		double thetaS = fastSS()/Es;	// unbiased segregating sites estimator
		double thetan = ((n-1.0)/(n))*numberOfSingletons();
		wattersons[0] = thetaS;
		wattersons[1] = thetaPi;
		wattersons[2] = thetan;
		return wattersons;	
	}
	
	public double[] FuLiTest(){
		double thetaEst[] = wattersonEstimates();
		double S = fastSS();
		double Dstar = 0.0; double Fstar=0.0;
		double an = 0.0; double bn = 0.0; double an1 = 0.0;
		for(double i=1;i<n;i++){
			an += 1.0/i;
			bn += 1.0/(i*i);
		}
		for(double i=1;i==n;i++){
			an1 += 1.0/i;
		}
		// Dstar terms **********
		double betaterm1 = 1.0/(Math.pow(an, 2.0)+bn);
		double betaterm2 = bn/Math.pow(an,2.0);
		double betaterm3 = (2.0/n)*(1.0+(1.0/an)-an+(an/n));
		double betaD = betaterm1*(betaterm2-betaterm3-(1.0/(n*n)));
		double alphaD = (1.0/an)*(((n-1.0)/n) - (1.0/an)) - betaD;
		// calcualtion of Dstar
		Dstar = (thetaEst[0]-thetaEst[2])/(Math.sqrt((alphaD*S)+(betaD*S*S)));
		// Fstar terms  ***********
		double betafterm1 = 1.0/((n*n)+bn);
		double betafterm2 = ((2.0*n*n*n) + ((n*n)*110.0) - (255.0*n) + 153.0)/(n*n*9.0 * (n-1.0));
		double betafterm3 = (2.0*(n-1.0)*an)/Math.pow(n,2.0);
		double betaF = betafterm1 * (betafterm2 + betafterm3 + ((8.0*bn)/n));
		double alphafterm1 = ((4.0*n*n) + (19.0*n) + 3.0 - ((12.0*(n+1.0))*an1))/(3.0*n*(n-1.0));
		double alphaF = ((1.0/an) * alphafterm1) - betaF;
		// calculation of Fstar
		Fstar = (thetaEst[1]-thetaEst[2])/(Math.sqrt((alphaF*S)+(betaF*S*S)));
		double []DF = {Dstar,Fstar};
		return DF;
	}

}














