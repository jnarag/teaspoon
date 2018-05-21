package teaspoon.adaptation;

import teaspoon.app.utils.TeaspoonMethods;

public class DiversityStats {

	public final int[][] integer_matrix;
	public final int[] integer_ancestral;
	public final boolean[] bad_sites_list;
	TeaspoonMethods preprocess = new TeaspoonMethods();
	MethodWithoutOutgroup preprocess2 = new MethodWithoutOutgroup();
	public final Double n;
	public final Double numsites;

	public DiversityStats() {
		// default no-arg constructor
		throw new RuntimeException(
		"ERROR: please input the intger matrix and the ancestral matrix or just integer matrix");
	}

	public DiversityStats(int[][] m) {
		// if no outgroup then creates just a bad sites list
		integer_matrix = m;
		integer_ancestral = null;
		bad_sites_list = preprocess2.bad_sites_list(m);	
		n = new Double(m.length);
		numsites = (double)m[0].length;
	}


	public DiversityStats(int[][] m, int[] a) {
		// creates good matricies which contain no gaps or sequencing errors
		// creates integer and amino acid matricies
		integer_matrix = preprocess.get_object(m, a, "base");
		integer_ancestral = preprocess.get_ancestral_object(m,a,"base");
		bad_sites_list = preprocess.bad_sites_list(m, a);	
		n = new Double(m.length);
		numsites = (double)m[0].length;
	}




	public double numberOfSegregatingSites(){
		double numberSegregatingSites = 0.0;
		for(int site=0;site<integer_matrix[0].length;site++){
			if(bad_sites_list[site]==false){
				int[] basesPresent = preprocess.which_bases(integer_matrix, site);	// check which bases are present
				int total = basesPresent[0] +basesPresent[1] +basesPresent[2] + basesPresent[3];	// check if polymorphic
				if(total != 1){		// if polymorphic
					numberSegregatingSites++;	// add to segregating sites
				}
			}
		}
		return numberSegregatingSites;
	}
	
	public double fastSS(){
		double numberSegregatingSites = 0.0;
		int[] ss = new int[integer_matrix[0].length];
		int l = integer_matrix[0].length;
		for(int i=1;i<integer_matrix.length;i++){
			for(int j=0;j<l;j++){
				if(ss[j]==0 && integer_matrix[i][j] != integer_matrix[0][j] && bad_sites_list[j]==false ){
					ss[j]=1;
					numberSegregatingSites++;
				}			
			}
		}
		return numberSegregatingSites;
	}
	
	


	public double theta(){
		double Es = 0.0;
		for(double i=1.0;i< n;i++){
			Es += 1.0/(i);		// harmonic series estimated from coalescent
		}
		double theta = fastSS()/Es;
//		double theta = numberOfPairwiseDifferences();
		return theta;  // theta for each individual site
	}
	
	public double thetaPerSite(){
		double Es = 0.0;
		for(double i=1.0;i< n;i++){
			Es += 1.0/(i);		// harmonic series estimated from coalescent
		}
		double theta = fastSS()/Es;
		double numsites=0;
		for(int i=0;i<bad_sites_list.length;i++){
			if(bad_sites_list[i]==false){
				numsites++;
			}
		}
		
//		double theta = numberOfPairwiseDifferences();
		return theta/numsites;  // theta for each individual site
	}

	public double numberOfPairwiseDifferences(){
		double numberPairwiseDifferences = 0;
		double n = integer_matrix.length;
		double[] n1 = new double[4];
		int[] types = {1,2,3,4};
		for(int site=0;site<integer_matrix[0].length;site++){
			if(bad_sites_list[site]==false){		// if site is good
				n1 = new double[4];
				int[] basesPresent = preprocess.which_bases(integer_matrix, site);	// check which bases are present
				int total = basesPresent[0] + basesPresent[1] + basesPresent[2] + basesPresent[3];	// check if polymorphic
				double numdif = 0.0;
				if(total==2){
					n1[0] = preprocess.num_of_base(integer_matrix,integer_matrix[0][site], site);
					numdif = n1[0]*(n-n1[0]);
				}
				if(total==3){
					int count=0;
					for(int i=0;i<basesPresent.length;i++){
						if(basesPresent[i]==1){
							n1[count] = preprocess.num_of_base(integer_matrix,types[i], site);
							count++;
						}
					}
					numdif = (n1[0]*(n-n1[0])) + (n1[1]*(n-n1[0]-n1[1]));					
				}
				if(total==4){
					int count=0;
					for(int i=0;i<basesPresent.length;i++){
						if(basesPresent[i]==1){
							n1[count] = preprocess.num_of_base(integer_matrix,types[i], site);
							count++;
						}
					}
					numdif = (n1[0]*(n-n1[0])) + (n1[1]*(n-n1[0]-n1[1])) + (n1[2]*(n-n1[0]-n1[1]-n1[2])) ; 
				}
				numberPairwiseDifferences += numdif;
			}
		}		
		// for avg pairwise differences divide by the total number of possible comparisons
		// which is n choose 2
		//	double totalComparisons = Math.exp(preprocess.factorial(n.floatValue()))/(2*(Math.exp(preprocess.factorial(n.floatValue()-2f))));
		double totalComparisons = (n*(n-1.0))/2.0;
		double ans = numberPairwiseDifferences/totalComparisons;
		return ans;
	}

	public double[] wattersonEstimates(){
		double[] wattersons = new double[2];
		double Es = 0.0;
		for(double i=1.0;i<n;i++){
			Es += 1.0/(i);		// harmonic series estimated from coalescent
		}
		double thetaPi = numberOfPairwiseDifferences();	// unbiased pairwise differences estimator
		double thetaS = numberOfSegregatingSites()/Es;	// unbiased segregating sites estimator
		wattersons[0] = thetaS;
		wattersons[1] = thetaPi;
		return wattersons;	
	}

	public double TajimasD(){
		double[] thetaEst = wattersonEstimates(); // index 0 is segregating sites, 1 pairwise differences
		double Es = 0.0;
		for(double i=1.0;i< n;i++){
			Es += 1.0/(i);		// harmonic series estimated from coalescent
		}
		double S = thetaEst[0]*Es;
		double a1 = 0.0; double a2 = 0.0;
		for(double i=1.0;i<n;i++){
			a1 += (1.0/i);
			a2 += (1.0/(Math.pow(i,2)));
		}
		double b1 = (n+1.0)/(3.0*(n-1.0));
		double b2 = (2.0*(Math.pow(n,2.0) + n + 3.0))/(9.0*n*(n-1.0));
		double c1 = b1 - (1.0/a1);
		double c2 = b2 - ((n+2.0)/(a1*n)) + (a2/Math.pow(a1, 2.0));
		double e1 = c1/a1;
		double e2 = c2/(Math.pow(a1,2.0) + a2);
		double D = (thetaEst[1] - thetaEst[0])/(Math.sqrt((e1*S)+(e2*S*(S-1.0))));
		return D;

	}

	public double[] TajimasDValues(){
		double[] thetaEst = wattersonEstimates(); // index 0 is segregating sites, 1 pairwise differences
		double Es = 0.0;
		for(double i=1.0;i< n;i++){
			Es += 1.0/(i);		// harmonic series estimated from coalescent
		}
		double S = thetaEst[0]*Es;
		double a1 = 0.0; double a2 = 0.0;
		for(double i=1.0;i<n;i++){
			a1 += (1.0/i);
			a2 += (1.0/(Math.pow(i,2)));
		}
		double b1 = (n+1.0)/(3.0*(n-1.0));
		double b2 = (2.0*(Math.pow(n,2.0) + n + 3.0))/(9.0*n*(n-1.0));
		double c1 = b1 - (1.0/a1);
		double c2 = b2 - ((n+2.0)/(a1*n)) + (a2/Math.pow(a1, 2.0));
		double e1 = c1/a1;
		double e2 = c2/(Math.pow(a1,2.0) + a2);
		double D = (thetaEst[1] - thetaEst[0])/(Math.sqrt((e1*S)+(e2*S*(S-1.0))));
		TajimasTable confidenceTable = new TajimasTable();
		double lowConfidence = 0;
		double highConfidence = 0;
		for(int i=0;i<confidenceTable.numOfSamples.length;i++){
			if(n==confidenceTable.numOfSamples[i]){
				lowConfidence = confidenceTable.lowConfidence[i];
				highConfidence = confidenceTable.highConfidence[i];
				break;
			} else{
				lowConfidence = highConfidence = Double.NaN;
			}	
		}
		double[] valuearray = new double[3];
		if(new Double(D).isNaN()){
			valuearray[0] = 0.0;	
		} else {
		valuearray[0]= D;
		}
		valuearray[1]= highConfidence;
		valuearray[2]= lowConfidence;
		return valuearray;	
	}

	public String[] headers(){
		String[] headers = new String[3];
		headers[0] = "D";
		headers[1] = "highConfidence";
		headers[2] = "lowConfidence";
		return headers;
	}

}























