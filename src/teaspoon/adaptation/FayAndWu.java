package teaspoon.adaptation;

public class FayAndWu extends DiversityStats {
// fay and wu weights heterozygosity of derived variants as opposed to the ancestral variants
	public FayAndWu() {
		// default no-arg constructor
		throw new RuntimeException("ERROR: please input the intger matrix and the ancestral matrix or just integer matrix");
	}
	
	public FayAndWu(int[][] m,int[] a){
		super(m,a);
	}
	
	
	public double[] Estimators(){
		double[] counter = new double[integer_matrix.length];
		for(int j=0;j<integer_matrix[0].length;j++){
			double temp = 0;
			if(bad_sites_list[j] == false){
				temp = preprocess.num_of_base(integer_matrix, integer_ancestral[j], j);	// if not equal to ancestral
				temp = n-temp;
				if(temp != 0){
					counter[(int) temp-1]++;
				}
			}
		}
	
		double thetaw = 0; // number of segregating sites
		double part1 = 0;
		double part2 = 0;
		// loop goes from 1 -> n-1 where n is the sample size
		// counter goes from 0 -> n-2. which traverses all entries except the last - equvalent to n-1
		for(int k=1;k<integer_matrix.length;k++){
			part1 += counter[k-1];
			part2 += 1.0/k;
		}
		thetaw = part1*(1.0/part2);
		double thetak = 0; // number of pairwise differences
		for(int k=1;k<integer_matrix.length;k++){
			thetak += (2.0*counter[k-1]*k*(n-k))/(n*(n-1.0));
		}
		double thetah = 0; // weighted
		for(int k=1;k<integer_matrix.length;k++){
			thetah += (2.0*counter[k-1]*k*k)/(n*(n-1.0));
		}
		
		
		double[] est = {thetaw,thetak,thetah};
		return est;
	}
	
	public double[] FayWu(){
		double[] est = Estimators();
		double D = est[1] - est[0];
		double H = est[1] - est[2];
		double[] DH = {D,H};
		return DH;
	}
}












