package teaspoon.adaptation;

import java.util.ArrayList;


public class Method_no_Outgroup {
	private static double[] cof = {76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};



//	**********************************************************************	
//	**********************************************************************	
//	**********************************************************************		
//	Processing routines
//	Count number of bases in either i or j of sequence integer matrix
//	takes character values of base	
	public double num_of_base(int[][] matrix, char base, int site){
		double count = 0;
		int int_base = 5;	// unidentified base
		if (base == 'A'){
			int_base = 1;
		} else if (base == 'C'){
			int_base = 2;
		} else if (base == 'G'){
			int_base = 3;
		} else if (base == 'T'){
			int_base = 4;
		} else {
			int_base = 5;
		}		
		for (int i=0; i< matrix.length; i++){
			if (matrix[i][site] == int_base){
				count=count+1;						// counter
			}
		}
		return count;
	}
//	overloaded to take integer values of base
	public double num_of_base(int[][] matrix, int base, int site){
		double count = 0;	
		for (int i=0; i< matrix.length; i++){
			if (matrix[i][site] == base){
				count=count+1;						// counter
			}
		}
		return count;

	}
//	Find which bases are present
//	returns integer array 1=true 0=false
	public int[] which_bases(int[][] matrix, int site){
		int[] whichbases = new int[4];
		for (int k=0; k< whichbases.length; k++){
			whichbases[k] = 0;
		}
		for(int i=0; i< matrix.length; i++){
			if(matrix[i][site] == 1){
				whichbases[0] = 1;
			} else if(matrix[i][site] == 2){
				whichbases[1] = 1;
			} else if(matrix[i][site] == 3){
				whichbases[2] = 1;
			} else if(matrix[i][site] == 4){
				whichbases[3] = 1;
			}
		}
		return whichbases;
	}

//	finds the gamma function of a real number
	public static double gammaln(float xx){
		double x,y,tmp,ser;
		int j;
		y=x=xx;
		tmp=x+5.5;
		tmp -= (x+0.5)*Math.log(tmp);
		ser=1.000000000190015;
		for (j=0;j <= 5;j++) {ser += cof[j]/++y;}
		return -tmp+ Math.log(2.5066282746310005*ser/x);
	}

//	Uses the definition of the gamma fucntion to find the factorial	
	public double factorial(float n){
		if(n<=1){return 0.0;}		
		double value = gammaln(n + 1.0f);
		return value;
	}




//	Find bad sites - Ones with gaps and N's - input integer objects
	public int[] find_bad_sites(int[][] integer_matrix){
		boolean flag = false;
//		int sumbase=0; 
		ArrayList<Integer> badsites = new ArrayList<Integer>();
		for (int i = 0; i< integer_matrix[0].length; i++){	// for sites
			// remove any sites with gaps or invalid characters
			if(num_of_base(integer_matrix,5,i)>0){
				flag=true;
			}
//			*************** IMPORTANT if we dont want 3 4 poly sites uncomment this ******************			
//			// counts number of bases
//			sumbase=0;
//			int[] base = which_bases(integer_matrix,i);
//			for(int x=0;x<base.length;x++){
//				sumbase = sumbase + base[x];	// count polymorphisms for pos 1
//			}
//
//			check for 3 4 polymorphism			
//			if(sumbase>2){		// if 3 or 4 polymorphic
//				flag = true;	// site is bad
//			}

//			if site flagged for any of the above reasons label as a bad site
			if(flag == true){
				Integer temp = new Integer(i);
				badsites.add(temp);					// add this site to a flagged list
			}
			flag = false;							// reset flag status
		}
		int[] badlist = new int[badsites.size()];					
		for(int i = 0; i < badsites.size(); i++) {
			Integer temp2 = badsites.get(i);	// cast flagged sites to string
			badlist[i] = temp2.intValue();
		}
		return badlist;
	}


//	boolean array showing if a site is bad or not	
	public boolean[] bad_sites_list(int[][] good_integer){
		int[] badsites = find_bad_sites(good_integer);
		boolean[] bad = new boolean[good_integer[0].length];
		for(int i=0;i<good_integer[0].length;i++){
			bad[i]=false;
			for(int x=0;x<badsites.length;x++){
				if(i==badsites[x]){
					bad[i]=true;
					break;
				}
			}
		}
		return bad;
	}
//  number of bad sites
	public int number_of_bad_sites(boolean[] bad_sites_list){
		int count=0;
		for(int i=0;i<bad_sites_list.length;i++){
			if(bad_sites_list[i]){
				count=count+1;
			}
		}
		return count;
	}


}

