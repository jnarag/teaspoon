package teaspoon.adaptation;

import teaspoon.app.utils.TeaspoonMethods;

public class Site_estimation {
	public final int[][] integer_matrix;
	public final int[] integer_ancestral;
	public final int[] codon_ancestral;
	public  boolean[] bad_sites_list;
	private final static double Small_pos_num = 1.0e-250;
	private static double s;
	private final static double Large_neg_num = -1.0e250;
	private final static double integratoraccuracy = 1.0e-6;

	public final String[] AA =	{"K","N","K","N","T","T","T","T","R","S","R","S","I","I","M","I","Q","H","Q","H","P","P","P","P",
			"R","R","R","R","L","L","L","L","E","D","E","D","A","A","A","A","G","G","G","G","V","V","V","V",
			"X","Y","X","Y","S","S","S","S","X","C","W","C","L","F","L","F","?","-","?" };
//	constructors	
	public Site_estimation(){
		// default no-arg constructor
		throw new RuntimeException("please input the raw intger matrix and the ancestral matrix");
	}


	public Site_estimation(int[][] m, int[] a){
		// creates good matricies which contain no gaps or sequencing errors
		// creates integer and amino acid matricies
		TeaspoonMethods creatematrix = new TeaspoonMethods();
		integer_matrix = creatematrix.get_object(m, a, "base");
		integer_ancestral = creatematrix.get_ancestral_object(m,a,"base");
		codon_ancestral = creatematrix.get_ancestral_object(m, a, "codon");	
		bad_sites_list = creatematrix.InvalidSites(integer_matrix, integer_ancestral);		
	}
//	*********************************************************************************
//	Finds the site frequency in a given number of bins	
	TeaspoonMethods preprocess = new TeaspoonMethods();
//	***********************************************************************
//	***********************************************************************
//	Funtions to evaluate the integral and the coefficient	
//	Find coefficient for beta-binomial integration
	public double findcoefficient(float n, float k){
		return Math.floor(0.5 + Math.exp( preprocess.factorial(n+1.0f) - preprocess.factorial(k) - preprocess.factorial(n-k) ));
	}

//	evaluate the integral using the trapezoid rule	
	public static double findintegral(double lowbound,double highbound,double n, double k){
		double s,st,ost,os;
		ost = os = Large_neg_num;
		if(highbound <= lowbound){throw new RuntimeException("ERROR: upper bound is smaller than lower bound");}
		for(int j=1;j<=50;j++){
			st = trpd(j,lowbound,highbound,n,k);
			s = (4.0*st-ost)/3.0;
//			check accuracy
			if (Math.abs(s-os) < integratoraccuracy*Math.abs(os)) {return s;}
			os = s;
			ost = st;
		}
		throw new RuntimeException("ERROR: Too many iterations");
	}
//	trapezoid rule	
	public static double trpd(int qqq, double low, double high, double n, double x){
		double www,tnm,sum,del,funclow,funchigh,funcx;
		int it;
		if (qqq==1.0){
			funclow = integralfunc(low,n,x);
			if(funclow == Small_pos_num){throw new RuntimeException("ERROR: Infinite Integral?");}
			funchigh = integralfunc(high,n,x);
			if(funchigh == Small_pos_num){throw new RuntimeException("ERROR: Infinite Integral?");}
			s =  0.5 * (high-low) * (funclow + funchigh);
			return s;
		} else{
			it=1;
//			Bit-shift left is equivalent to a multiplication by 2
			for(int j=1;j<qqq-1;j++){it <<=1;}
			tnm = it; // number of points
			del = (high-low)/tnm;
			www = low + 0.5*del;
			sum =0.0;
			for (int j=1;j<=it;j++,www+=del){
				funcx = integralfunc(www,n,x);
				if(funcx == Small_pos_num){throw new RuntimeException("ERROR: Infinite Integral?");}
				sum += funcx;
			}
			s = 0.5*(s+(high-low)*sum/tnm);	//New s uses previous s value
			return s; // return value
		}


	}
//	function definition
	public static double integralfunc(double p, double n, double x){
		return (Math.pow(p,x)*Math.pow((1.0-p), (n-x)));
	}	
//	**********************************************************************		

//	the program
//	**********************************************************************	
	// estimates the site frequencies - P(u<pi<v|N,D) - for a given site
	public double estimation(double u,double v, int site){
//		______________________________ initial statements _________________________________________
		double notderived=0.0;
		double coeff; double integral = 0.0; double D = 0.0;
		double N = integer_matrix.length;
		// find frequencies of observed
//		______________________________ solve the formula ___________________________________________	
//		base frequencies
		// find number of non derived nucleotides
		notderived = preprocess.num_of_base(integer_matrix, integer_ancestral[site], site);
		D = N-notderived; 
//		______________________________ Numerator,deominator,integral________________________________
		Double Derived = new Double(D);		      // derived 
		Double Number = new Double(N);
		float floatD = Derived.floatValue();	  // float not derived
		float floatN = Number.floatValue();
		// numerator
		coeff = findcoefficient(floatN,floatD);
		// Integral
		integral = findintegral(u,v,N,Derived.doubleValue());

		Double fx = new Double(coeff*integral);	// evaluate f(x)
		return fx.doubleValue();
	}

//	finds the identity of a site as silent replacement or ambiguous

//	**********************************************************************	
//	Obtain site identity
//	To check and see if a site is silent or replacement, takes raw unchanged integer and ancertral seqeunces
//	***********************************************************************
	// finds identity as silent or replacement

//	***************************************************************************	
//	****************** pybus-bhatt method *******************************************	

	public int getcodonnumber(int pos1, int pos2, int pos3){
		if (pos1<5 && pos2<5 && pos3<5) { // normal codon
			return (16*(pos1-1)) + (4*(pos2-1)) + (1*(pos3-1));
		}
		else if (pos1==5 && pos2==5 && pos3==5) { // codon is NNN
			return 64;
		}
		else if (pos1==5 && pos2==5 && pos3==5) { // codon is ---
			return 65;
		}
		return 66; // codon is unrecognised
	}

	public double is_silent(int index,int codon, int site, int pos){
		double value = 1.0;
		int[] ansbases = preprocess.codonSplitter(codon_ancestral[codon]);
		int anscodonNumber = getcodonnumber(ansbases[0],ansbases[1],ansbases[2]);
		int derivedcodonNumber = 65;
		if(pos == 1){	
			derivedcodonNumber = getcodonnumber(integer_matrix[index][site],ansbases[1],ansbases[2]);
		} else if(pos == 2){
			derivedcodonNumber = getcodonnumber(ansbases[0],integer_matrix[index][site],ansbases[2]);
		} else if(pos == 3){
			derivedcodonNumber = getcodonnumber(ansbases[0],ansbases[1],integer_matrix[index][site]);
		}
		if(AA[anscodonNumber].equals(AA[derivedcodonNumber]) ){
			value= 1.0;
		} else{
			value = 0.0;
		}
		return value;
	}


	public boolean isSameBase(int site){
		int[] which = preprocess.which_bases(integer_matrix, site);
		int total = which[0] + which[1] + which[2] + which[3];
		if(total == 1){
			return true;
		} else {
			return false;
		}

	}

	// returns true if there are no derived states
	public boolean not_variable(int site){
		double num_of_ans_base_in_main = preprocess.num_of_base(integer_matrix, integer_ancestral[site], site);
		if(num_of_ans_base_in_main == (double) integer_matrix.length){		// if all bases in the main alignment are the same as the ancestral
			return true;		// there are no derived states at this site
		} else{
			return false;
		}
	}
	// returns true if site is the same as ancestral
	public boolean same_as_ans(int site){
		if(integer_matrix[0][site] == integer_ancestral[site]){
			return true;
		} else {
			return false;
		}
	}

	public double derived_and_silent(int site, int codon,int pos){
		int index = 0;
		for(int i=0;i<integer_matrix.length;i++){
			if(integer_matrix[i][site] != integer_ancestral[site]){	
				index = i;
				break;
			}
		}
		return is_silent(index,codon,site,pos);

	}

	public double[] find_identity(int site, int codon){
		double[] identity = new double[3];
		int[] ancestralbases = preprocess.codonSplitter(codon_ancestral[codon]);
		int codonNumber =  getcodonnumber(ancestralbases[0],ancestralbases[1],ancestralbases[2]);
		// pos 1
		double pos1_identity = 0.0;
		if (codonNumber==8 || codonNumber==10 || codonNumber==24 ||
				codonNumber==25 || codonNumber==26 || codonNumber==27 ||
				codonNumber==28 || codonNumber==29 || codonNumber==30 ||
				codonNumber==31 || codonNumber==60 || codonNumber==62 ||
				// new ones
				codonNumber==18 || codonNumber==32 || codonNumber==34) {
			pos1_identity = 0.5;
		}
		// pos 2
		double pos2_identity = 0.0;
		if (codonNumber==48 || codonNumber==50 || codonNumber==56) {
			pos2_identity = 0.5;
		}
		// pos 3		
		double pos3_identity = 1.0;
		if (codonNumber==14 || codonNumber==58) {
			pos3_identity = 0.0;
		} else if (codonNumber==0 || codonNumber==1 || codonNumber==2 || codonNumber==3 ||
				codonNumber==9 || codonNumber==11 || codonNumber==12 || codonNumber==13 ||
				codonNumber==15 || codonNumber==16 || codonNumber==17 || codonNumber==18 ||
				codonNumber==19 || codonNumber==32 || codonNumber==33 || codonNumber==34 ||
				codonNumber==35 || codonNumber==48 || codonNumber==49 || codonNumber==50 ||
				codonNumber==51 || codonNumber==56 || codonNumber==57 || codonNumber==59 ||
				codonNumber==61 || codonNumber==63) {
			pos3_identity = 0.5;
		}

		// checking ambiguity
		// pos1*******
		if(pos1_identity == 0.5){
			if(not_variable(site)) {
				if(same_as_ans(site)) {
					pos1_identity = 1.0;
				}
				pos1_identity = is_silent(0,codon,site,1);
			} else {
				pos1_identity = derived_and_silent(site,codon,1);
			}
		}
		// pos 2*******
		if(pos2_identity == 0.5){
			if(not_variable(site+1)) {
				if(same_as_ans(site+1)) {
					pos2_identity = 1.0;
				}
				pos2_identity = is_silent(0,codon,site+1,2);
			} else {
				pos2_identity = derived_and_silent(site+1,codon,2);
			}
		}
		// pos3********
		if(pos3_identity == 0.5){
			if(not_variable(site+2)) {
				if(same_as_ans(site+2)) {
					pos3_identity = 1.0;
				}
				pos3_identity = is_silent(0,codon,site+2,3);
			} else {
				pos3_identity = derived_and_silent(site+2,codon,3);
			}
		}

		identity[0] = pos1_identity;
		identity[1] = pos2_identity;
		identity[2] = pos3_identity;

		return identity;
	}

//	third model	
//	second method. gives same results but takes average for 3 and 4 site poly	
	public double[][] find_identityMK(int site, int codon){
		double[][] identity = new double[3][2];
		int[] ancestralbases = preprocess.codonSplitter(codon_ancestral[codon]);
		int codonNumber =  getcodonnumber(ancestralbases[0],ancestralbases[1],ancestralbases[2]);
		int silcount = 0;
		// pos 1
		for(int i=0;i<integer_matrix.length;i++){	// loop through all samples
			int derivedcodonNumber = 0;
			double pos1_identity = 0.0;	//position 1
			if (codonNumber==8 || codonNumber==10 || codonNumber==24 ||
					codonNumber==25 || codonNumber==26 || codonNumber==27 ||
					codonNumber==28 || codonNumber==29 || codonNumber==30 ||
					codonNumber==31 || codonNumber==60 || codonNumber==62 ||
					// new ones
					codonNumber==18 || codonNumber==32 || codonNumber==34) {
				pos1_identity = 0.5;
			}
			if(pos1_identity==0.5){
				derivedcodonNumber = getcodonnumber(integer_matrix[i][site],ancestralbases[1],ancestralbases[2]);
				if(AA[codonNumber].equals(AA[derivedcodonNumber]) ){
					pos1_identity=1.0;
				} else {
					pos1_identity=0.0;
				}
			}
			if(pos1_identity==1.0){
				silcount++;
			}
		}
		identity[0][0] = silcount/integer_matrix.length;
		identity[0][1] = 1.0 - identity[0][0]; 
		// pos 2
		silcount = 0;
		for(int i=0;i<integer_matrix.length;i++){	// loop through all samples
			int derivedcodonNumber = 0;
			double pos2_identity = 0.0;	//position 1


			if (codonNumber==48 || codonNumber==50 || codonNumber==56) {
				pos2_identity = 0.5;
			}
			if(pos2_identity==0.5){
				derivedcodonNumber = getcodonnumber(ancestralbases[0],integer_matrix[i][site+1],ancestralbases[2]);
				if(AA[codonNumber].equals(AA[derivedcodonNumber]) ){
					pos2_identity=1.0;
				} else {
					pos2_identity=0.0;
				}
			}
			if(pos2_identity==1.0){
				silcount++;
			}	
		}
		identity[1][0] = silcount/integer_matrix.length;
		identity[1][1] = 1.0 - identity[1][0];
		// pos 3	
		silcount = 0;
		for(int i=0;i<integer_matrix.length;i++){	// loop through all samples
			int derivedcodonNumber = 0;
			double pos3_identity = 1.0;	//position 1

			if (codonNumber==14 || codonNumber==58) {
				pos3_identity= 0.0;
			} else if (codonNumber==0 || codonNumber==1 || codonNumber==2 || codonNumber==3 ||
					codonNumber==9 || codonNumber==11 || codonNumber==12 || codonNumber==13 ||
					codonNumber==15 || codonNumber==16 || codonNumber==17 || codonNumber==18 ||
					codonNumber==19 || codonNumber==32 || codonNumber==33 || codonNumber==34 ||
					codonNumber==35 || codonNumber==48 || codonNumber==49 || codonNumber==50 ||
					codonNumber==51 || codonNumber==56 || codonNumber==57 || codonNumber==59 ||
					codonNumber==61 || codonNumber==63) {
				pos3_identity = 0.5;
			}
			if(pos3_identity==0.5){
				derivedcodonNumber = getcodonnumber(ancestralbases[0],ancestralbases[1],integer_matrix[i][site+2]);
				if(AA[codonNumber].equals(AA[derivedcodonNumber]) ){
					pos3_identity=1.0;
				} else {
					pos3_identity=0.0;
				}
			}
			if(pos3_identity==1.0){
				silcount++;
			}
		}
		identity[2][0] = silcount/integer_matrix.length;
		identity[2][1] = 1.0 - identity[2][0];
		return identity;

	}


//	***************************************************************************	
//	***************************************************************************	

//	finds site frequency for a range between u and v. for all sites

//	***********************************************************************
//	finds the site frequency for every site in a given range u-v	
	public double[] find_site_frequency(double u,double v){
		double EstimatesPos1 = 0.0;double EstimatesPos2 = 0.0; double EstimatesPos3 = 0.0;
		double[] temp = new double[3];
		double rho = 0.0;
		double sigma = 0.0;
		double[] finalans = new double[6];
		double total = 0.0;		

		for(int site=0,codon=0; site<integer_matrix[0].length-2;site=site+3,codon++){
			temp = find_identity(site,codon);
			if(bad_sites_list[site]==false){
				EstimatesPos1 = estimation(u,v,site);		// find site freq for pos1
				total += EstimatesPos1;
				sigma += EstimatesPos1*temp[0];
				rho += EstimatesPos1*(1.0-temp[0]);
			}
			if(bad_sites_list[site+1]==false){
				EstimatesPos2 = estimation(u,v,site+1);		// find site freq for pos1
				total += EstimatesPos2;
				sigma += EstimatesPos2*temp[1];
				rho += EstimatesPos2*(1.0-temp[1]);
			}
			if(bad_sites_list[site+2]==false){
				EstimatesPos3 = estimation(u,v,site+2);		// find site freq for pos1
				total += EstimatesPos3;
				sigma += EstimatesPos3*temp[2];
				rho +=  EstimatesPos3*(1.0-temp[2]);
			}			
		}

		finalans[0] = sigma;
		finalans[1] = rho;
		finalans[2] = total;
		finalans[3] = rho/sigma;
		return finalans;
	}

}
