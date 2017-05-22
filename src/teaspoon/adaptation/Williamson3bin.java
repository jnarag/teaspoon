package teaspoon.adaptation;

import java.util.Random;
public class Williamson3bin extends SiteEstMulti{
	Random generator = new Random();
	double low_A;
	double mid_A;
	double high_A;
	double fix_A;
	double low_R;
	double low_S;
	double mid_R;
	double mid_S;
	double high_R;
	double high_S;
	double fix_R;
	double fix_S;
	double Nr;
	double Adapt;
	double Prop;
	double L;
	double[] gapcount;
	double[] gapSitecount;

	public Williamson3bin() {
		// default no-arg constructor
		throw new RuntimeException("ERROR: please input the intger matrix and the ancestral matrix");
	}

	public Williamson3bin(int[][] m, int[] a){
		super(m,a);
	}

	public Williamson3bin(int[][] m, int[] a,boolean[] badsites){
		super(m,a,badsites);
	}


	public double[][] williamsonEstimator(double[] identity,int sitelocation, double[] low, double[] mid, double[] high){
		double fs = 0;double fn = 0; double as = 0;double an = 0; double ds = 0; double dn = 0;double ns=0;double nn=0;
		SiteInfo Info = SiteInformation(sitelocation); 

		//		CASE 1
		if(Info.Case==1){ // invariant
			// do nothing as site has no information
		}
		//		CASE 2
		if(Info.Case==2){// fixed
			fs += 1.0*identity[0] ;
			fn += 1.0*identity[1] ;
		}
		//		CASE 3
		if(Info.Case==3){// 1 state derived and ans
			double l=0;
			double h=0;
			double m=0;
			double y=0; // number derived nucleotides
			for(int i=0;i<4;i++){
				if(Info.data[i].NObs>high[0] && Info.data[i].NObs<high[1] && Info.data[i].inans==false){
					h++;
					y++;
				}
				else if(Info.data[i].NObs>mid[0] && Info.data[i].NObs<=mid[1] && Info.data[i].inans==false){
					m++;
					y++;
				}
				else if(Info.data[i].NObs>low[0] && Info.data[i].NObs<=low[1] && Info.data[i].inans==false){
					l++;
					y++;
				}
			}
			double propfixed = 0.0;
			double proppoly = 1.0;
			if(y==0){
				System.out.println("Debug Problem variant sites");
				System.out.println(sitelocation);
			}
			// fixed catagory		
			fs = propfixed*identity[0] ;
			fn = propfixed*identity[1] ;
			//polymorphic catagory
			h=h/y;					
			l=l/y;
			m=m/y;			
			as = proppoly*h*identity[0];
			an = proppoly*h*identity[1];
			ns = proppoly*m*identity[0];
			nn = proppoly*m*identity[1];
			ds = proppoly*l*identity[0];
			dn = proppoly*l*identity[1];
		}
		// CASE 4
		if(Info.Case==4){// 2 state derived no ans
			double l=0;
			double h=0;
			double m=0;
			double y=0; // number derived nucleotides
			for(int i=0;i<4;i++){
				if(Info.data[i].NObs>high[0] && Info.data[i].NObs<high[1] && Info.data[i].inans==false){
					h++;
					y++;
				}
				else if(Info.data[i].NObs>mid[0] && Info.data[i].NObs<=mid[1] && Info.data[i].inans==false){
					m++;
					y++;
				}
				else if(Info.data[i].NObs>low[0] && Info.data[i].NObs<=low[1] && Info.data[i].inans==false){
					l++;
					y++;
				}
			}
			double propfixed = 0.5;
			double proppoly = 0.5;
			// fixed catagory		
			fs = propfixed*identity[0] ;
			fn = propfixed*identity[1] ;
			//polymorphic catagory
			h=h/y;					
			l=l/y;
			m=m/y;			
			as = proppoly*h*identity[0];
			an = proppoly*h*identity[1];
			ns = proppoly*m*identity[0];
			nn = proppoly*m*identity[1];
			ds = proppoly*l*identity[0];
			dn = proppoly*l*identity[1];

		}
		//		CASE 5
		if(Info.Case==5){// 2 state derived and ans
			double l=0;
			double h=0;
			double m=0;
			double y=0; // number derived nucleotides
			for(int i=0;i<4;i++){
				if(Info.data[i].NObs>high[0] && Info.data[i].NObs<high[1] && Info.data[i].inans==false){
					h++;
					y++;
				}
				else if(Info.data[i].NObs>mid[0] && Info.data[i].NObs<=mid[1] && Info.data[i].inans==false){
					m++;
					y++;
				}
				else if(Info.data[i].NObs>low[0] && Info.data[i].NObs<=low[1] && Info.data[i].inans==false){
					l++;
					y++;
				}
			}
			double propfixed = 0.0;
			double proppoly = 1.0;
			// fixed catagory		
			fs = propfixed*identity[0] ;
			fn = propfixed*identity[1] ;
			//polymorphic catagory
			h=h/y;					
			l=l/y;
			m=m/y;			
			as = proppoly*h*identity[0];
			an = proppoly*h*identity[1];
			ns = proppoly*m*identity[0];
			nn = proppoly*m*identity[1];
			ds = proppoly*l*identity[0];
			dn = proppoly*l*identity[1];
		}
		//		CASE 6
		if(Info.Case==6){// 3 state derived no ans
			double l=0;
			double h=0;
			double m=0;
			double y=0; // number derived nucleotides
			for(int i=0;i<4;i++){
				if(Info.data[i].NObs>high[0] && Info.data[i].NObs<high[1] && Info.data[i].inans==false){
					h++;
					y++;
				}
				else if(Info.data[i].NObs>mid[0] && Info.data[i].NObs<=mid[1] && Info.data[i].inans==false){
					m++;
					y++;
				}
				else if(Info.data[i].NObs>low[0] && Info.data[i].NObs<=low[1] && Info.data[i].inans==false){
					l++;
					y++;
				}
			}
			double propfixed = (1.0/3.0);
			double proppoly = (2.0/3.0);
			// fixed catagory		
			fs = propfixed*identity[0] ;
			fn = propfixed*identity[1] ;
			//polymorphic catagory
			h=h/y;					
			l=l/y;
			m=m/y;			
			as = proppoly*h*identity[0];
			an = proppoly*h*identity[1];
			ns = proppoly*m*identity[0];
			nn = proppoly*m*identity[1];
			ds = proppoly*l*identity[0];
			dn = proppoly*l*identity[1];
		}
		//		CASE 7
		if(Info.Case==7){// 3 state derived and ans
			double l=0;
			double h=0;
			double m=0;
			double y=0; // number derived nucleotides
			for(int i=0;i<4;i++){
				if(Info.data[i].NObs>high[0] && Info.data[i].NObs<high[1] && Info.data[i].inans==false){
					h++;
					y++;
				}
				else if(Info.data[i].NObs>mid[0] && Info.data[i].NObs<=mid[1] && Info.data[i].inans==false){
					m++;
					y++;
				}
				else if(Info.data[i].NObs>low[0] && Info.data[i].NObs<=low[1] && Info.data[i].inans==false){
					l++;
					y++;
				}
			}
			double propfixed = 0.0;
			double proppoly = 1.0;
			// fixed catagory		
			fs = propfixed*identity[0] ;
			fn = propfixed*identity[1] ;
			//polymorphic catagory
			h=h/y;					
			l=l/y;
			m=m/y;			
			as = proppoly*h*identity[0];
			an = proppoly*h*identity[1];
			ns = proppoly*m*identity[0];
			nn = proppoly*m*identity[1];
			ds = proppoly*l*identity[0];
			dn = proppoly*l*identity[1];
		}
		double[] store1 = {dn, nn, an, fn};
		double[] store2 = {ds, ns, as, fs};
		double[][] value = new double[2][store1.length];
		for(int i=0;i<store1.length;i++){
			value[0][i] = store1[i];
			value[1][i] = store2[i];
		}
		return value;
	}

	public double[][] williamsonEstimator_IncludeMainGaps(double[] identity,int sitelocation, double[] low, double[] mid, double[] high){
		double fs = 0;double fn = 0; double as = 0;double an = 0; double ds = 0; double dn = 0;double ns=0;double nn=0;
		SiteInfo Info = SiteInformation_IncludeMainGaps(sitelocation);
		//		CASE 1
		if(Info.Case==1){ // invariant
			// do nothing as site has no information
		}
		//		CASE 2
		if(Info.Case==2){// fixed
			fs += 1.0*identity[0] ;
			fn += 1.0*identity[1] ;
		}
		//		CASE 3
		if(Info.Case==3){// 1 state derived and ans
			double l=0;
			double h=0;
			double m=0;
			double y=0; // number derived nucleotides
			for(int i=0;i<4;i++){
				if(Info.data[i].NObs>high[0] && Info.data[i].NObs<high[1] && Info.data[i].inans==false){
					h++;
					y++;
				}
				else if(Info.data[i].NObs>mid[0] && Info.data[i].NObs<=mid[1] && Info.data[i].inans==false){
					m++;
					y++;
				}
				else if(Info.data[i].NObs>low[0] && Info.data[i].NObs<=low[1] && Info.data[i].inans==false){
					l++;
					y++;
				}
			}
			double propfixed = 0.0;
			double proppoly = 1.0;
			if(y==0){
				System.out.println("Debug Problem variant sites");
				System.out.println(sitelocation);
			}
			// fixed catagory		
			fs = propfixed*identity[0] ;
			fn = propfixed*identity[1] ;
			//polymorphic catagory
			h=h/y;					
			l=l/y;
			m=m/y;			
			as = proppoly*h*identity[0];
			an = proppoly*h*identity[1];
			ns = proppoly*m*identity[0];
			nn = proppoly*m*identity[1];
			ds = proppoly*l*identity[0];
			dn = proppoly*l*identity[1];
		}
		// CASE 4
		if(Info.Case==4){// 2 state derived no ans
			double l=0;
			double h=0;
			double m=0;
			double y=0; // number derived nucleotides
			for(int i=0;i<4;i++){
				if(Info.data[i].NObs>high[0] && Info.data[i].NObs<high[1] && Info.data[i].inans==false){
					h++;
					y++;
				}
				else if(Info.data[i].NObs>mid[0] && Info.data[i].NObs<=mid[1] && Info.data[i].inans==false){
					m++;
					y++;
				}
				else if(Info.data[i].NObs>low[0] && Info.data[i].NObs<=low[1] && Info.data[i].inans==false){
					l++;
					y++;
				}
			}
			double propfixed = 0.5;
			double proppoly = 0.5;
			// fixed catagory		
			fs = propfixed*identity[0] ;
			fn = propfixed*identity[1] ;
			//polymorphic catagory
			h=h/y;					
			l=l/y;
			m=m/y;			
			as = proppoly*h*identity[0];
			an = proppoly*h*identity[1];
			ns = proppoly*m*identity[0];
			nn = proppoly*m*identity[1];
			ds = proppoly*l*identity[0];
			dn = proppoly*l*identity[1];

		}
		//		CASE 5
		if(Info.Case==5){// 2 state derived and ans
			double l=0;
			double h=0;
			double m=0;
			double y=0; // number derived nucleotides
			for(int i=0;i<4;i++){
				if(Info.data[i].NObs>high[0] && Info.data[i].NObs<high[1] && Info.data[i].inans==false){
					h++;
					y++;
				}
				else if(Info.data[i].NObs>mid[0] && Info.data[i].NObs<=mid[1] && Info.data[i].inans==false){
					m++;
					y++;
				}
				else if(Info.data[i].NObs>low[0] && Info.data[i].NObs<=low[1] && Info.data[i].inans==false){
					l++;
					y++;
				}
			}
			double propfixed = 0.0;
			double proppoly = 1.0;
			// fixed catagory		
			fs = propfixed*identity[0] ;
			fn = propfixed*identity[1] ;
			//polymorphic catagory
			h=h/y;					
			l=l/y;
			m=m/y;			
			as = proppoly*h*identity[0];
			an = proppoly*h*identity[1];
			ns = proppoly*m*identity[0];
			nn = proppoly*m*identity[1];
			ds = proppoly*l*identity[0];
			dn = proppoly*l*identity[1];
		}
		//		CASE 6
		if(Info.Case==6){// 3 state derived no ans
			double l=0;
			double h=0;
			double m=0;
			double y=0; // number derived nucleotides
			for(int i=0;i<4;i++){
				if(Info.data[i].NObs>high[0] && Info.data[i].NObs<high[1] && Info.data[i].inans==false){
					h++;
					y++;
				}
				else if(Info.data[i].NObs>mid[0] && Info.data[i].NObs<=mid[1] && Info.data[i].inans==false){
					m++;
					y++;
				}
				else if(Info.data[i].NObs>low[0] && Info.data[i].NObs<=low[1] && Info.data[i].inans==false){
					l++;
					y++;
				}
			}
			double propfixed = (1.0/3.0);
			double proppoly = (2.0/3.0);
			// fixed catagory		
			fs = propfixed*identity[0] ;
			fn = propfixed*identity[1] ;
			//polymorphic catagory
			h=h/y;					
			l=l/y;
			m=m/y;			
			as = proppoly*h*identity[0];
			an = proppoly*h*identity[1];
			ns = proppoly*m*identity[0];
			nn = proppoly*m*identity[1];
			ds = proppoly*l*identity[0];
			dn = proppoly*l*identity[1];
		}
		//		CASE 7
		if(Info.Case==7){// 3 state derived and ans
			double l=0;
			double h=0;
			double m=0;
			double y=0; // number derived nucleotides
			for(int i=0;i<4;i++){
				if(Info.data[i].NObs>high[0] && Info.data[i].NObs<high[1] && Info.data[i].inans==false){
					h++;
					y++;
				}
				else if(Info.data[i].NObs>mid[0] && Info.data[i].NObs<=mid[1] && Info.data[i].inans==false){
					m++;
					y++;
				}
				else if(Info.data[i].NObs>low[0] && Info.data[i].NObs<=low[1] && Info.data[i].inans==false){
					l++;
					y++;
				}
			}
			double propfixed = 0.0;
			double proppoly = 1.0;
			// fixed catagory		
			fs = propfixed*identity[0] ;
			fn = propfixed*identity[1] ;
			//polymorphic catagory
			h=h/y;					
			l=l/y;
			m=m/y;			
			as = proppoly*h*identity[0];
			an = proppoly*h*identity[1];
			ns = proppoly*m*identity[0];
			nn = proppoly*m*identity[1];
			ds = proppoly*l*identity[0];
			dn = proppoly*l*identity[1];
		}
		double[] store1 = {dn, nn, an, fn};
		double[] store2 = {ds, ns, as, fs};
		double[][] value = new double[2][store1.length];
		for(int i=0;i<store1.length;i++){
			value[0][i] = store1[i];
			value[1][i] = store2[i];
		}
		return value;
	}

	public void williamson3bin_method(double neutralratio,double[] low, double[] mid, double[] high){
		double fs = 0;double fn = 0; double as = 0;double an = 0; double ds = 0; double dn = 0;double ns=0;double nn=0;
		for (int site = 0, codon = 0; site < integer_matrix[0].length - 2; site = site + 3, codon++) {
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { // check codon is not bad
				//		double[][] identity = find_identityMK(site, codon);
				double[][] identity = NGmethodWilliamson(site, codon);
				/************************	Position 1	************************************************/
				double[] pos1id = new double[2];
				pos1id[0] = identity[0][0];
				pos1id[1] = identity[0][1];
				double[][] pos1 = williamsonEstimator(pos1id,site,low,mid,high);
				dn+=pos1[0][0];nn+=pos1[0][1];an+=pos1[0][2];fn+=pos1[0][3];
				ds+=pos1[1][0];ns+=pos1[1][1];as+=pos1[1][2];fs+=pos1[1][3];

				/************************	Position 2	************************************************/
				double[] pos2id = new double[2];
				pos2id[0] = identity[1][0];
				pos2id[1] = identity[1][1];
				double[][] pos2 = williamsonEstimator(pos2id, site+1,low,mid,high);
				dn+=pos2[0][0];nn+=pos2[0][1];an+=pos2[0][2];fn+=pos2[0][3];
				ds+=pos2[1][0];ns+=pos2[1][1];as+=pos2[1][2];fs+=pos2[1][3];
				/************************	Position 3	************************************************/
				double[] pos3id = new double[2];
				pos3id[0] = identity[2][0];
				pos3id[1] = identity[2][1];
				double[][] pos3 = williamsonEstimator(pos3id,site+2,low,mid,high);
				dn+=pos3[0][0];nn+=pos3[0][1];an+=pos3[0][2];fn+=pos3[0][3];
				ds+=pos3[1][0];ns+=pos3[1][1];as+=pos3[1][2];fs+=pos3[1][3];

			}
		}
		double[] finalmat = new double[4];
		finalmat[0]= fn - ((fs)*(neutralratio)) ;  //adaptive fixations
		finalmat[1]= an - ((as)*(neutralratio)) ;	//adaptive substitutions
		finalmat[2]= dn - ((ds)*(neutralratio));	//adaptive substitutions
		finalmat[3]= nn - ((ns)*(neutralratio));

		

		for(int i=0;i<finalmat.length;i++){
			if(finalmat[i]<0){
				finalmat[i]=0;
			}
		}
		double totalAdapt = finalmat[0]+finalmat[1];		
		double propAdapt = 1.0 - ((as+fs)/(an+fn))*neutralratio;
		this.low_A = finalmat[2];
		this.mid_A = finalmat[3];
		this.high_A = finalmat[1];
		this.fix_A = finalmat[0];
		this.low_R=dn;
		this.low_S=ds;
		this.mid_R =nn;
		this.mid_S =ns;
		this.high_R =an;
		this.high_S =as;
		this.fix_R =fn;
		this.fix_S =fs;
		this.Adapt=totalAdapt;
		this.Prop=propAdapt;
		this.L=integer_matrix.length;
	}
	// 
	public void williamson3bin_method_delta(double neutralratio,double mu_ns,double[] low, double[] mid, double[] high){
		double fs = 0;double fn = 0; double as = 0;double an = 0; double ds = 0; double dn = 0;double ns=0;double nn=0;
		for (int site = 0, codon = 0; site < integer_matrix[0].length - 2; site = site + 3, codon++) {
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { // check codon is not bad
				//		double[][] identity = find_identityMK(site, codon);
				double[][] identity = NGmethodWilliamson(site, codon);
				/************************	Position 1	************************************************/
				double[] pos1id = new double[2];
				pos1id[0] = identity[0][0];
				pos1id[1] = identity[0][1];
				double[][] pos1 = williamsonEstimator(pos1id,site,low,mid,high);
				dn+=pos1[0][0];nn+=pos1[0][1];an+=pos1[0][2];fn+=pos1[0][3];
				ds+=pos1[1][0];ns+=pos1[1][1];as+=pos1[1][2];fs+=pos1[1][3];

				/************************	Position 2	************************************************/
				double[] pos2id = new double[2];
				pos2id[0] = identity[1][0];
				pos2id[1] = identity[1][1];
				double[][] pos2 = williamsonEstimator(pos2id, site+1,low,mid,high);
				dn+=pos2[0][0];nn+=pos2[0][1];an+=pos2[0][2];fn+=pos2[0][3];
				ds+=pos2[1][0];ns+=pos2[1][1];as+=pos2[1][2];fs+=pos2[1][3];
				/************************	Position 3	************************************************/
				double[] pos3id = new double[2];
				pos3id[0] = identity[2][0];
				pos3id[1] = identity[2][1];
				double[][] pos3 = williamsonEstimator(pos3id,site+2,low,mid,high);
				dn+=pos3[0][0];nn+=pos3[0][1];an+=pos3[0][2];fn+=pos3[0][3];
				ds+=pos3[1][0];ns+=pos3[1][1];as+=pos3[1][2];fs+=pos3[1][3];

			}
		}
		double[] finalmat = new double[4];

		finalmat[0]= fn - ((fs)*(neutralratio)*(1.0+(1.0/mu_ns))) ;  //adaptive fixations
		finalmat[1]= an - ((as)*(neutralratio)*(1.0+(1.0/mu_ns))) ;	//adaptive substitutions
		finalmat[2]= dn - ((ds)*(neutralratio)*(1.0+(1.0/mu_ns)));	//adaptive substitutions
		finalmat[3]= nn - ((ns)*(neutralratio)*(1.0+(1.0/mu_ns)));

		

		for(int i=0;i<finalmat.length;i++){
			if(finalmat[i]<0){
				finalmat[i]=0;
			}
		}
		double totalAdapt = finalmat[0]+finalmat[1];		
		double propAdapt = 1.0 - ((as+fs)/(an+fn))*neutralratio;
		this.low_A = finalmat[2];
		this.mid_A = finalmat[3];
		this.high_A = finalmat[1];
		this.fix_A = finalmat[0];
		this.low_R=dn;
		this.low_S=ds;
		this.mid_R =nn;
		this.mid_S =ns;
		this.high_R =an;
		this.high_S =as;
		this.fix_R =fn;
		this.fix_S =fs;
		this.Adapt=totalAdapt;
		this.Prop=propAdapt;
		this.L=integer_matrix.length;
	}

	
	// overloaded for no neutral ratio
	public double williamson3bin_method(double[] low, double[] mid, double[] high){
		double fs = 0;double fn = 0; double as = 0;double an = 0; double ds = 0; double dn = 0;double ns=0;double nn=0;
		for (int site = 0, codon = 0; site < integer_matrix[0].length - 2; site = site + 3, codon++) {
		//	double[][] identity = find_identityMK(site, codon);
			double[][] identity = NGmethodWilliamson(site, codon);
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { // check codon is not bad
				/************************	Position 1	************************************************/
				double[] pos1id = new double[2];
				pos1id[0] = identity[0][0];
				pos1id[1] = identity[0][1];
				double[][] pos1 = williamsonEstimator(pos1id,site,low,mid,high);
				dn+=pos1[0][0];nn+=pos1[0][1];an+=pos1[0][2];fn+=pos1[0][3];
				ds+=pos1[1][0];ns+=pos1[1][1];as+=pos1[1][2];fs+=pos1[1][3];

				/************************	Position 2	************************************************/
				double[] pos2id = new double[2];
				pos2id[0] = identity[1][0];
				pos2id[1] = identity[1][1];
				double[][] pos2 = williamsonEstimator(pos2id, site+1,low,mid,high);
				dn+=pos2[0][0];nn+=pos2[0][1];an+=pos2[0][2];fn+=pos2[0][3];
				ds+=pos2[1][0];ns+=pos2[1][1];as+=pos2[1][2];fs+=pos2[1][3];
				/************************	Position 3	************************************************/
				double[] pos3id = new double[2];
				pos3id[0] = identity[2][0];
				pos3id[1] = identity[2][1];
				double[][] pos3 = williamsonEstimator(pos3id,site+2,low,mid,high);
				dn+=pos3[0][0];nn+=pos3[0][1];an+=pos3[0][2];fn+=pos3[0][3];
				ds+=pos3[1][0];ns+=pos3[1][1];as+=pos3[1][2];fs+=pos3[1][3];

			}
		}
		double neutralratio  = nn/ns;
//		double[][] finalmat = new double[4][3];
//		finalmat[3][0]= fn - ((fs)*(neutralratio)) ;  //adaptive fixations
//		finalmat[3][1]= an - ((as)*(neutralratio)) ;	//adaptive substitutions
//		finalmat[3][2]= dn - ((ds)*(neutralratio));	//adaptive substitutions
//
//
//		for(int i=0;i<finalmat[0].length;i++){
//			if(finalmat[3][i]<0){
//				//	finalmat[3][i]=0;
//			}
//		}

        double[] finalmat = new double[4];
        finalmat[0]= fn - ((fs)*(neutralratio)) ;  //adaptive fixations
        finalmat[1]= an - ((as)*(neutralratio)) ;	//adaptive substitutions
        finalmat[2]= dn - ((ds)*(neutralratio));	//adaptive substitutions
        finalmat[3]= nn - ((ns)*(neutralratio));    // this will always be zero...



        for(int i=0;i<finalmat.length;i++){
            if(finalmat[i]<0){
                finalmat[i]=0;
            }
        }
        double totalAdapt = finalmat[0]+finalmat[1];
        double propAdapt = 1.0 - ((as+fs)/(an+fn))*neutralratio;
        this.low_A = finalmat[2];
        this.mid_A = finalmat[3];
        this.high_A = finalmat[1];
        this.fix_A = finalmat[0];
//        this.low_A = finalmat[3][2];
//        this.mid_A = finalmat[3];
//        this.high_A = finalmat[1];
//        this.fix_A = finalmat[0];
        this.low_R=dn;
        this.low_S=ds;
        this.mid_R =nn;
        this.mid_S =ns;
        this.high_R =an;
        this.high_S =as;
        this.fix_R =fn;
        this.fix_S =fs;
        this.Adapt=totalAdapt;
        this.Prop=propAdapt;
        this.L=integer_matrix.length;
        this.Nr = neutralratio;

		return totalAdapt;


	}


	//overloaded to include main gaps
	public void williamson3bin_method_IncludeMainGaps(double neutralratio,double[] low, double[] mid, double[] high){
		double fs = 0;double fn = 0; double as = 0;double an = 0; double ds = 0; double dn = 0;double ns=0;double nn=0;
		for (int site = 0, codon = 0; site < integer_matrix[0].length - 2; site = site + 3, codon++) {
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { // check codon is not bad
				double[][] identity = NGmethodWilliamson_IncludeMainGaps(site, codon);
				/************************	Position 1	************************************************/
				double[] pos1id = new double[2];
				pos1id[0] = identity[0][0];
				pos1id[1] = identity[0][1];
				double[][] pos1 = williamsonEstimator_IncludeMainGaps(pos1id,site,low,mid,high);
				dn+=pos1[0][0];nn+=pos1[0][1];an+=pos1[0][2];fn+=pos1[0][3];
				ds+=pos1[1][0];ns+=pos1[1][1];as+=pos1[1][2];fs+=pos1[1][3];

				/************************	Position 2	************************************************/
				double[] pos2id = new double[2];
				pos2id[0] = identity[1][0];
				pos2id[1] = identity[1][1];
				double[][] pos2 = williamsonEstimator_IncludeMainGaps(pos2id, site+1,low,mid,high);
				dn+=pos2[0][0];nn+=pos2[0][1];an+=pos2[0][2];fn+=pos2[0][3];
				ds+=pos2[1][0];ns+=pos2[1][1];as+=pos2[1][2];fs+=pos2[1][3];
				/************************	Position 3	************************************************/
				double[] pos3id = new double[2];
				pos3id[0] = identity[2][0];
				pos3id[1] = identity[2][1];
				double[][] pos3 = williamsonEstimator_IncludeMainGaps(pos3id,site+2,low,mid,high);
				dn+=pos3[0][0];nn+=pos3[0][1];an+=pos3[0][2];fn+=pos3[0][3];
				ds+=pos3[1][0];ns+=pos3[1][1];as+=pos3[1][2];fs+=pos3[1][3];

			}
		}
		double[][] finalmat = new double[4][3];
		finalmat[3][0]= fn - ((fs)*(neutralratio)) ;  //adaptive fixations
		finalmat[3][1]= an - ((as)*(neutralratio)) ;	//adaptive substitutions
		finalmat[3][2]= dn - ((ds)*(neutralratio));	//adaptive substitutions



		for(int i=0;i<finalmat[0].length;i++){
			if(finalmat[3][i]<0){
				finalmat[3][i]=0;
			}
		}
		double totalAdapt = finalmat[3][0]+finalmat[3][1];		
		double propAdapt = 1.0 - ((as+fs)/(an+fn))*neutralratio;
		this.low_R=dn;
		this.low_S=ds;
		this.mid_R =nn;
		this.mid_S =ns;
		this.high_R =an;
		this.high_S =as;
		this.fix_R =fn;
		this.fix_S =fs;
		this.Adapt=totalAdapt;
		this.Prop=propAdapt;
		this.L=integer_matrix.length;
	}
	// overloaded to exclude codons
	public double williamson3bin_method_ExcludeCodons(double neutralratio, int[] excludelist,double[] low, double[] mid, double[] high){
		double fs = 0;double fn = 0; double as = 0;double an = 0; double ds = 0; double dn = 0;double ns=0;double nn=0;
		for (int site = 0, codon = 0; site < integer_matrix[0].length - 2; site = site + 3, codon++) {
			double[][] identity = find_identityMK(site, codon);


			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false && excludelist[codon]>0) { // check codon is not bad
				/************************	Position 1	************************************************/
				double[] pos1id = new double[2];
				pos1id[0] = identity[0][0];
				pos1id[1] = identity[0][1];
				double[][] pos1 = williamsonEstimator(pos1id,site,low,mid,high);
				dn+=pos1[0][0];nn+=pos1[0][1];an+=pos1[0][2];fn+=pos1[0][3];
				ds+=pos1[1][0];ns+=pos1[1][1];as+=pos1[1][2];fs+=pos1[1][3];

				/************************	Position 2	************************************************/
				double[] pos2id = new double[2];
				pos2id[0] = identity[1][0];
				pos2id[1] = identity[1][1];
				double[][] pos2 = williamsonEstimator(pos2id, site+1,low,mid,high);
				dn+=pos2[0][0];nn+=pos2[0][1];an+=pos2[0][2];fn+=pos2[0][3];
				ds+=pos2[1][0];ns+=pos2[1][1];as+=pos2[1][2];fs+=pos2[1][3];
				/************************	Position 3	************************************************/
				double[] pos3id = new double[2];
				pos3id[0] = identity[2][0];
				pos3id[1] = identity[2][1];
				double[][] pos3 = williamsonEstimator(pos3id,site+2,low,mid,high);
				dn+=pos3[0][0];nn+=pos3[0][1];an+=pos3[0][2];fn+=pos3[0][3];
				ds+=pos3[1][0];ns+=pos3[1][1];as+=pos3[1][2];fs+=pos3[1][3];

			}
		}
		double[][] finalmat = new double[4][3];
		finalmat[3][0]= fn - ((fs)*(neutralratio)) ;  //adaptive fixations
		finalmat[3][1]= an - ((as)*(neutralratio)) ;	//adaptive substitutions
		finalmat[3][2]= dn - ((ds)*(neutralratio));	//adaptive substitutions

		for(int i=0;i<finalmat[0].length;i++){
			if(finalmat[3][i]<0){
				finalmat[3][i]=0;
			}
		}
		double totalAdapt = finalmat[3][0]+finalmat[3][1];
		//	return nn;
		return totalAdapt;

	}


	public int[] samplingArray(int blocksize, int length){
		//	double numblocks = length-blocksize+1; //overlapping
		double numblocks = length/blocksize; // non overlapping
		int[] sampler = new int[(int) numblocks];
		for(int x=0;x<sampler.length;x++){
			int randint = generator.nextInt((int)numblocks);//generator.nextInt((int) numblocks-1);
			sampler[x] = randint;
		}
		return sampler;
	}

	public BlockStruct[][] makeSeqBlocks(int blocksize,int length){
		double numblocks = length/blocksize;
		BlockStruct[][] blockmat = new BlockStruct[integer_matrix.length][(int) numblocks];
		int[] temp = new int[blocksize];
		for (int site = 0,x=0; site < length - (blocksize-1); site = site + blocksize,x++) {
			for(int i=0;i<integer_matrix.length;i++){
				int k=0;
				for(int j=site;j<site+blocksize;j++){
					temp[k] = integer_matrix[i][j];
					k++;
				}
				blockmat[i][x] = new BlockStruct(blocksize);
				for(int t=0;t<temp.length;t++){
					blockmat[i][x].sites[t] = temp[t];
				}
			}
		} 

		return blockmat;
	}

	public BlockStruct[] makeAnsBlocks(int blocksize,int length){
		double numblocks = length/blocksize;
		BlockStruct[] blockmat = new BlockStruct[(int) numblocks];

		int[] temp = new int[blocksize];
		for (int site = 0,x=0; site < length - (blocksize-1); site = site + blocksize,x++) {
			int k=0;
			for(int j=site;j<site+blocksize;j++){
				temp[k] = integer_ancestral[j];
				k++;
			}
			blockmat[x] = new BlockStruct(blocksize);
			for(int t=0;t<temp.length;t++){
				blockmat[x].sites[t] = temp[t];
			}
		}
		return blockmat;
	}

	public BlockStruct[][] makeSeqBlocksOverlapping(int blocksize,int length){
		int numblocks = length-blocksize+1;
		BlockStruct[][] blockmat = new BlockStruct[integer_matrix.length][(int) numblocks];
		int[] temp = new int[blocksize];
		for (int site = 0,x=0; site < length - (blocksize-1); site=site+3,x++) {
			for(int i=0;i<integer_matrix.length;i++){
				int k=0;
				for(int j=site;j<site+blocksize;j++){
					temp[k] = integer_matrix[i][j];
					k++;
				}
				blockmat[i][x] = new BlockStruct(blocksize);
				for(int t=0;t<temp.length;t++){
					blockmat[i][x].sites[t] = temp[t];
				}
			}
		}
		return blockmat;
	}

	public BlockStruct[] makeAnsBlocksOverlapping(int blocksize,int length){
		double numblocks = length-blocksize+1;
		BlockStruct[] blockmat = new BlockStruct[(int) numblocks];

		int[] temp = new int[blocksize];
		for (int site = 0,x=0; site < length - (blocksize-1); site = site+3,x++) {
			int k=0;
			for(int j=site;j<site+blocksize;j++){
				temp[k] = integer_ancestral[j];
				k++;
			}
			blockmat[x] = new BlockStruct(blocksize);
			for(int t=0;t<temp.length;t++){
				blockmat[x].sites[t] = temp[t];
			}
		}
		return blockmat;
	}

	public Store CreateBlocks(int blocksize,int length, int[] Sampling){
		Store MatrixStore = new Store();
		MatrixStore.BlockMatrix = makeSeqBlocks(blocksize,length);
		MatrixStore.AnsBlockMatrix = makeAnsBlocks(blocksize,length);		
		MatrixStore.RandomisedBlockMatrix = makeSeqBlocks(blocksize,length);
		MatrixStore.RandomisedAnsBlockMatrix = makeAnsBlocks(blocksize,length);	
		MatrixStore.SamplingArray =  Sampling;

		for(int j=0;j<MatrixStore.BlockMatrix[0].length;j++){
			MatrixStore.RandomisedAnsBlockMatrix[j] = MatrixStore.AnsBlockMatrix[MatrixStore.SamplingArray[j]];
			for(int i=0;i<MatrixStore.BlockMatrix.length;i++){
				MatrixStore.RandomisedBlockMatrix[i][j] = MatrixStore.BlockMatrix[i][MatrixStore.SamplingArray[j]];
			}
		}

		MatrixStore.RandomisedIntegerMatrix = new int[MatrixStore.RandomisedBlockMatrix.length][MatrixStore.RandomisedBlockMatrix[0].length*blocksize];
		MatrixStore.RandomisedIntegerAncestral = new int[MatrixStore.RandomisedAnsBlockMatrix.length*blocksize];
		int k=0;
		for(int j=0;j<MatrixStore.RandomisedAnsBlockMatrix.length;j++){
			for(int x=0;x<MatrixStore.RandomisedAnsBlockMatrix[j].sites.length;x++){
				MatrixStore.RandomisedIntegerAncestral[k] =  MatrixStore.RandomisedAnsBlockMatrix[j].sites[x];
				k++;
			}
		}


		for(int i=0;i<MatrixStore.RandomisedBlockMatrix.length;i++){
			k=0;
			for(int j=0;j<MatrixStore.RandomisedBlockMatrix[0].length;j++){
				for(int x=0;x<MatrixStore.RandomisedBlockMatrix[i][j].sites.length;x++){
					MatrixStore.RandomisedIntegerMatrix[i][k] =  MatrixStore.RandomisedBlockMatrix[i][j].sites[x];
					k++;
				}
			}
		}

		return MatrixStore;
	}

	public Store CreateBlocks(int blocksize,int length){
		Store MatrixStore = new Store();
		MatrixStore.BlockMatrix = makeSeqBlocks(blocksize,length);
		MatrixStore.AnsBlockMatrix = makeAnsBlocks(blocksize,length);		
		MatrixStore.RandomisedBlockMatrix = makeSeqBlocks(blocksize,length);
		MatrixStore.RandomisedAnsBlockMatrix = makeAnsBlocks(blocksize,length);	
		MatrixStore.SamplingArray =  samplingArray(blocksize,length);

		for(int j=0;j<MatrixStore.BlockMatrix[0].length;j++){
			MatrixStore.RandomisedAnsBlockMatrix[j] = MatrixStore.AnsBlockMatrix[MatrixStore.SamplingArray[j]];
			for(int i=0;i<MatrixStore.BlockMatrix.length;i++){
				MatrixStore.RandomisedBlockMatrix[i][j] = MatrixStore.BlockMatrix[i][MatrixStore.SamplingArray[j]];
			}
		}

		MatrixStore.RandomisedIntegerMatrix = new int[MatrixStore.RandomisedBlockMatrix.length][MatrixStore.RandomisedBlockMatrix[0].length*blocksize];
		MatrixStore.RandomisedIntegerAncestral = new int[MatrixStore.RandomisedAnsBlockMatrix.length*blocksize];
		int k=0;
		for(int j=0;j<MatrixStore.RandomisedAnsBlockMatrix.length;j++){
			for(int x=0;x<MatrixStore.RandomisedAnsBlockMatrix[j].sites.length;x++){
				MatrixStore.RandomisedIntegerAncestral[k] =  MatrixStore.RandomisedAnsBlockMatrix[j].sites[x];
				k++;
			}
		}


		for(int i=0;i<MatrixStore.RandomisedBlockMatrix.length;i++){
			k=0;
			for(int j=0;j<MatrixStore.RandomisedBlockMatrix[0].length;j++){
				for(int x=0;x<MatrixStore.RandomisedBlockMatrix[i][j].sites.length;x++){
					MatrixStore.RandomisedIntegerMatrix[i][k] =  MatrixStore.RandomisedBlockMatrix[i][j].sites[x];
					k++;
				}
			}
		}

		return MatrixStore;
	}

	public double[][] NGmethodWilliamson(int site, int codonsite){
		double[][] identity = new double[3][2];
		int[] ancestralbases = preprocess.codon_split(codon_ancestral[codonsite]);
		int[] mainbases = new int[3];
		int[] count = new int[3];
		for(int i=0;i<integer_matrix.length;i++){
			// adds codon bases
			mainbases[0]=integer_matrix[i][site];
			mainbases[1]=integer_matrix[i][site+1];
			mainbases[2]=integer_matrix[i][site+2];
			double[] tmp = NGpathway(ancestralbases,mainbases);
			if(tmp[0]!=2.0){
				identity[0][0] += tmp[0];
				count[0]++;
			}
			if(tmp[1]!=2.0){
				identity[1][0] += tmp[1];
				count[1]++;
			}
			if(tmp[2]!=2.0){
				identity[2][0] += tmp[2];
				count[2]++;
			}
		}
		identity[0][0]=identity[0][0]/count[0];
		identity[1][0]=identity[1][0]/count[1];
		identity[2][0]=identity[2][0]/count[2];
		identity[0][1] = 1.0-identity[0][0];
		identity[1][1] = 1.0-identity[1][0];
		identity[2][1] = 1.0-identity[2][0];
		return identity;
	}

	public double[][] NGmethodWilliamson_IncludeMainGaps(int site, int codonsite){
		double[][] identity = new double[3][2];
		int[] ancestralbases = preprocess.codon_split(codon_ancestral[codonsite]);
		int[] mainbases = new int[3];
		int[] count = new int[3];
		for(int i=0;i<integer_matrix.length;i++){
			// adds codon bases
			mainbases[0]=integer_matrix[i][site];
			mainbases[1]=integer_matrix[i][site+1];
			mainbases[2]=integer_matrix[i][site+2];
			double[] tmp = NGpathway_IncludeMainGaps(ancestralbases,mainbases);
			if(tmp[0]!=2.0 && tmp[0]!=3.0 ){
				identity[0][0] += tmp[0];
				count[0]++;
			}
			if(tmp[1]!=2.0 && tmp[1]!=3.0 ){
				identity[1][0] += tmp[1];
				count[1]++;
			}
			if(tmp[2]!=2.0 && tmp[2]!=3.0 ){
				identity[2][0] += tmp[2];
				count[2]++;
			}
		}
		
		SiteInfo ss = SiteInformation_IncludeMainGaps(site);
		// test routine
		if(ss.Case>1 && count[0]==0.0){
		System.out.println(site+"\t"+count[0]+"\t"+ss.Case+"\t"+ ss.totalNumBases+"\t" + ss.data[0].NObs+"\t" + ss.data[1].NObs+"\t" + ss.data[2].NObs+"\t" + ss.data[3].NObs);
		}
		
		identity[0][0]=identity[0][0]/count[0];
		identity[1][0]=identity[1][0]/count[1];
		identity[2][0]=identity[2][0]/count[2];
		identity[0][1] = 1.0-identity[0][0];
		identity[1][1] = 1.0-identity[1][0];
		identity[2][1] = 1.0-identity[2][0];
	//	if(count[0]==0){System.out.println(identity[0][0]+"\t"+site);}//identity[0][0]=0;identity[0][1]=0;};
	//	if(count[1]==0){identity[1][0]=0;identity[1][1]=0;};
	//	if(count[2]==0){identity[2][0]=0;identity[2][1]=0;};
		return identity;
	}

	public double[] NGpathway_IncludeMainGaps(int[] a, int[] b){ //a = ancestral b=main// 1 silent, 0 replacement, 2 invariant
		int degen = 0;

		int AnscodonNumber =  getcodonnumber(a[0],a[1],a[2]);
		int codonNumber =  getcodonnumber(b[0],b[1],b[2]);
		boolean[] diff = whichDiff(a,b);


		double[] identity = new double[3];
		identity[0]=3;identity[1]=3;identity[2]=3;
		if(b[0]<5 || b[1]<5 || b[2]<5){  // EXCESSIVELY IMPORTANT - if any gaps identity is 3
			identity[0]=2;identity[1]=2;identity[2]=2;  // if no gaps set identity to assume invariance i.e 2
			for(int i=0;i<3;i++){if(diff[i]){degen++;}} // find how many sites are degererate
			// 1 site different ****************************************************************************************
			if(degen==1){
				if(diff[0]){identity[0]=SilentOrReplacement(AnscodonNumber,codonNumber);}
				if(diff[1]){identity[1]=SilentOrReplacement(AnscodonNumber,codonNumber);}
				if(diff[2]){identity[2]=SilentOrReplacement(AnscodonNumber,codonNumber);}
			}
			// 2 sites different ****************************************************************************************
			if(degen==2){ 
				//			Pathway 1  (a) X11-X21-X22 and  (b) X11-X12-X22	  *Both pathways occur with equal probability*
				if(diff[0]==false){
					// (a) X11-X21-X22
					int a_codonNumber =  getcodonnumber(a[0],b[1],a[2]);
					identity[1]=SilentOrReplacement(AnscodonNumber,a_codonNumber);
					identity[2]=SilentOrReplacement(a_codonNumber,codonNumber);
					// (b) X11-X12-X22
					int b_codonNumber =  getcodonnumber(a[0],a[1],b[2]);
					identity[2]+=SilentOrReplacement(AnscodonNumber,b_codonNumber);
					identity[1]+=SilentOrReplacement(b_codonNumber,codonNumber);
					identity[1]=identity[1]/2;identity[2]=identity[2]/2; //average over all pathways
				}
				//			Pathway 2  (a) 1X1-2X1-2X2 and (b) 1X1-1X2-2X2		
				if(diff[1]==false){
					// (a) 1X1-2X1-2X2
					int a_codonNumber =  getcodonnumber(b[0],a[1],a[2]);
					identity[0]=SilentOrReplacement(AnscodonNumber,a_codonNumber);
					identity[2]=SilentOrReplacement(a_codonNumber,codonNumber);

					// (b) 1X1-1X2-2X2	
					int b_codonNumber =  getcodonnumber(a[0],a[1],b[2]);
					identity[2]+=SilentOrReplacement(AnscodonNumber,b_codonNumber);
					identity[0]+=SilentOrReplacement(b_codonNumber,codonNumber);
					identity[0]=identity[0]/2;identity[2]=identity[2]/2; //average over all pathways
				}
				//			Pathway 3   (a) 11X-21X-22X and  (b) 11X-12X-22X	
				if(diff[2]==false){
					// (a) 11X-21X-22X 
					int a_codonNumber =  getcodonnumber(b[0],a[1],a[2]);
					identity[0]=SilentOrReplacement(AnscodonNumber,a_codonNumber);
					identity[1]=SilentOrReplacement(a_codonNumber,codonNumber);
					// (b) 11X-12X-22X
					int b_codonNumber =  getcodonnumber(a[0],b[1],a[2]);
					identity[1]+=SilentOrReplacement(AnscodonNumber,b_codonNumber);
					identity[0]+=SilentOrReplacement(b_codonNumber,codonNumber);
					identity[0]=identity[0]/2;identity[1]=identity[1]/2; //average over all pathways
				}
			}

			// 3 sites different ****************************************************************************************
			if(degen==3){
				// Pathway 1 111-211-221-222	*x -(a)-(b)- x*        
				int a_codonNumber = getcodonnumber(b[0],a[1],a[2]);
				int b_codonNumber = getcodonnumber(b[0],b[1],a[2]);
				identity[0]=SilentOrReplacement(AnscodonNumber,a_codonNumber);
				identity[1]=SilentOrReplacement(a_codonNumber,b_codonNumber);
				identity[2]=SilentOrReplacement(b_codonNumber,codonNumber);
				// Pathway 2 111-211-212-222       
				a_codonNumber = getcodonnumber(b[0],a[1],a[2]);
				b_codonNumber = getcodonnumber(b[0],a[1],b[2]);
				identity[0]+=SilentOrReplacement(AnscodonNumber,a_codonNumber);
				identity[2]+=SilentOrReplacement(a_codonNumber,b_codonNumber);
				identity[1]+=SilentOrReplacement(b_codonNumber,codonNumber);				
				// Pathway 3 111-121-221-222    
				a_codonNumber = getcodonnumber(a[0],b[1],a[2]);
				b_codonNumber = getcodonnumber(b[0],b[1],a[2]);
				identity[1]+=SilentOrReplacement(AnscodonNumber,a_codonNumber);
				identity[0]+=SilentOrReplacement(a_codonNumber,b_codonNumber);
				identity[2]+=SilentOrReplacement(b_codonNumber,codonNumber);		
				// Pathway 4  111-121-122-222     
				a_codonNumber = getcodonnumber(a[0],b[1],a[2]);
				b_codonNumber = getcodonnumber(a[0],b[1],b[2]);
				identity[1]+=SilentOrReplacement(AnscodonNumber,a_codonNumber);
				identity[2]+=SilentOrReplacement(a_codonNumber,b_codonNumber);
				identity[0]+=SilentOrReplacement(b_codonNumber,codonNumber);		
				// Pathway 5 111-112-212-222      
				a_codonNumber = getcodonnumber(a[0],a[1],b[2]);
				b_codonNumber = getcodonnumber(b[0],a[1],b[2]);
				identity[2]+=SilentOrReplacement(AnscodonNumber,a_codonNumber);
				identity[0]+=SilentOrReplacement(a_codonNumber,b_codonNumber);
				identity[1]+=SilentOrReplacement(b_codonNumber,codonNumber);
				// Pathway 6 111-112-122-222       
				a_codonNumber = getcodonnumber(a[0],a[1],b[2]);
				b_codonNumber = getcodonnumber(a[0],b[1],b[2]);
				identity[2]+=SilentOrReplacement(AnscodonNumber,a_codonNumber);
				identity[1]+=SilentOrReplacement(a_codonNumber,b_codonNumber);
				identity[0]+=SilentOrReplacement(b_codonNumber,codonNumber);
				identity[0]=identity[0]/6.0;identity[1]=identity[1]/6.0;identity[2]=identity[2]/6.0;
			}
		}

		return identity;
	}

	public void gapInfo(){
		double[] invalidcount = new double[integer_matrix.length];
		double count=0;
		for(int i=0;i<integer_matrix.length;i++){
			count=0;
			for(int j=0;j<integer_matrix[0].length;j++){
				if(integer_matrix[i][j]==5){
					count++;
				}
			}
			invalidcount[i] = 100*(count/(double) integer_matrix[0].length);
		}
		this.gapcount=invalidcount;
	}

	public void gapSiteInfo(){
		double[] invalidcount = new double[integer_matrix[0].length];
		double count=0;
		for(int i=0;i<integer_matrix[0].length;i++){
			count=0;
			for(int j=0;j<integer_matrix.length;j++){
				if(integer_matrix[j][i]==5){
					count++;
				}
			}
			invalidcount[i] = (count/(double) integer_matrix.length);
		}
		this.gapSitecount=invalidcount;
	}
	public int[] BaseCount(){
		int count=0;
		int[] invalidcount = new int[integer_matrix[0].length];
		for(int i=0;i<integer_matrix[0].length;i++){
			count=0;
			for(int j=0;j<integer_matrix.length;j++){
				if(integer_matrix[j][i]<5){
					count++;
				}
			}
			invalidcount[i] = (count);
		}
		return invalidcount;

	}

	public void badsites_IncludeMainGaps(){
		boolean[] badsites = new boolean[integer_matrix[0].length];
		for (int j=0; j<integer_matrix[0].length; j++){
			int tot=0;
			if(integer_ancestral[j]==5){badsites[j]=true;}
			for(int i=0;i<integer_matrix.length;i++){
				if(integer_matrix[i][j]==5){tot++;}
			}
			if(tot==integer_matrix.length){badsites[j]=true;} // can change this to remove sites with insufficient data
		}
		this.bad_sites_list = badsites;
	}

	public void codonGapDegeneracy(){
		int one = 0;
		int two=0;
		int three=0;
		int zero=0;
		for(int i=0;i<integer_matrix.length;i++){
			for (int j=0; j<integer_matrix[0].length; j=j+3){
				int[] c = new int[3];
				if(integer_matrix[i][j]==5){c[0]=1;}
				if(integer_matrix[i][j+1]==5){c[1]=1;}
				if(integer_matrix[i][j+2]==5){c[2]=1;}
				int tot = c[0]+c[1]+c[2];
				if(tot==0){
					zero++;
				}
				if(tot==1){
					one++;
				}
				if(tot==2){
					two++;
				}
				if(tot==3){
					three++;
				}
			}
		}
		System.out.println("There are : "+zero+" Full codons");
		System.out.println("There are : "+one+" One gap codons");
		System.out.println("There are : "+two+" two gap codons");
		System.out.println("There are : "+three+" three gap codons");
	}



}
