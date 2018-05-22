package teaspoon.adaptation;

import java.util.Random;

/**
 * teaspoon.adaptation.Williamson class contains two main methods:
 *
 * 1) williamson_method
 * 2) eyre_walker_method
 *
 * Both of these methods return a finalmat matrix object, but NB of different sizes
 *
 * williamson_method finalmat double[4][3]
 *
 * [0][0] - ds = fixed synonymous mutations
 * [0][1] - rs = rare synonymous mutations
 * [0][2] - cs = common synonymous mutations
 * [1][0] - dn = fixed non-synonymous mutations
 * [1][1] - rn = rare synonymous mutations
 * [1][2] - cn = common non-synonymous mutations
 * [2][0] - replacement-silentProb ratio for fixed mutations
 * [2][1] - replacement-silentProb ratio for rare mutations
 * [2][2] - replacement-silentProb ration for common mutations
 * [3][0] - adaptive fixations
 * [3][1] - adaptive substitutions (from mutation)
 * [3][2] - total number of adaptations
 *
 *
 * eyre_walker_method finalmat double[6][2]
 *
 * [0][0] - fixed silentProb count
 * [0][1] - polymorphic silentProb count
 * [1][0] - fixed replacement count
 * [1][1] - polymorphic replacement count
 * [2][0] - silentProb-replacement ratio (neutral/polymorphic)
 * [2][1] - silentProb-replacement ration (fixed sites)
 * [3][0] -
 * [3][1] - proportion of sites that are adaptive?
 * [4][0] -
 * [4][1] -  r_1-(s_1/r_1)(r_<1/(s_<1 + 1))
 * [5][0] -
 * [5][1] - proportion of sites that are adaptive?
 *
 *
 */


public class Williamson extends SiteEstMulti{
	Random generator = new Random();
	public Williamson() {
		// default no-arg constructor
		throw new RuntimeException("ERROR: please input the intger matrix and the ancestral matrix");
	}

	public Williamson(int[][] m, int[] a){
		super(m,a);

	}


    // return matrix finalmat (double [][]) which contains 6 objects of length 3.
    //it calls upon MK test
    //

	public double[][] williamson_method(){
		double ds = 0;double dn = 0; double cs = 0;double cn = 0; double rs = 0; double rn = 0;
		McDonaldKreitman mk = new McDonaldKreitman(integer_matrix,integer_ancestral);
		for (int site = 0, codon = 0; site < integer_matrix[0].length - 2; site = site + 3, codon++) {
			double[][] identity = find_identityMK(site, codon);
			
			/********************************************************************************************************************************************************/
			// position 1
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) { // check codon is not bad
				SiteInformation Info = mk.SiteInformation(site);
				if(Info.polymorphismCase==1){ // invariant
					// do nothing as site has no information
				}
				if(Info.polymorphismCase==2){// fixed
					ds += 1.0*identity[0][0] ;
					dn += 1.0*identity[0][1] ;
				}
				if(Info.polymorphismCase==3){// 1 state derived and ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}
					rare=rare/y;
					rs += rare*identity[0][0];
					rn += rare*identity[0][1];							

					common=common/y;
					cs += common*identity[0][0];
					cn += common*identity[0][1];
				}
				if(Info.polymorphismCase==4){// 2 state derived no ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}

					ds += 0.5*identity[0][0] ;
					dn += 0.5*identity[0][1] ;

					rare=rare/y;
					common=common/y;
					if(common==0.0){
						rs += 0.5*rare*identity[0][0];
						rn += 0.5*rare*identity[0][1];		
					} else {
						rs += 0.25*rare*identity[0][0];
						rn += 0.25*rare*identity[0][1];								
					}

					if(rare==0.0){
						cs += 0.5*common*identity[0][0];
						cn += 0.5*common*identity[0][1];
					} else {
						cs += 0.25*rare*identity[0][0];
						cn += 0.25*rare*identity[0][1];	
					}


				}
				if(Info.polymorphismCase==5){// 2 state derived and ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}
					rare=rare/y;
					rs += rare*identity[0][0];
					rn += rare*identity[0][1];							

					common=common/y;
					cs += common*identity[0][0];
					cn += common*identity[0][1];		
				}
				if(Info.polymorphismCase==6){// 3 state derived no ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}

					ds += (1.0/3.0)*identity[0][0] ;
					dn += (1.0/3.0)*identity[0][1] ;
					rare=rare/y;
					common=common/y;
					if(common==0.0){
						rs += (2.0/3.0)*rare*identity[0][0];
						rn += (2.0/3.0)*rare*identity[0][1];		
					} else {
						rs += (1.0/3.0)*rare*identity[0][0];
						rn += (1.0/3.0)*rare*identity[0][1];								
					}

					if(rare==0.0){
						cs += (2.0/3.0)*common*identity[0][0];
						cn += (2.0/3.0)*common*identity[0][1];
					} else {
						cs += (1.0/3.0)*rare*identity[0][0];
						cn += (1.0/3.0)*rare*identity[0][1];	
					}
				}
				if(Info.polymorphismCase==7){// 3 state derived and ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}
					rare=rare/y;
					rs += rare*identity[0][0];
					rn += rare*identity[0][1];							

					common=common/y;
					cs += common*identity[0][0];
					cn += common*identity[0][1];

				}
				/********************************************************************************************************************************************************/
				// positions 2 *****************
				Info = mk.SiteInformation(site+1);
				if(Info.polymorphismCase==1){ // invariant
					// do nothing as site has no information
				}
				if(Info.polymorphismCase==2){// fixed
					ds += 1.0*identity[1][0] ;
					dn += 1.0*identity[1][1] ;
				}
				if(Info.polymorphismCase==3){// 1 state derived and ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}
					rare=rare/y;
					rs += rare*identity[1][0];
					rn += rare*identity[1][1];							

					common=common/y;
					cs += common*identity[1][0];
					cn += common*identity[1][1];
				}
				if(Info.polymorphismCase==4){// 2 state derived no ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}

					ds += 0.5*identity[1][0] ;
					dn += 0.5*identity[1][1] ;

					rare=rare/y;
					common=common/y;
					if(common==0.0){
						rs += 0.5*rare*identity[1][0];
						rn += 0.5*rare*identity[1][1];		
					} else {
						rs += 0.25*rare*identity[1][0];
						rn += 0.25*rare*identity[1][1];								
					}

					if(rare==0.0){
						cs += 0.5*common*identity[1][0];
						cn += 0.5*common*identity[1][1];
					} else {
						cs += 0.25*rare*identity[1][0];
						cn += 0.25*rare*identity[1][1];	
					}


				}
				if(Info.polymorphismCase==5){// 2 state derived and ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}
					rare=rare/y;
					rs += rare*identity[1][0];
					rn += rare*identity[1][1];							

					common=common/y;
					cs += common*identity[1][0];
					cn += common*identity[1][1];		
				}
				if(Info.polymorphismCase==6){// 3 state derived no ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}

					ds += (1.0/3.0)*identity[1][0] ;
					dn += (1.0/3.0)*identity[1][1] ;
					rare=rare/y;
					common=common/y;
					if(common==0.0){
						rs += (2.0/3.0)*rare*identity[1][0];
						rn += (2.0/3.0)*rare*identity[1][1];		
					} else {
						rs += (1.0/3.0)*rare*identity[1][0];
						rn += (1.0/3.0)*rare*identity[1][1];								
					}

					if(rare==0.0){
						cs += (2.0/3.0)*common*identity[1][0];
						cn += (2.0/3.0)*common*identity[1][1];
					} else {
						cs += (1.0/3.0)*rare*identity[1][0];
						cn += (1.0/3.0)*rare*identity[1][1];	
					}


				}
				if(Info.polymorphismCase==7){// 3 state derived and ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}
					rare=rare/y;
					rs += rare*identity[1][0];
					rn += rare*identity[1][1];							

					common=common/y;
					cs += common*identity[1][0];
					cn += common*identity[1][1];

				}


				/********************************************************************************************************************************************************/
				// positions 3 *****************
				Info = mk.SiteInformation(site+2);
				if(Info.polymorphismCase==1){ // invariant
					// do nothing as site has no information
				}
				if(Info.polymorphismCase==2){// fixed
					ds += 1.0*identity[2][0] ;
					dn += 1.0*identity[2][1] ;
				}
				if(Info.polymorphismCase==3){// 1 state derived and ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}
					rare=rare/y;
					rs += rare*identity[2][0];
					rn += rare*identity[2][1];							

					common=common/y;
					cs += common*identity[2][0];
					cn += common*identity[2][1];
				}
				if(Info.polymorphismCase==4){// 2 state derived no ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}

					ds += 0.5*identity[2][0] ;
					dn += 0.5*identity[2][1] ;

					rare=rare/y;
					common=common/y;
					if(common==0.0){
						rs += 0.5*rare*identity[2][0];
						rn += 0.5*rare*identity[2][1];		
					} else {
						rs += 0.25*rare*identity[2][0];
						rn += 0.25*rare*identity[2][1];								
					}

					if(rare==0.0){
						cs += 0.5*common*identity[2][0];
						cn += 0.5*common*identity[2][1];
					} else {
						cs += 0.25*rare*identity[2][0];
						cn += 0.25*rare*identity[2][1];	
					}

				}
				if(Info.polymorphismCase==5){// 2 state derived and ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}
					rare=rare/y;
					rs += rare*identity[2][0];
					rn += rare*identity[2][1];							

					common=common/y;
					cs += common*identity[2][0];
					cn += common*identity[2][1];		
				}
				if(Info.polymorphismCase==6){// 3 state derived no ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}

					ds += (1.0/3.0)*identity[2][0] ;
					dn += (1.0/3.0)*identity[2][1] ;
					rare=rare/y;
					common=common/y;
					if(common==0.0){
						rs += (2.0/3.0)*rare*identity[2][0];
						rn += (2.0/3.0)*rare*identity[2][1];		
					} else {
						rs += (1.0/3.0)*rare*identity[2][0];
						rn += (1.0/3.0)*rare*identity[2][1];								
					}

					if(rare==0.0){
						cs += (2.0/3.0)*common*identity[2][0];
						cn += (2.0/3.0)*common*identity[2][1];
					} else {
						cs += (1.0/3.0)*rare*identity[2][0];
						cn += (1.0/3.0)*rare*identity[2][1];	
					}


				}
				if(Info.polymorphismCase==7){// 3 state derived and ans
					double rare=0;
					double common=0;
					double y=0; // number derived nucleotides

					for(int i=0;i<4;i++){
						if(Info.siteData[i].numObs>0.0 && Info.siteData[i].numObs<0.5 && Info.siteData[i].inAncestral==false){
							rare++;
							y++;
						}
						if(Info.siteData[i].numObs>=0.5 && Info.siteData[i].numObs<1.0 && Info.siteData[i].inAncestral==false){
							common++;
							y++;
						}
					}
					rare=rare/y;
					rs += rare*identity[2][0];
					rn += rare*identity[2][1];							

					common=common/y;
					cs += common*identity[2][0];
					cn += common*identity[2][1];		
				}
			}
		}
		double[][] finalmat = new double[4][3];

		finalmat[0][0] = ds;
		finalmat[1][0] = dn;
		finalmat[0][1] = rs;
		finalmat[1][1] = rn;
		finalmat[0][2] = cs;
		finalmat[1][2] = cn;


		for(int i=0;i<finalmat[0].length;i++){
			if(finalmat[0][i] != 0){
				finalmat[2][i] = finalmat[1][i]/finalmat[0][i];  // replacement-silentProb ratio       
			} else{
				finalmat[2][i] = Double.NaN; 
			}		
		}
		if(rs==0 /*|| rn==0*/){
			finalmat[3][0] = Double.NaN; 
			finalmat[3][1] = Double.NaN; 
			finalmat[3][2] = Double.NaN; 
		} else {
			finalmat[3][0]= dn - ((ds)*(rn/rs)*(1.0+(1.0/rs)));   //adaptive fixations
			finalmat[3][1]= cn - ((cs)*(rn/rs)*(1.0+(1.0/rs))); 	//adaptive substitutions
		}
		for(int i=0;i<finalmat[0].length;i++){
			if(finalmat[3][i]<0){
				finalmat[3][i]=0;
			}
		}
		finalmat[3][2]= finalmat[3][0]+finalmat[3][1];// total adaptations


		return finalmat;
	}

	public double[][] eyrewalker_method(){
		double[] sigma = new double[2];
		double[] rho = new double[2];

		McDonaldKreitman mk = new McDonaldKreitman(integer_matrix,integer_ancestral);
		double[][] mat = mk.createContingencyNew();
		double[][] finalmat = new double[6][2];
		// less than 1
		finalmat[0][0] = mat[0][1];		// polymorphic silentProb count
		finalmat[1][0] = mat[1][1];		// polymorphic replacement count
		// equal to 1
		finalmat[0][1] = mat[0][0];		// fixed silentProb count
		finalmat[1][1] = mat[1][0];		// fixed replacement count

		sigma[0]=mat[0][1]; // polymorphic silentProb count
		sigma[1]=mat[0][0];	// fixed silentProb count
		rho[0]=mat[1][1]; //polymorphic replacement count
		rho[1]=mat[1][0]; //fixed replacement count

		finalmat[2][0] = rho[0]/(sigma[0]+1);  // silentProb-replacement ratio (neutral ratio)

		if(sigma[1]==0){
			finalmat[2][1] = Double.NaN;	
		}else {
			finalmat[2][1] = sigma[1]/(rho[1]); // silentProb-replacement ratio (fixed ratio)
		}

		finalmat[3][0] = Double.NaN;
		finalmat[3][1] = 1.0 - ((finalmat[2][1])*finalmat[2][0]);  // 1 - (s_1/s_<1)*(r_1/r_<1) or   s_1*r_1/(s_<1 * r_<1)
		if(finalmat[3][1]<0){
			finalmat[3][1]=Double.NaN;;
		}

		finalmat[4][1] = rho[1] - sigma[1]*(rho[0]/(sigma[0]+1));   //  r_1-(s_1/r_1)(r_<1/s_<1), where instead of r_<1/s_<1, r_<1/(s_<1 + 1) is used (to avoid when s_<1 is undefined
		if(finalmat[4][1]<0){
			finalmat[4][1]=0;
		}

		finalmat[5][1] = 1.0 - ((sigma[1]/rho[1])*(rho[0]/(sigma[0]+1)));   // proportion of adaptive mutations
		if(finalmat[5][1]<0){
			finalmat[5][1]=0;
		}

		return finalmat;	

	}




}

