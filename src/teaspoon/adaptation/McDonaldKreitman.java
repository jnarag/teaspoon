package teaspoon.adaptation;

import weka.core.ContingencyTables;

public class McDonaldKreitman extends SiteEstMulti {

	public McDonaldKreitman() {
		// default no-arg constructor
		throw new RuntimeException(
		"ERROR: please input the intger matrix and the ancestral matrix");
	}

	public McDonaldKreitman(int[][] m, int[] a) {
		super(m, a);
	}


	
	// original implementation method
	public double[][] createContingency() {
		double[][] contingencytable = new double[2][2];
		for (int site = 0, codon = 0; site < integer_matrix[0].length - 2; site = site + 3, codon++) {
			// counts number of bases that are the same as the ancestral base
			double[][] identity = find_identityMK(site, codon);
			// if site is not bad, and derived
			// position 1
			if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) {
				SiteInfo Info1 = SiteInformation(site);
				if(Info1.Case==1){ // invariant
					// do nothing as site has no information
				}
				else if(Info1.Case==2){// fixed // if nucleotide is fixed
					contingencytable[0][0] += identity[0][0] ;
					contingencytable[1][0] += identity[0][1] ;
				} else { // is polymorphic
					contingencytable[0][1] += identity[0][0] ;
					contingencytable[1][1] += identity[0][1] ;
					
				}
			// position 2
				SiteInfo Info2 = SiteInformation(site+1);
				if(Info2.Case==1){ // invariant
					// do nothing as site has no information
				}
				else if(Info2.Case==2){// fixed // if nucleotide is fixed
					contingencytable[0][0] += identity[1][0] ;
					contingencytable[1][0] += identity[1][1] ;
				} else { // is polymorphic
					contingencytable[0][1] += identity[1][0] ;
					contingencytable[1][1] += identity[1][1] ;					
				}
			// position 3
				SiteInfo Info3 = SiteInformation(site+2);
				if(Info3.Case==1){ // invariant
					// do nothing as site has no information
				}
				else if(Info3.Case==2){// fixed // if nucleotide is fixed
					contingencytable[0][0] += identity[2][0] ;
					contingencytable[1][0] += identity[2][1] ;
				} else { // is polymorphic
					contingencytable[0][1] += identity[2][0] ;
					contingencytable[1][1] += identity[2][1] ;
					
				}
			}
		}
		return contingencytable;
	}


	public double[][] createContingencyNew() {
		double[][] contingencytable = new double[2][2];

		for (int site = 0, codon = 0; site < integer_matrix[0].length - 2; site = site + 3, codon++) {
			double[][] identity = find_identityMK(site, codon);
			// position
				if (bad_sites_list[site] == false && bad_sites_list[site+1] == false && bad_sites_list[site+2] == false) {
					SiteInfo Info1 = SiteInformation(site);
					if(Info1.Case==1){ // invariant
						// do nothing as site has no information
					}
					if(Info1.Case==2){// fixed
						contingencytable[0][0] += 1.0*identity[0][0] ;
						contingencytable[1][0] += 1.0*identity[0][1] ;
					}
					if(Info1.Case==3){// 1 state derived and ans
						contingencytable[0][1] += 1.0*identity[0][0] ;
						contingencytable[1][1] += 1.0*identity[0][1] ;
					}
					if(Info1.Case==4){// 2 state derived no ans
						contingencytable[0][0] += 0.5*identity[0][0] ;
						contingencytable[1][0] += 0.5*identity[0][1] ;
						contingencytable[0][1] += 0.5*identity[0][0] ;
						contingencytable[1][1] += 0.5*identity[0][1] ;
					}
					if(Info1.Case==5){// 2 state derived and ans
						contingencytable[0][1] += 1.0*identity[0][0] ;
						contingencytable[1][1] += 1.0*identity[0][1] ;		
					}
					if(Info1.Case==6){// 3 state derived no ans
						contingencytable[0][0] += (1.0/3.0)*identity[0][0] ;
						contingencytable[1][0] += (1.0/3.0)*identity[0][1] ;
						contingencytable[0][1] += (2.0/3.0)*identity[0][0] ;
						contingencytable[1][1] += (2.0/3.0)*identity[0][1] ;

					}
					if(Info1.Case==7){// 3 state derived and ans
						contingencytable[0][1] += 1.0*identity[0][0] ;
						contingencytable[1][1] += 1.0*identity[0][1] ;	

					}
				
				
				// position 2
					SiteInfo Info2 = SiteInformation(site+1);
					if(Info2.Case==1){ // invariant
						// do nothing as site has no information
					}
					if(Info2.Case==2){// fixed
						contingencytable[0][0] += 1.0*identity[1][0] ;
						contingencytable[1][0] += 1.0*identity[1][1] ;
					}
					if(Info2.Case==3){// 1 state derived and ans
						contingencytable[0][1] += 1.0*identity[1][0] ;
						contingencytable[1][1] += 1.0*identity[1][1] ;
					}
					if(Info2.Case==4){// 2 state derived no ans
						contingencytable[0][0] += 0.5*identity[1][0] ;
						contingencytable[1][0] += 0.5*identity[1][1] ;
						contingencytable[0][1] += 0.5*identity[1][0] ;
						contingencytable[1][1] += 0.5*identity[1][1] ;
					}
					if(Info2.Case==5){// 2 state derived and ans
						contingencytable[0][1] += 1.0*identity[1][0] ;
						contingencytable[1][1] += 1.0*identity[1][1] ;		
					}
					if(Info2.Case==6){// 3 state derived no ans
						contingencytable[0][0] += (1.0/3.0)*identity[1][0] ;
						contingencytable[1][0] += (1.0/3.0)*identity[1][1] ;
						contingencytable[0][1] += (2.0/3.0)*identity[1][0] ;
						contingencytable[1][1] += (2.0/3.0)*identity[1][1] ;

					}
					if(Info2.Case==7){// 3 state derived and ans
						contingencytable[0][1] += 1.0*identity[1][0] ;
						contingencytable[1][1] += 1.0*identity[1][1] ;	
					}
				
				
				// position 3
					SiteInfo Info3 = SiteInformation(site+2);
					if(Info3.Case==1){ // invariant
						// do nothing as site has no information
					}
					if(Info3.Case==2){// fixed
						contingencytable[0][0] += 1.0*identity[2][0] ;
						contingencytable[1][0] += 1.0*identity[2][1] ;
					}
					if(Info3.Case==3){// 1 state derived and ans
						contingencytable[0][1] += 1.0*identity[2][0] ;
						contingencytable[1][1] += 1.0*identity[2][1] ;
					}
					if(Info3.Case==4){// 2 state derived no ans
						contingencytable[0][0] += 0.5*identity[2][0] ;
						contingencytable[1][0] += 0.5*identity[2][1] ;
						contingencytable[0][1] += 0.5*identity[2][0] ;
						contingencytable[1][1] += 0.5*identity[2][1] ;
					}
					if(Info3.Case==5){// 2 state derived and ans
						contingencytable[0][1] += 1.0*identity[2][0] ;
						contingencytable[1][1] += 1.0*identity[2][1] ;		
					}
					if(Info3.Case==6){// 3 state derived no ans
						contingencytable[0][0] += (1.0/3.0)*identity[2][0] ;
						contingencytable[1][0] += (1.0/3.0)*identity[2][1] ;
						contingencytable[0][1] += (2.0/3.0)*identity[2][0] ;
						contingencytable[1][1] += (2.0/3.0)*identity[2][1] ;
					}
					if(Info3.Case==7){// 3 state derived and ans
					contingencytable[0][1] += 1.0*identity[2][0] ;
					contingencytable[1][1] += 1.0*identity[2][1] ;	

					}
				}
			}
			return contingencytable;
	}

	public String[] labels() {
		String[] names = new String[12];
		names[0] = "Chi-squared probability: ";
		names[1] = "Chi-squared value: ";
		names[2] = "Cramer's V: ";
		names[3] = "Entropy conditioned on columns: ";
		names[4] = "Entropy conditioned on rows: ";
		names[5] = "Entropy conditioned on rows (with Laplace): ";
		names[6] = "Entropy of rows: ";
		names[7] = "Entropy of columns: ";
		names[8] = "Gain ratio: ";
		names[9] = "Negative log2 of multiple hypergeometric probability: ";
		names[10] = "Symmetrical uncertainty: ";
		names[11] = "Tau value: ";
		return names;

	}

	public double[] significanceTest(){
		double [][] matrix = createContingencyNew();
	//	double [][] matrix = createContingency();
		boolean yatescorrection = false;
		double n = matrix[0][0] + matrix[0][1] + matrix[1][0] + matrix[1][1];	// total
		// (a+b+c+d)
		if(n<200){yatescorrection = true;}	// for small sample sizes use yates
		// correction
		double[] values = new double[5];
		values[0]=matrix[0][0];
		values[1]=matrix[0][1];
		values[2]=matrix[1][0];
		values[3]=matrix[1][1];

		// using weka.core functions. Found at http://www.cs.waikhato.ac.nz/~ml/weka/
		values[4]= ContingencyTables.chiSquared(matrix, yatescorrection);
		
		/*	values[1]= ContingencyTables.chiVal(matrix, yatescorrection);
		values[2]= ContingencyTables.CramersV(matrix);
		values[3]= ContingencyTables.entropyConditionedOnColumns(matrix);
		values[4]= ContingencyTables.entropyConditionedOnRows(matrix);
		values[5]= ContingencyTables.entropyConditionedOnRows(matrix, matrix, 3);
		values[6]= ContingencyTables.entropyOverRows(matrix);
		values[7]= ContingencyTables.entropyOverColumns(matrix);
		values[8]= ContingencyTables.gainRatio(matrix);
		values[9]= ContingencyTables.log2MultipleHypergeometric(matrix);
		values[10]= ContingencyTables.symmetricalUncertainty(matrix);
		values[11]= ContingencyTables.tauVal(matrix); */

/*
		for (int site = 0, codon = 0; site < integer_matrix[0].length - 2; site = site + 3, codon++) {
			int[] whichbases = new int[4];
			for(int i=0; i< integer_matrix.length; i++){

				if(integer_matrix[i][site] == 1){
					whichbases[0] = 1;
				} else if(integer_matrix[i][site] == 2){
					whichbases[1] = 1;
				} else if(integer_matrix[i][site] == 3){
					whichbases[2] = 1;
				} else if(integer_matrix[i][site] == 4){
					whichbases[3] = 1;
				}
			}
			int total = whichbases[0] + whichbases[1] + whichbases[2] + whichbases[3]; 

			if(total==1){
				values[0]++;
			} if(total==2){
				values[1]++;
			} if(total==3){
				values[2]++;
			}if(total==4){
				values[3]++;
			}

		} */

		return values;
	}



}
