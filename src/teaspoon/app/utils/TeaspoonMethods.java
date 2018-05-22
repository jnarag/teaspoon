package teaspoon.app.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import teaspoon.adaptation.Analysis;
import teaspoon.adaptation.BhattMethod;
import teaspoon.adaptation.DataSet;
import teaspoon.adaptation.SequenceInfo;
import teaspoon.adaptation.Williamson3bin;


public class TeaspoonMethods {
	private static double[] cof = {76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};



	//	**********************************************************************	
	//	**********************************************************************	
	//	**********************************************************************		
	//	Processing routines
	//	Count number of bases in either i or j of sequence integer matrix
	//	takes character values of base	
	public static double num_of_base(int[][] matrix, char base, int site){
		double count = 0.0;
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
				count++;						// counter
			}
		}
		return count;
	}

	//	overloaded to take integer values of base
	public static double num_of_base(int[][] matrix, int base, int site){
		double count = 0.0;
		for (int i=0; i< matrix.length; i++){
			if (matrix[i][site] == base){
				count++;						// counter
			}
		}
		return count;

	}
	//	Find which bases are present
	//	returns integer array 1=true 0=false
	public static int[] which_bases(int[][] matrix, int site){
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

	public static int[] consensusArray(int[][] integer_matrix){
		int[] consensus = new int[integer_matrix[0].length];
		double[] counter = new double[5];
		for(int site=0;site<integer_matrix[0].length;site++){
			// count numbers
			counter[0] = TeaspoonMethods.num_of_base(integer_matrix, 1, site);
			counter[1] = TeaspoonMethods.num_of_base(integer_matrix, 2, site);
			counter[2] = TeaspoonMethods.num_of_base(integer_matrix, 3, site);
			counter[3] = TeaspoonMethods.num_of_base(integer_matrix, 4, site);
			counter[4] = TeaspoonMethods.num_of_base(integer_matrix, 5, site);
			int length = counter.length;
			double max = counter[0];
			int position = 0;
			for (int i = 1; i < length; i++) {
				if (counter[i]>max) {
					max = counter[i];	// update max
					position = (i+1);
				}
			}
			//after the loop, min contains the minimum value,
			//position contains its position inside the array

			consensus[site] = position;
		}
		return consensus;
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
	public static double factorial(float n){
		if(n<=1.0){return 0.0;}
		double value = gammaln(n + 1.0f);
		return value;
	}



	//	**********************************************************************	
	//	ALL METHODS BELOW ARE OVERLOADED TO WORK FOR 2D AND 1D ARRAYS	
	//	**********************************************************************	
	//	**********************************************************************	
	//	Integer matrix	
	public static int[][] make_integer(int[][] codon_matrix){
		int[][] integer_matrix = new int[codon_matrix.length][codon_matrix[0].length*3];
		for(int i=0;i<integer_matrix.length;i++){
			for (int j=0 ,k=0; j< integer_matrix[0].length; j=j+3, k++){
				int[] split = codonSplitter(codon_matrix[i][k]);
				integer_matrix[i][j] = split[0];
				integer_matrix[i][j+1] = split[1];
				integer_matrix[i][j+2] = split[2];
			}
		}
		return integer_matrix;

	}

	public static int[] makeIntegerFromCodon(int[] codon_array){
		int[] integer_array = new int[codon_array.length*3];

		for (int j=0 ,k=0; j< integer_array.length; j=j+3, k++){
			//System.out.println(codon_array[k]);
			int[] split = codonSplitter(codon_array[k]);
			integer_array[j] = split[0];
			integer_array[j+1] = split[1];
			integer_array[j+2] = split[2];
		}

		return integer_array;

	}
	//	Codon Matrix 
	public static int[][] make_codon (int[][] integer_matrix){
		int[][] codon_matrix = new int[integer_matrix.length][integer_matrix[0].length/3];
		String codon;
		for (int i=0; i<integer_matrix.length; i++){
			for (int j=0 ,k=0; j< integer_matrix[0].length; j=j+3, k++){
				String pos1 = Integer.toString(integer_matrix[i][j]);	// position 1
				String pos2 = Integer.toString(integer_matrix[i][j+1]); // position 2
				String pos3 = Integer.toString(integer_matrix[i][j+2]); // position 3
				codon = new StringBuffer().append(pos1).append(pos2).append(pos3).toString();	// Join 3 bases to for a codon
				codon_matrix[i][k] = Integer.parseInt(codon); // make integer
			}
		}
		return codon_matrix;
	}
	public static int[] makeCodon(int[] integer_array){
		int[] codon_array = new int[integer_array.length/3];
		String codon;
		for (int j=0 ,k=0; j< integer_array.length; j=j+3, k++){
			String pos1 = Integer.toString(integer_array[j]);	// position 1
			String pos2 = Integer.toString(integer_array[j+1]);	// position 2
			String pos3 = Integer.toString(integer_array[j+2]);	// position 3
			codon = new StringBuffer().append(pos1).append(pos2).append(pos3).toString();// Join 3 bases to for a codon
			codon_array[k] = Integer.parseInt(codon);	// make integer
		}
		return codon_array;
	}

	//	**********************************************************************	
	//	**********************************************************************		
	//	removes bad ancestral numCodons (REMOVE INVALID CHARACTER)
	public static int[] remove_bad_ancestral_codons(int[] integer_array){
		ArrayList<Integer> badcodons = new ArrayList<Integer>();
		boolean flag=false;
		for(int i=0,k=0; i<integer_array.length-3;i=i+3,k++){
			// check if invalid character
			if(integer_array[i]==5 || integer_array[i+1]==5 || integer_array[i+2]==5){
				flag=true;
			}
			if(flag){
				Integer temp = new Integer(k);
				badcodons.add(temp);					// add this codon to a flagged list
			}
			flag = false;
		}
		int[] badlist = new int[badcodons.size()];
		for(int i = 0; i < badcodons.size(); i++) {
			Integer temp2 = badcodons.get(i);	// cast flagged sites to string
			badlist[i] = temp2.intValue();
		}
		return badlist;
	}

	public static int[][] good_codons(int[][] codon_matrix, int[] badcodons){
		if(badcodons.length==0){ return codon_matrix;}
		int[][] good_codon_matrix = new int[codon_matrix.length][(codon_matrix[0].length)-(badcodons.length)];
		boolean flag = false;
		int m = 0;
		for (int j = 0; j < codon_matrix[0].length; j++){
			flag = false;
			for(int k=0; k< badcodons.length;k++){
				if(j == badcodons[k]){
					flag = true;						// if codon is flagged as bad
					break;
				}
			}
			if(flag == false){
				for (int i = 0; i < codon_matrix.length; i++){
					good_codon_matrix[i][m] = codon_matrix[i][j];  // only add unflaged numCodons
				}
				m++;
			}
		}

		return good_codon_matrix;
	}
	public static int[] good_codons(int[] codon_array, int[] badcodons){
		if(badcodons.length==0){ return codon_array;}
		int[] good_codon_array = new int[codon_array.length-(badcodons.length)];
		boolean flag = false;
		int m = 0;
		for (int j = 0; j < codon_array.length; j++){
			flag = false;
			for(int k=0; k< badcodons.length;k++){
				if(j == badcodons[k]){
					flag = true;						// if codon is flagged as bad
					break;
				}
			}
			if(flag == false){
				good_codon_array[m] = codon_array[j];	// only add unflaged numCodons
				m = m + 1;
			}
		}
		return good_codon_array;
	}

	// Normal Integer matrix
	public static int[][] good_integer(int[][] codon_matrix, int[] badcodons){
		int[][] good_integer_matrix = new int[codon_matrix.length][((codon_matrix[0].length)-(badcodons.length))*3];
		int[][] good_codon_matrix = good_codons(codon_matrix,badcodons);
		int m = 0;
		for (int i = 0; i< good_codon_matrix.length; i++){
			m=0;
			for (int j = 0; j < good_codon_matrix[0].length; j++){
				int[] indv_codons = codonSplitter(good_codon_matrix[i][j]);
				for (int k=0; k< indv_codons.length; k++){
					good_integer_matrix[i][m] = indv_codons[k];
					m = m + 1;
				}
			}
		}
		return good_integer_matrix;
	}
	
	public static int[] good_integer(int[] codon_array, int[] badcodons){
		int[] good_integer_array = new int[((codon_array.length)-(badcodons.length))*3];
		int[] good_codon_array = good_codons(codon_array,badcodons);
		int m = 0;
		for (int j = 0; j < good_codon_array.length; j++){
			int[] indv_codons = codonSplitter(good_codon_array[j]);
			for (int k=0; k< indv_codons.length; k++){
				good_integer_array[m] = indv_codons[k];
				m = m + 1;
			}
		}
		return good_integer_array;
	}
	/*
	 *  Splits numCodons
	 */
	public static int[] codonSplitter(int codon){
		int[] indv_codons = new int[3];
		String temp = Integer.toString(codon);
		String[] temp_array = temp.split("");
		for (int i=0, k=0; i< temp_array.length; i++, k++){
			indv_codons[k] = Integer.parseInt(temp_array[i],10);
//			System.out.println(indv_codons[k]);
//			System.out.println(Integer.parseInt(temp_array[i]));

		}
		return indv_codons;
	}

	//	Find bad sites - Ones with gaps and numReplicates's - input integer objects
	public static int[] find_bad_sites(int[][] integer_matrix, int[] integer_array){
		boolean flag = false;
		ArrayList<Integer> badsites = new ArrayList<Integer>();
		for (int i = 0; i< integer_matrix[0].length; i++){	// for sites
			// remove any sites with gaps or invalid characters
			if(num_of_base(integer_matrix,5,i)>0){
				flag=true;
			}
			if(integer_array[i]==5){  // check anscetor does not have invalid sites
				flag =true;
			}
			//			if site flagged for any of the above reasons label as a bad site
			if(flag == true){
				Integer temp = new Integer(i);
				badsites.add(temp);					// add this site to a flagged list
			}
			flag = false;							// reset flag status
		}
		int[] badlist = new int[badsites.size()];
		for(int i = 0; i < badsites.size(); i++) {
			Integer temp2 = badsites.get(i);	// cast flagged sites to int
			badlist[i] = temp2.intValue();
		}
		return badlist;
	}

	public static boolean[] InvalidSites(int[][] integer_matrix, int[] integer_array){
		boolean flag = false;
		boolean[] badlist = new boolean[integer_matrix[0].length];
		for (int i = 0; i< integer_matrix[0].length; i++){	// for sites
			// flag any sites with gaps or invalid characters
			if(num_of_base(integer_matrix,5,i)>0){
				flag=true;  // check sequence alignment
			}
			if(integer_array[i]>4){  // check anscetor does not have invalid sites
				flag =true; // check ancestral sequence
			}
			//			if site flagged for any of the above reasons label as a bad site
			badlist[i]=flag;
			flag = false;							// reset flag status
		}
		return badlist;
	}


	//	boolean array showing if a site is bad or not	
	public static boolean[] bad_sites_list(int[][] good_integer, int[] good_ancestral){
		int[] badsites = find_bad_sites(good_integer,good_ancestral);
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
	//	number of bad sites
	public static int number_of_bad_sites(boolean[] bad_sites_list){
		int count=0;
		for(int i=0;i<bad_sites_list.length;i++){
			if(bad_sites_list[i]){
				count++;
			}
		}
		return count;
	}

	//	**********************************************************************		
	//	The accesor comands to access get all the matricies (Integer, Codon and Amino Acid)
	//	The choices for the user are: "base", "codon" and "AA". Only good numCodons are used	


	public static int[][] get_object(int[][] integer_matrix,int[] ancestral_array, String type){
		int[][] codonmat = make_codon(integer_matrix);
		int[] badlist = remove_bad_ancestral_codons(ancestral_array);
		if(type == "base"){					// returns base integer matrix
			if(badlist.length==0){
				return integer_matrix;
			} else{
				int[][] goodinteger = good_integer(codonmat,badlist);
				return goodinteger;
			}
		} else if(type == "codon"){			// returns codon matrix
			if(badlist.length ==0){
				return codonmat;
			} else{
				int[][] goodcodon = good_codons(codonmat,badlist);
				return goodcodon;
			}
		} else { throw new RuntimeException("Illegal choice, choose 'base','codon'");
		}
	}
	//	This overload takes the badlist and the int array and returns what is requited	
	public static int[] get_ancestral_object(int[][] integer_matrix,int[] ancestral_array, String type){
		int[] codonarray = makeCodon(ancestral_array);
		int[] badlist = remove_bad_ancestral_codons(ancestral_array);
		if(type == "base"){
			int[] goodinteger = good_integer(codonarray,badlist);
			return  goodinteger;
		} else if(type == "codon"){
			int[] goodcodon = good_codons(codonarray,badlist);
			return goodcodon;
		} else { throw new RuntimeException("Illegal choice, choose 'base','codon'");
		}
	}


	public static double[] EmpiricalFreqs(int[][] integer_matrix){
		double k=0;
		double[] val = new double[4];
		for(int i=0;i<integer_matrix.length;i++){
			for(int j=0;j<integer_matrix[0].length;j++){
				if(integer_matrix[i][j]==1){
					k++;val[0]++;
				}
				if(integer_matrix[i][j]==2){
					k++;val[1]++;
				}
				if(integer_matrix[i][j]==3){
					k++;val[2]++;
				}
				if(integer_matrix[i][j]==4){
					k++;val[3]++;
				}
			}
		}
		for(int i=0;i<val.length;i++){
			val[i] = val[i]/k;
		}
		return val;
	}

	// subsetter takes an array with numbers 1,2,3 which denote which group numCodons belong too. Then extracts any number
	// of groups
	public static int[][] subsetter(int[][] m, int[] list, int which){
		int[][] codon = make_codon(m);
		int length=0;

		for(int x=0;x<list.length;x++){
			if(list[x]==which){ //**********
				length++;
			}
		}
		int[][] subset = new int[codon.length][length];
		int k=0;
		//System.out.println(length);
		for(int x=0;x<list.length;x++){

			//System.out.println(x);
			if(list[x]==which){ //**********

				for(int i=0;i<codon.length;i++){
					subset[i][k] = codon[i][x];
					//System.out.println(i+":"+k+":"+x);
				}
				k++;
			}
		}

		int[][] subset_integer = make_integer(subset);
		return subset_integer;
	}

	public static int[] subsetter(int[] m, int[] list, int which){
		int[] codon = makeCodon(m);
		int length=0;
		for(int x=0;x<list.length;x++){
			if(list[x]==which){ //**********
				length++;
			}
		}
		int[] subset = new int[length];
		int k=0;
		for(int x=0;x<list.length;x++){
			//System.out.println(">"+x);

			if(list[x]==which){ //**********
				subset[k] = codon[x];
				//System.out.println(x+" : "+codon[x]);
				k++;
			}
		}

		int[] subset_integer = makeIntegerFromCodon(subset);
		return subset_integer;
	}

	public static void writeSubsetAlignment(String file,String output1,String output0, int[] Array){
		DataSet data = new DataSet(file);
		int count1=0;
		int count0=0;
		for(int i=0;i<Array.length;i++){
			if(Array[i]==1){
				count1++;
			}else if(Array[i]==0){
				count0++;
			}
		}
		if(count1+count0!=Array.length){
			throw new RuntimeException("Error! Array contains characters other than 1 and 0. please check.");
		}

		ArrayList<SequenceInfo> M0 = new ArrayList<SequenceInfo>();
		ArrayList<SequenceInfo> M1 = new ArrayList<SequenceInfo>();
		int[][] mat = make_codon(data.getIntegerMatrix());
		if(mat[0].length!=Array.length){
			System.out.println(mat[0].length+"\t"+Array.length);
			throw new RuntimeException("Error! Array length does not equal Matrix length");
		}
		int[][] mat0 = new int[mat.length][count0];
		int[][] mat1 = new int[mat.length][count1];


		for(int i=0;i<mat.length;i++){
			int k1=0,k0=0;
			for(int j=0;j<mat[0].length;j++){
				if(Array[j]==1){
					mat1[i][k1]=mat[i][j];
					k1++;
				}
				if(Array[j]==0) {
					mat0[i][k0]=mat[i][j];
					k0++;
				}
			}
		}
		int[][] matint1=make_integer(mat1);
		int[][] matint0=make_integer(mat0);


		for(int i=0;i<data.getIntegerMatrix().length;i++){
			SequenceInfo S0 = new SequenceInfo();
			SequenceInfo S1 = new SequenceInfo();
			S0.setTaxon(data.getTaxonMatrix()[i]);
			S0.setSequence(matint0[i]);
			S1.setTaxon(data.getTaxonMatrix()[i]);
			S1.setSequence(matint1[i]);
			M0.add(S0);
			M1.add(S1);
		}
		data.exportFASTA(output1, M1);
		data.exportFASTA(output0, M0);


	}


	public static int[][] subMatrix(int[][] sequenceMatrix, int start, int end, boolean gapLimit) {


		int window_length = end-start;
		List<int[]> sub_matrix = new ArrayList<int[]>();

		int limit = (int) Math.round(window_length*0.05);     // tolerates 5% of total sites with gaps

		int n = 0;
		List<String> cnames = new ArrayList<>();
		for (int[] aSequenceMatrix : sequenceMatrix) {



			int[] sub_sub_mat = Arrays.copyOfRange(aSequenceMatrix, start, end);

			if(gapLimit) {
				if (gapCheck(limit, sub_sub_mat)) {
					sub_matrix.add(sub_sub_mat);
				}
			}
			else{
				sub_matrix.add(sub_sub_mat);

			}

			n++;
		}

		int[][] sub_arrayOfarray = new int[sub_matrix.size()][window_length];

		if(sub_matrix.size()>30000) {

			List<Integer> indices = new ArrayList<Integer>();
			sub_arrayOfarray = new int[30000][window_length];
			List<Double> new_sub_sampleTimes = new ArrayList<Double>();
			for(int i = 0; i< sub_matrix.size(); i++) {

				indices.add(i);
			}

			Collections.shuffle(indices);

			for(int i = 0; i < 30000; i++) {

				sub_arrayOfarray[i] = sub_matrix.get(indices.get(i));
				//new_sub_sampleTimes.add(sub_sampleTimes.get(indices.get(i)));

			}
			//sub_sampleTimes.clear();
			//sub_sampleTimes.addAll(new_sub_sampleTimes);
		}
		else{

			sub_arrayOfarray = new int[sub_matrix.size()][window_length];


			for(int i=0; i < sub_matrix.size(); i++) {

				sub_arrayOfarray[i] = sub_matrix.get(i);
			}
		}


		return sub_arrayOfarray;
	}

	private static boolean gapCheck(int limit, int[] codonArray) {

		boolean gapLessThanLimit = true;
		int count = 0;
		for (int aCodonArray : codonArray) {

			if (aCodonArray == 5) { // 5 is 'gap'
				count++;
			}

			//only one gap allowed
			if (count > limit) {
				gapLessThanLimit = false;
				break;
			}

		}
		return gapLessThanLimit;
	}


	public static void record(StringBuffer sb, String dataset, double[] parameters, BhattMethod bm) {



		if(parameters.length==3) {
			double t = parameters[0];
			double d = parameters[1];
			double index = parameters[2];

			sb.append(dataset + "," + d + "," + (bm.getReplacementSubstitutionsCountArray()[(int) index] + bm.getSilentSubstitutionsCountArray()[(int) index]) + "," +
					bm.getSilentSubstitutionsCountArray()[(int) index] + "," + bm.getReplacementSubstitutionsCountArray()[(int) index] + "," + bm.getReplacementToSilentRatesRatio()[(int) index] + "," + bm.getNonNeutralSubstitutions()[(int) index] + "\n");
		}
		else{
			double t = parameters[0];
			double window = parameters[1];
			double d = parameters[2];
			double index = parameters[3];

			sb.append(dataset + "," + d + "," + window + "," + (bm.getReplacementSubstitutionsCountArray()[(int) index] + bm.getSilentSubstitutionsCountArray()[(int) index]) + "," +
					bm.getSilentSubstitutionsCountArray()[(int) index] + "," + bm.getReplacementSubstitutionsCountArray()[(int) index] + "," + bm.getReplacementToSilentRatesRatio()[(int) index] + "," + bm.getNonNeutralSubstitutions()[(int) index] + "\n");


		}

	}

	public static void record(StringBuffer sb, String dataset, double[] parameters, Williamson3bin w3b) {



		if(parameters.length==3) {
			double t = parameters[0];
			double d = parameters[1];
			double index = parameters[2];

			double r = 0; double s = 0; double nns = 0; //replacement, silentProb, non-neutral substitutions;
			if(index==0) {

				r = w3b.getLowR();
				s = w3b.getLowS();
				nns = w3b.getLowA();
			}
			else if(index==1) {
				r = w3b.getMidR();
				s = w3b.getMidS();
				nns = w3b.getMidA();
			}
			else if(index==2) {
				r = w3b.getHighR();
				s = w3b.getHighS();
				nns = w3b.getAdaptiveMutations(); //adaptive mutations in high-freq and fixn

			}

			double ratio = r/s;
			sb.append(dataset + "," + d + "," + (r + s) + "," + s + "," + r + "," + ratio + "," +nns + "\n");
		}
		else{
//			double t = parameters[0];
//			double window = parameters[1];
//			double d = parameters[2];
//			double index = parameters[3];
//
//			sb.append(dataset + "," + d + "," + window + "," + (bm.ReplacementCountArray[(int) index] + bm.SilentCountArray[(int) index]) + "," +
//					bm.SilentCountArray[(int) index] + "," + bm.ReplacementCountArray[(int) index] + "," + bm.ReplacementSilentRatio[(int) index] + "," + bm.NonNeutralSubstitutions[(int) index] + "\n");
//

		}

	}


	private static void correlation(String file) {
		String [][] rawinput = new String[25][];
		double [][] matrix = new double[101][24];
		try {
			BufferedReader reader1 = new BufferedReader(new FileReader(file));

			int i = 0;
			while(reader1.ready()) {

				String s = reader1.readLine();
				rawinput[i] = s.split(",");
				i++;

			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}


		for(int j = 0; j < 24; j++) {
			for(int i = 0; i<101; i++) {

				double value = Double.parseDouble(rawinput[j+1][i+1].trim());
				matrix[i][j] = value;
			}

		}


		double[] x = matrix[matrix.length-1];
		double[] bs_r = new double[100];
		for(int i=0; i<matrix.length-1; i++) {


			double r = new PearsonsCorrelation().correlation(x, matrix[i]);
			//System.out.println(r);
			bs_r[i] = r;

		}
		DescriptiveStatistics ds = new DescriptiveStatistics(bs_r);

		double lq = ds.getPercentile(25);
		double uq = ds.getPercentile(75);
		System.out.println(lq+","+uq);


	}

	private static void quartiles_per_group(String file, int groups) {
		String [][] rawinput = new String[25][];
		//double [][] matrix = new double[100][24];

		List<Map<Integer, double[]>> matrix = new ArrayList<Map<Integer, double[]>>();

		try {
			BufferedReader reader1 = new BufferedReader(new FileReader(file));

			int i = 0;
			while(reader1.ready()) {

				String s = reader1.readLine();
				rawinput[i] = s.split(",");
				//System.out.println(rawinput[i].length);
				i++;


			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		Map<Integer, String[][]> group_propAdapt = new HashMap<Integer, String[][]>();

		for(int i=0; i<groups; i++) {

			int count = 0;
			for(int j=0; j<24; j++) {

				if(Integer.parseInt(rawinput[j+1][101].trim()) == (i+1)) {

					count++;

				}
			}

			group_propAdapt.put(i+1, new String[count][]);

			int g = 0;
			for(int j=0; j<24; j++) {


				if(Integer.parseInt(rawinput[j+1][101].trim()) == (i+1)) {

					group_propAdapt.get(i + 1)[g] = rawinput[j + 1];
					g++;
				}


			}

			System.out.println("group "+(i+1));
			for(int k=0; k<group_propAdapt.get(i+1).length; k++){
				System.out.println("> "+(k+1));

				for(int j=0; j<100; j++){

//                    System.out.println(i+1);
//                    System.out.println(j+1);


					System.out.println((j+1)+":"+group_propAdapt.get(i+1)[k][j+1]);
				}
			}

		}

	}


	public static void readCodonMap(String csvfile, int n_codons, Analysis analysis) throws IOException {

		int [] map = new int[n_codons];
		try {
			BufferedReader reader = new BufferedReader(new FileReader(csvfile));

			int i = 0;
			while(reader.ready()) {

				String line = reader.readLine().trim();
				String [] parts = line.split(",");
				map[i] = Integer.parseInt(parts[1].trim());
				i++;
			}

			analysis.setCodonMap(map);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

	// codon_partition_name and the corresponding codon_group should be in the same order!

	public static Map<String, Integer> setWhichMap(String[] datasets_to_include, int[] which_key) {

		Map<String, Integer> which = new HashMap<>();

		for(int d=0; d < datasets_to_include.length; d++) {

			which.put(datasets_to_include[d], which_key[d]);
		}

		return which;
	}

//	public void fillCodonList(String[] datasets, int[] map, int[] keys, Analysis analysis) {
//
//		//if(datasets.length == maps.length) {
//
//		for (int i = 0; i < datasets.length; i++) {
//
//			//System.out.println(datasets[i]);
//			assert false;
//			analysis.codonList.put(datasets[i], map);
//            analysis.which.put(datasets[i], keys[i]);
//
//
//		}
//
//
//
//
//				//}
//	}



	private void writeoutBootstrapResults(String filename, String type, Analysis analysis) {

		assert false;

		String[] datasets = analysis.getDatasets();
		int no_datasets  = datasets.length;
		int no_timepoints;
		int bootstraps = analysis.getBootstraps();
		TeaspoonValues[][] value_matrix = analysis.getValue_matrix();

		String b_string = "";
		for(int i=0; i<bootstraps;i++){
			b_string += "B"+Integer.toString(i+1)+",";
		}
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(filename)));

			double lq = 0;
			double uq = 0;
			double med = 0;
			double std = 0;


			writer.write("Dataset,Time,No_Windows,Median,LQ,UQ,STD,,"+b_string+"\n");
			for (int d = 0; d < no_datasets; d++) {


				//System.out.println(bs.datasets[d]);
				no_timepoints = analysis.getTimepoints_per_dataset()[d];

				for (int t = 0; t < no_timepoints; t++) {


					DescriptiveStatistics window = new DescriptiveStatistics(value_matrix[t][d].getNumWindows());
					String s = String.valueOf(datasets[d])+","+value_matrix[t][d].getRow()+","+window.getMean();
					DescriptiveStatistics bstraps = null;

					switch (type.toLowerCase()) {

						case "s_high":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getSilentSubsHigh());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							//think about using arrayrealvectors


							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getSilentSubsHigh()[b];

							}

							writer.write(s+"\n");

							break;

						case "r_high":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getReplacementSubsHigh());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();
							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getReplacementSubsHigh()[b];

							}
							writer.write(s+"\n");

							break;

						case "a":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getAdaptations());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();
							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";
							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getAdaptations()[b];

							}

							writer.write(s+"\n");

							break;

						case "s_mid":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getSilentPolymorphsMidClass());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();
							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getSilentPolymorphsMidClass()[b];

							}

							writer.write(s+"\n");

							break;

						case "r_mid":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getReplacementPolymorphsMidClass());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getReplacementPolymorphsMidClass()[b];

							}

							writer.write(s+"\n");

							break;

						case "s_low":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getSilentPolymorphsLowClass());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();
							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getSilentPolymorphsLowClass()[b];

							}

							writer.write(s+"\n");

							break;

						case "r_low":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getReplacementPolymorphsLowClass());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getReplacementPolymorphsLowClass()[b];

							}

							writer.write(s+"\n");

							break;
						case "s_high/r_high":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getSilentToReplacementPolymorphsHighClassRatio());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getSilentToReplacementPolymorphsHighClassRatio()[b];

							}

							writer.write(s+"\n");

							break;

						case "r_high/s_high":

							bstraps = new DescriptiveStatistics(value_matrix[t][d].getReplacementToSilentPolymorphsHighClassRatio());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getReplacementToSilentPolymorphsHighClassRatio()[b];

							}

							writer.write(s+"\n");

							break;

						case "neutral_ratio":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getNeutralRatio());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getNeutralRatio()[b];

							}

							writer.write(s+"\n");

							break;

						case "r_low/s_low":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getReplacementToSilentPolymorphsLowClassRatio());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getReplacementToSilentPolymorphsLowClassRatio()[b];

							}

							writer.write(s+"\n");

							break;



						case "theta":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getTheta());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getTheta()[b];

							}

							writer.write(s+"\n");

							break;

						case "tajimasd":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getTajimasD());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getTajimasD()[b];

							}

							writer.write(s+"\n");

							break;

						case "s":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getWattersonS());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getWattersonS()[b];

							}

							writer.write(s+"\n");

							break;

						case "pi":
							bstraps = new DescriptiveStatistics(value_matrix[t][d].getWattersonPi());
							lq = bstraps.getPercentile(25);
							uq = bstraps.getPercentile(75);
							med = bstraps.getPercentile(50);
							std = bstraps.getStandardDeviation();

							s +=  "," + med + "," + lq + "," + uq + "," + std+ ",";

							for (int b = 0; b < bootstraps; b++) {

								s += "," + value_matrix[t][d].getWattersonPi()[b];

							}

							writer.write(s+"\n");

							break;
					}


				}

				System.out.println();

			}
			writer.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}



}
