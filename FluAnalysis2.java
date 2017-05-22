package teaspoon;

import teaspoon.adaptation.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Date;
import java.util.Scanner;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
public class FluAnalysis2 {
	Random generator = new Random();
	public FluAnalysis2(){

	}




	public double[] GAPcounter(String Location, String Table,String Gene){

		DataSet d = new DataSet(Location,Table,Gene);
		double[] GapStore = new double[d.tablelength];
		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();
		int i=0;
		while(It.hasNext()){
			SequenceInfo element = It.next();
			GapStore[i] = element.Gap;
			i++;
		}
		return GapStore;
	}

	public String[] Acc(String Location, String Table,String Gene){
		DataSet d = new DataSet(Location,Table,Gene);
		String[] GapStore = new String[d.tablelength];


		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();
		int i=0;
		while(It.hasNext()){
			SequenceInfo element = It.next();
			GapStore[i] = element.Accession;
			i++;
		}
		return GapStore;
	}



	public int[] Datecounter(String Location, String Table, String Gene,String type){
		DataSet d = new DataSet(Location,Table,Gene);
		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();
		double x=0;
		int[] dates = new int[2011-1930];
		while(It.hasNext()){
			SequenceInfo element = It.next();
			x=element.Year-1931;

			if(element.Genotype.equals(type) && element.Gap<5.0){

				dates[(int)x]++;
			}
		}
		x=1931;
		return dates;

	}

	public int[] Datecounter2(String Location, String Table, String Gene,String type){
		DataSet d = new DataSet(Location,Table,Gene);
		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();
		double x=0;
		int[] dates = new int[2011-1930];
		while(It.hasNext()){
			SequenceInfo element = It.next();
			x=element.Year-1931;
			if(element.Genotype.equals(type) || element.Genotype.equals(type+"a") || element.Genotype.equals(type+"b") || element.Genotype.equals(type+"c") || element.Genotype.equals(type+"d") || element.Genotype.equals(type+"e") || element.Genotype.equals(type+"f") || element.Genotype.equals(type+"g")){
				if(element.Gap<5.0){
					dates[(int)x]++;
				}
			}
		}
		x=1931;
		return dates;

	}

	public int[] DatecounterContinentSplit(String Location, String Table, String Gene,String type){
		DataSet d = new DataSet(Location,Table,Gene);
		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();
		double x=0;
		int[] dates = new int[2011-1930];
		while(It.hasNext()){
			SequenceInfo element = It.next();
			x=element.Year-1931;
			if(element.Genotype.equals(type) || element.Genotype.equals(type+"a") || element.Genotype.equals(type+"b") || element.Genotype.equals(type+"c") || element.Genotype.equals(type+"d") || element.Genotype.equals(type+"e") || element.Genotype.equals(type+"f") || element.Genotype.equals(type+"g")){
				if(element.Gap<5.0){
					if(element.Continent.equals("EU")){
						dates[(int)x]++;
					}
				}
			}
		}
		x=1931;
		return dates;

	}

	public int[] Datecounter(String Location, String Table, String Gene,String[] type){
		DataSet d = new DataSet(Location,Table,Gene);
		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();
		double x=0;
		int[] dates = new int[2011-1930];
		while(It.hasNext()){
			SequenceInfo element = It.next();
			x=element.Year-1931;
			for(int i=0;i<type.length;i++){
				if(element.Genotype.equals(type[i]) && element.Gap<5.0){
					//			if(element.Genotype.equals(type[i])){	
					dates[(int)x]++;
				}
			}

		}
		x=1931;
		return dates;

	}


	//	runner which splits sequences 
	public void Runner(String Location, String Table, String Gene,String[] type,String[] Dates,int start,int stop,double nr,int which){
		// inputs if which is 0 this prints SR sounds if 1 prints adap
		// if stop is -1 use full length sequences, if -2 and -3 this does the N1 subsetting
		int isNA=stop;
		DataSet d = new DataSet(Location,Table,Gene);
		final double[] low = {0.0,0.15};
		final double[] mid = {0.15,0.75};
		final double[] high = {0.75,1.0};	

		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();

		/***************** Fix ancestral sequence *****************************/

		ArrayList<SequenceInfo> outans = new ArrayList<SequenceInfo>();
		while(It.hasNext()){		
			SequenceInfo element = It.next();
			boolean flag = false;
			for(int f=0;f<type.length;f++){
				if(element.Genotype.equals(type[f])){ // fiter genotypes
					flag=true;
				}
			}
			String delimiter = "-";
			String[] temp = Dates[0].split(delimiter);
			boolean flag2 = false; // second flag checks through dates
			for(int n=0;n<temp.length;n++){
				if(element.Year==Double.valueOf(temp[n])){
					flag2=true;
				}
			}

			if(flag2==true && element.Gap<5.0 && flag==true){  // filter dates and gap counts and types
				outans.add(element);
			}
		}

		int[][] ans1 = new int[outans.size()][d.DataMainFrame.get(0).Sequence.length];
		for(int i=0;i<outans.size();i++){
			ans1[i] = outans.get(i).Sequence;
		}
		int[] ans = consensusArray(ans1);


		/************************** Specify MAIN SEQUENCE************************************************************ */
		double dateAvg = 0; // calculate date average
		for(int Z=1;Z<Dates.length;Z++){  // BIG loop through dates ************************ change date in loops
			It =  d.DataMainFrame.iterator();
			double count=0;
			dateAvg = 0;
			ArrayList<int[]> Sequences = new ArrayList<int[]>(); 
			ArrayList<String> Names = new ArrayList<String>(); 
			ArrayList<SequenceInfo> out = new ArrayList<SequenceInfo>();
			while(It.hasNext()){		
				SequenceInfo element = It.next();
				boolean flag = false;  // first flag checks through genotypes
				for(int f=0;f<type.length;f++){
					if(element.Genotype.equals(type[f])){ // fiter genotypes
						flag=true;
					}
				}
				String delimiter = "-";
				String[] temp = Dates[Z].split(delimiter);
				boolean flag2 = false; // second flag checks through dates
				for(int n=0;n<temp.length;n++){
					if(element.Year==Double.valueOf(temp[n])){
						flag2=true;
					}
				}
				// adds the sequences we want according to date,genotype and gap uncertainty
				if(element.Gap<5.0 && flag==true && flag2==true){	// filter sequences with gap count more than five
					Sequences.add(element.Sequence); //saves sequence
					Names.add(element.Accession); //saves name
					out.add(element); //for output if needed
					dateAvg+=element.DecimalDate; //calculates average date
					count++;
				}



			}
			// make integer matricies dates ********************************************************************************			
			dateAvg = dateAvg/count;
			int[][] integer_matrix = new int[Sequences.size()][Sequences.get(0).length];
			String[] accessions = new String[Names.size()];
			// add sequences
			Iterator<int[]> It2 =  Sequences.iterator(); // converts object into a sequence matrix
			int l=0;
			while(It2.hasNext()){
				int[] tmp = It2.next();
				integer_matrix[l] = tmp ;
				l++;
			}
			Iterator<String> It3 =  Names.iterator(); //converts object into the accessions matrix
			l=0;
			while(It3.hasNext()){
				String tmp = It3.next();
				accessions[l] = tmp ;
				l++;
			}

			if(stop<=-1){			
				stop = integer_matrix[0].length;			
			}
			// adjust matrix size
			int[][] mat = new int[integer_matrix.length][stop-start];  //start=51, stop=1035 for ha1, ha2 is start 1036 - end
			int[] matans = new int[stop-start];

			for(int i=0;i<integer_matrix.length;i++){
				int p=0;
				for(int j=start;j<stop;j++){		
					mat[i][p] = integer_matrix[i][j];
					matans[p] = ans[j];
					p++;
				}
			}

			Neuraminidase Nd = new Neuraminidase();

			if(isNA==-2){ //this is N1
				Methods m = new Methods();
				int[][] mat2 = m.Subsetter(mat, Nd.badlisth1n1_NoStopCodon, 1);
				int[] matans2 = m.Subsetter(matans, Nd.badlisth1n1_NoStopCodon, 1);
				Williamson3bin ws = new Williamson3bin(mat2,matans2);

				ws.williamson3bin_method(nr,low, mid, high);

				if(which==1){
					System.out.print(dateAvg  + "\t" + ws.integer_matrix.length + "\t" + ws.low_R+"\t"+ws.low_S + "\t" + ws.Neut_R+"\t"+ws.Neut_S + "\t" + ws.High_R+"\t"+ws.High_S+"\n");
				} else {
					System.out.print(dateAvg + "\t" +  ws.LowA+ "\t" +  ws.MidA+ "\t" +  ws.Adapt + "\n" );

				}

			}
			else if(isNA==-3){ //this is N1
				Methods m = new Methods();
				int[][] mat3 = m.Subsetter(mat, Nd.badlisth1n1_NoStopCodon, 0);
				int[] matans3 = m.Subsetter(matans, Nd.badlisth1n1_NoStopCodon, 0);
				Williamson3bin ws = new Williamson3bin(mat3,matans3);

				ws.williamson3bin_method(nr,low, mid, high);

				if(which==1){
					System.out.print(dateAvg  + "\t" + ws.integer_matrix.length + "\t" + ws.low_R+"\t"+ws.low_S + "\t" + ws.Neut_R+"\t"+ws.Neut_S + "\t" + ws.High_R+"\t"+ws.High_S+"\n");
				} else {
					System.out.print(dateAvg + "\t" +  ws.LowA+ "\t" +  ws.MidA+ "\t" +  ws.Adapt + "\n" );

				}
			} else { // this is for everything else
				Williamson3bin ws = new Williamson3bin(mat,matans);

				ws.williamson3bin_method(nr,low, mid, high);

				if(which==1){
					System.out.print(dateAvg  + "\t" + ws.integer_matrix.length + "\t" + ws.low_R+"\t"+ws.low_S + "\t" + ws.Neut_R+"\t"+ws.Neut_S + "\t" + ws.High_R+"\t"+ws.High_S+"\n");
				} else {
					System.out.print(dateAvg + "\t" +  ws.LowA+ "\t" +  ws.MidA+ "\t" +  ws.Adapt + "\n" );

				}

			}

		}

	}
	
	
	public void WriteSequencesGeno(String outloc,String Location, String Table, String Gene,String[] type,String[] Dates,int start,int stop,double nr,int Num,int[] list,int oneorzero,String Geno){
		int WhichSub=stop;
		DataSet d = new DataSet(Location,Table,Gene);


		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();
		
		ArrayList<SequenceInfo> out = new ArrayList<SequenceInfo>();

		for(int Z=0;Z<Dates.length;Z++){  // BIG loop through dates ************************ change date in loops
			It =  d.DataMainFrame.iterator();
			double count=0;
			ArrayList<String> Names = new ArrayList<String>(); 
			while(It.hasNext()){		
				SequenceInfo element = It.next();
				boolean flag = false;  // first flag checks through genotypes
				for(int f=0;f<type.length;f++){
					if(element.Genotype.equals(type[f])){ // fiter genotypes
						flag=true;
					}
				}
				String delimiter = "-";
				String[] temp = Dates[Z].split(delimiter);
				boolean flag2 = false; // second flag checks through dates
				for(int n=0;n<temp.length;n++){
					if(element.Year==Double.valueOf(temp[n])){
						flag2=true;
					}
				}
				// adds the sequences we want according to date,genotype and gap uncertainty
				if(element.Gap<5.0 && flag==true && flag2==true){	// filter sequences with gap count more than five
					if(WhichSub==-2){ //this is N1
						Methods m = new Methods();
						element.Sequence=m.Subsetter(element.Sequence, list, oneorzero);
					} 					
					Names.add(element.Accession); //saves name
					System.out.println(element.Accession);
					out.add(element); //for output if needed
					count++;
				}
			}

		}
	//	d.exportFASTA(outloc+"_"+Geno+"_"+Gene+"_"+Dates[Z]+".fa", out);

	}
	
	
	public void WriteSequences(String outloc,String Location, String Table, String Gene,String[] type,String[] Dates,int start,int stop,double nr,int Num,int[] list,int oneorzero,String Geno){
		int WhichSub=stop;
		WilliamsonOutputClass woc = new WilliamsonOutputClass(Dates.length,Num);
		DataSet d = new DataSet(Location,Table,Gene);
		final double[] low = {0.0,0.15};
		final double[] mid = {0.15,0.75};
		final double[] high = {0.75,1.0};	

		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();
		
		ArrayList<SequenceInfo> outans = new ArrayList<SequenceInfo>();
		while(It.hasNext()){		
			SequenceInfo element = It.next();
			boolean flag = false;
			for(int f=0;f<type.length;f++){
				if(element.Genotype.equals(type[f])){ // fiter genotypes
					flag=true;
				}
			}
			String delimiter = "-";
			String[] temp = Dates[0].split(delimiter);
			boolean flag2 = false; // second flag checks through dates
			for(int n=0;n<temp.length;n++){
				if(element.Year==Double.valueOf(temp[n])){
					flag2=true;
				}
			}

			if(flag2==true && element.Gap<5.0 && flag==true){  // filter dates and gap counts and types
				
				if(WhichSub==-2){ //this is N1
					Methods m = new Methods();
					element.Sequence=m.Subsetter(element.Sequence, list, oneorzero);
				} 
				outans.add(element);

			}
		}
		d.exportFASTA(outloc+"_"+Geno+"_"+Gene+"_"+Dates[0]+".fa", outans);
		
		for(int Z=1;Z<Dates.length;Z++){  // BIG loop through dates ************************ change date in loops
			It =  d.DataMainFrame.iterator();
			double count=0;
			ArrayList<int[]> Sequences = new ArrayList<int[]>(); 
			ArrayList<String> Names = new ArrayList<String>(); 
			ArrayList<SequenceInfo> out = new ArrayList<SequenceInfo>();
			while(It.hasNext()){		
				SequenceInfo element = It.next();
				boolean flag = false;  // first flag checks through genotypes
				for(int f=0;f<type.length;f++){
					if(element.Genotype.equals(type[f])){ // fiter genotypes
						flag=true;
					}
				}
				String delimiter = "-";
				String[] temp = Dates[Z].split(delimiter);
				boolean flag2 = false; // second flag checks through dates
				for(int n=0;n<temp.length;n++){
					if(element.Year==Double.valueOf(temp[n])){
						flag2=true;
					}
				}
				// adds the sequences we want according to date,genotype and gap uncertainty
				if(element.Gap<5.0 && flag==true && flag2==true){	// filter sequences with gap count more than five
					if(WhichSub==-2){ //this is N1
						Methods m = new Methods();
						element.Sequence=m.Subsetter(element.Sequence, list, oneorzero);
					} 					
					Names.add(element.Accession); //saves name
					out.add(element); //for output if needed
					count++;
				}
			}
			d.exportFASTA(outloc+"_"+Geno+"_"+Gene+"_"+Dates[Z]+".fa", out);

		}
		
	}


	// this is the main output file // needs to be double checked
	public WilliamsonOutputClass RunnerBS(String Location, String Table, String Gene,String[] type,String[] Dates,int start,int stop,double nr,int Num,int[] list,int oneorzero){
		int WhichSub=stop;
		WilliamsonOutputClass woc = new WilliamsonOutputClass(Dates.length,Num);
		DataSet d = new DataSet(Location,Table,Gene);
		final double[] low = {0.0,0.15};
		final double[] mid = {0.15,0.75};
		final double[] high = {0.75,1.0};	

		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();

		/***************** Fix ancestral sequence *****************************/

		ArrayList<SequenceInfo> outans = new ArrayList<SequenceInfo>();
		while(It.hasNext()){		
			SequenceInfo element = It.next();
			boolean flag = false;
			for(int f=0;f<type.length;f++){
				if(element.Genotype.equals(type[f])){ // fiter genotypes
					flag=true;
				}
			}
			String delimiter = "-";
			String[] temp = Dates[0].split(delimiter);
			boolean flag2 = false; // second flag checks through dates
			for(int n=0;n<temp.length;n++){
				if(element.Year==Double.valueOf(temp[n])){
					flag2=true;
				}
			}

			if(flag2==true && element.Gap<5.0 && flag==true){  // filter dates and gap counts and types
				outans.add(element);
			}
		}


		int[][] ans1 = new int[outans.size()][d.DataMainFrame.get(0).Sequence.length];
		for(int i=0;i<outans.size();i++){
			ans1[i] = outans.get(i).Sequence;
		}
		int[] ans = consensusArray(ans1);

		for(int run=0;run<Num;run++){
			// boot strap stuff
			int[] sampler = new int[ans.length/3];
			int numsites=0;

			if(WhichSub==-1){
				sampler = new int[ans.length/3];			
				for(int x=0;x<sampler.length;x++){
					int randint = generator.nextInt(sampler.length-1);			
					sampler[x] = randint;
				}
			} else if(WhichSub==-2){
				for(int j=0;j<list.length;j++){if(list[j]==oneorzero){numsites++;}}
				sampler = new int[(numsites)];
				for(int x=0;x<sampler.length;x++){
					int randint = generator.nextInt(sampler.length-1);			
					sampler[x] = randint;
				}
			}

			/************************** Specify MAIN SEQUENCE************************************************************ */
			double dateAvg = 0; // calculate date average
			for(int Z=0;Z<Dates.length;Z++){  // BIG loop through dates ************************ change date in loops
				It =  d.DataMainFrame.iterator();
				double count=0;
				dateAvg = 0;
				ArrayList<int[]> Sequences = new ArrayList<int[]>(); 
				ArrayList<String> Names = new ArrayList<String>(); 
				ArrayList<SequenceInfo> out = new ArrayList<SequenceInfo>();
				while(It.hasNext()){		
					SequenceInfo element = It.next();
					boolean flag = false;  // first flag checks through genotypes
					for(int f=0;f<type.length;f++){
						if(element.Genotype.equals(type[f])){ // fiter genotypes
							flag=true;
						}
					}
					String delimiter = "-";
					String[] temp = Dates[Z].split(delimiter);
					boolean flag2 = false; // second flag checks through dates
					for(int n=0;n<temp.length;n++){
						if(element.Year==Double.valueOf(temp[n])){
							flag2=true;
						}
					}
					// adds the sequences we want according to date,genotype and gap uncertainty
					if(element.Gap<5.0 && flag==true && flag2==true){	// filter sequences with gap count more than five
						Sequences.add(element.Sequence); //saves sequence
						Names.add(element.Accession); //saves name
						out.add(element); //for output if needed
						dateAvg+=element.DecimalDate; //calculates average date
						count++;
					}



				}
				// make integer matricies dates ********************************************************************************			
				dateAvg = dateAvg/count;
				int[][] integer_matrix = new int[Sequences.size()][Sequences.get(0).length];
				String[] accessions = new String[Names.size()];
				// add sequences
				Iterator<int[]> It2 =  Sequences.iterator(); // converts object into a sequence matrix
				int l=0;
				while(It2.hasNext()){
					int[] tmp = It2.next();
					integer_matrix[l] = tmp ;
					l++;
				}
				Iterator<String> It3 =  Names.iterator(); //converts object into the accessions matrix
				l=0;
				while(It3.hasNext()){
					String tmp = It3.next();
					accessions[l] = tmp ;
					l++;
				}

				if(stop<=-1){
					stop = integer_matrix[0].length;
				}
				// adjust matrix size
				int[][] mat = new int[integer_matrix.length][stop-start];  //start=51, stop=1035 for ha1, ha2 is start 1036 - end
				int[] matans = new int[stop-start];

				for(int i=0;i<integer_matrix.length;i++){
					int p=0;
					for(int j=start;j<stop;j++){		
						mat[i][p] = integer_matrix[i][j];
						matans[p] = ans[j];
						p++;
					}
				}
				if(WhichSub==-2){ //this is N1
					Methods m = new Methods();
					int[][] mat2 = m.Subsetter(mat, list, oneorzero);
					int[] matans2 = m.Subsetter(matans, list, oneorzero);
					Williamson3bin ww = new Williamson3bin(mat2,matans2); //******
					Store s = ww.CreateBlocks(3,mat2[0].length,sampler); //******
					Williamson3bin w = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);

					w.williamson3bin_method(nr,low, mid, high);


					if(run==0 && Z==0){
						woc.Low_S.add(Double.valueOf(0.0));
						woc.Low_R.add(Double.valueOf(0.0));
						woc.Mid_S.add(Double.valueOf(0.0));
						woc.Mid_R.add(Double.valueOf(0.0));
						woc.High_S.add(Double.valueOf(0.0));
						woc.High_R.add(Double.valueOf(0.0));
						woc.Fix_S.add(Double.valueOf(0.0));
						woc.Fix_R.add(Double.valueOf(0.0));
						woc.dateavg.add(dateAvg);
						woc.LowA.add(Double.valueOf(0.0));
						woc.MidA.add(Double.valueOf(0.0));
						woc.HighA.add(Double.valueOf(0.0));
						woc.FixA.add(Double.valueOf(0.0));
						woc.TotalA.add(Double.valueOf(0.0));
						woc.numSamples.add(mat.length);
						woc.neut=nr;
						woc.L=matans2.length;
					}
					else if(run==0){
						Williamson3bin ws = new Williamson3bin(mat2,matans2);
						ws.williamson3bin_method(nr,low, mid, high);
						woc.Low_S.add(ws.low_S);
						woc.Low_R.add(ws.low_R);
						woc.Mid_S.add(ws.Neut_S);
						woc.Mid_R.add(ws.Neut_R);
						woc.High_S.add(ws.High_S);
						woc.High_R.add(ws.High_R);
						woc.Fix_S.add(ws.Fix_S);
						woc.Fix_R.add(ws.Fix_R);
						woc.dateavg.add(dateAvg);
						woc.LowA.add(ws.LowA);
						woc.MidA.add(ws.MidA);
						woc.HighA.add(ws.HighA);
						woc.FixA.add(ws.FixA);
						woc.TotalA.add(ws.Adapt);
						woc.numSamples.add(mat.length);
						woc.neut=nr;
						woc.L=matans2.length;
					}

					woc.BS.get(Z)[run]=w.Adapt;

				}else { // this is for everything else
					Williamson3bin ww = new Williamson3bin(mat,matans); //******
					Store s = ww.CreateBlocks(3,mat[0].length,sampler); //******
					Williamson3bin w = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
					w.williamson3bin_method(nr,low, mid, high);


					if(run==0 && Z==0){
						woc.Low_S.add(Double.valueOf(0.0));
						woc.Low_R.add(Double.valueOf(0.0));
						woc.Mid_S.add(Double.valueOf(0.0));
						woc.Mid_R.add(Double.valueOf(0.0));
						woc.High_S.add(Double.valueOf(0.0));
						woc.High_R.add(Double.valueOf(0.0));
						woc.Fix_S.add(Double.valueOf(0.0));
						woc.Fix_R.add(Double.valueOf(0.0));
						woc.dateavg.add(dateAvg);
						woc.LowA.add(Double.valueOf(0.0));
						woc.MidA.add(Double.valueOf(0.0));
						woc.HighA.add(Double.valueOf(0.0));
						woc.FixA.add(Double.valueOf(0.0));
						woc.TotalA.add(Double.valueOf(0.0));
						woc.numSamples.add(mat.length);
						woc.neut=nr;
						woc.L=matans.length;
					}
					else if(run==0){
						Williamson3bin ws = new Williamson3bin(mat,matans);
						ws.williamson3bin_method(nr,low, mid, high);
						woc.Low_S.add(ws.low_S);
						woc.Low_R.add(ws.low_R);
						woc.Mid_S.add(ws.Neut_S);
						woc.Mid_R.add(ws.Neut_R);
						woc.High_S.add(ws.High_S);
						woc.High_R.add(ws.High_R);
						woc.Fix_S.add(ws.Fix_S);
						woc.Fix_R.add(ws.Fix_R);
						woc.dateavg.add(dateAvg);
						woc.LowA.add(ws.LowA);
						woc.MidA.add(ws.MidA);
						woc.HighA.add(ws.HighA);
						woc.FixA.add(ws.FixA);
						woc.TotalA.add(ws.Adapt);
						woc.numSamples.add(mat.length);
						woc.neut=nr;
						woc.L=matans.length;
					}

					woc.BS.get(Z)[run]=w.Adapt;

				}

			}
		}


		for(int l=0;l<woc.BS.size();l++){
			Arrays.sort(woc.BS.get(l));
		}
		double lc= 0.05*Num; double hc= 0.95*Num; 
		for(int l=0;l<woc.BS.size();l++){				
			woc.LC.add(woc.BS.get(l)[(int)lc]);
			woc.HC.add(woc.BS.get(l)[(int)hc]);
		}
		return woc;
	}

	public WilliamsonOutputClass GenericRunner(String Location,String[] Dates,double nr,int Num){
		DataSet data = new DataSet(Location);
		System.out.println("Analysis on "+Location);
		//ancestral sequence
		WilliamsonOutputClass woc = new WilliamsonOutputClass(Dates.length,Num);
		final double[] low = {0.0,0.15};
		final double[] mid = {0.15,0.75};
		final double[] high = {0.75,1.0};	

		///////******************************************************* Generate Ans**************************
		ArrayList<int[]> Ans = new ArrayList<int[]>();
		String[] temp = Dates[0].split("_");
		for(int i=0;i<data.taxon_matrix.length;i++){
			String line = data.taxon_matrix[i];
			for(int x=0;x<temp.length;x++){
				if(line.matches(".*"+String.valueOf(temp[x]))){
					Ans.add(data.integer_matrix[i]);
				}
			}
		}
		int[][] ansMat = new int[Ans.size()][Ans.get(0).length];
		for(int i=1;i<Ans.size();i++){
			ansMat[i]=Ans.get(i);
		}
		Read_main r = new Read_main();
		int[] ans = r.consensusArray(ansMat); // this is the ancestral sequence
		///////******************************************************* Generate Main**************************
		// DO main sequence
		for(int run=0;run<Num;run++){

			for(int d=0;d<Dates.length;d++){	// cycle through all time points
				temp = Dates[d].split("_"); // split clumped dates
				double dateAvg=0; //initialise date counter
				double k=0;
				ArrayList<int[]> Seq = new ArrayList<int[]>();
				for(int i=0;i<data.taxon_matrix.length;i++){
					String line = data.taxon_matrix[i];
					for(int x=0;x<temp.length;x++){
						if(line.matches("(?i).*_"+String.valueOf(temp[x])+".*")){
							String[] split = line.split("_");
							dateAvg+=(double) Double.valueOf(split[split.length-1]);
							k++;
							Seq.add(data.integer_matrix[i]);
						}
					}
				}
				int[][] seq = new int[Seq.size()][Seq.get(0).length];
				for(int i=0;i<Seq.size();i++){
					seq[i]=Seq.get(i);
				}
				// update the teaspoon.adaptation.Williamson Output class
				if(run==0){
					Williamson3bin ws = new Williamson3bin(seq,ans);
					ws.williamson3bin_method(nr,low, mid, high);
					woc.Low_S.add(ws.low_S);
					woc.Low_R.add(ws.low_R);
					woc.Mid_S.add(ws.Neut_S);
					woc.Mid_R.add(ws.Neut_R);
					woc.High_S.add(ws.High_S);
					woc.High_R.add(ws.High_R);
					woc.Fix_S.add(ws.Fix_S);
					woc.Fix_R.add(ws.Fix_R);
					woc.dateavg.add(dateAvg/k);
					woc.LowA.add(ws.LowA);
					woc.MidA.add(ws.MidA);
					woc.HighA.add(ws.HighA);
					woc.FixA.add(ws.FixA);
					woc.TotalA.add(ws.Adapt);
					woc.numSamples.add(seq.length);
					woc.neut=nr;
					woc.L=ans.length;	
				}else{
					///////******************************************************* Generate Bootstraps**************************
					// generate a new sample set
					int[] sampler = new int[ans.length/3];			
					for(int x=0;x<sampler.length;x++){
						int randint = generator.nextInt(sampler.length-1);			
						sampler[x] = randint;
					}
					Williamson3bin ww = new Williamson3bin(seq,ans); //******
					Store s = ww.CreateBlocks(3,seq[0].length,sampler); //******
					Williamson3bin w = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
					w.williamson3bin_method(nr,low, mid, high);
					woc.BS.get(d)[run]=w.Adapt;
				}
			}		
		}
		return woc;

	}

	public void RunnerTS(String Location, String Table, String Gene,String[] type,String[] Dates,int start,int stop,double nr,int which){


		DataSet d = new DataSet(Location,Table,Gene);
		final double[] low = {0.0,0.15};
		final double[] mid = {0.15,0.75};
		final double[] high = {0.75,1.0};	

		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();

		/************************** Specify MAIN SEQUENCE************************************************************ */
		double dateAvg = 0; // calculate date average
		for(int Z=1;Z<Dates.length;Z++){  // BIG loop through dates ************************ change date in loops
			It =  d.DataMainFrame.iterator();
			double count=0;
			dateAvg = 0;
			ArrayList<int[]> Sequences = new ArrayList<int[]>(); 
			ArrayList<String> Names = new ArrayList<String>(); 
			ArrayList<SequenceInfo> out = new ArrayList<SequenceInfo>();

			ArrayList<int[]> Sequences2 = new ArrayList<int[]>(); 
			ArrayList<String> Names2 = new ArrayList<String>(); 
			ArrayList<SequenceInfo> out2 = new ArrayList<SequenceInfo>();
			while(It.hasNext()){		
				SequenceInfo element = It.next();
				// find main sequences
				boolean flag = false;  // first flag checks through genotypes
				for(int f=0;f<type.length;f++){
					if(element.Genotype.equals(type[f])){ // fiter genotypes
						flag=true;
					}
				}
				String delimiter = "-";
				String[] temp = Dates[Z].split(delimiter);
				boolean flag2 = false; // second flag checks through dates
				for(int n=0;n<temp.length;n++){
					if(element.Year==Double.valueOf(temp[n])){
						flag2=true;
					}
				}
				// adds the sequences we want according to date,genotype and gap uncertainty
				if(element.Gap<5.0 && flag==true && flag2==true){	// filter sequences with gap count more than five
					Sequences.add(element.Sequence); //saves sequence
					Names.add(element.Accession); //saves name
					out.add(element); //for output if needed
					dateAvg+=element.DecimalDate; //calculates average date
					count++;
				}
				// find ancestral sequences
				flag = false;  // first flag checks through genotypes
				for(int f=0;f<type.length;f++){
					if(element.Genotype.equals(type[f])){ // fiter genotypes
						flag=true;
					}
				}
				delimiter = "-";
				temp = Dates[Z-1].split(delimiter);
				flag2 = false; // second flag checks through dates
				for(int n=0;n<temp.length;n++){
					if(element.Year==Double.valueOf(temp[n])){
						flag2=true;
					}
				}
				// adds the sequences we want according to date,genotype and gap uncertainty
				if(element.Gap<5.0 && flag==true && flag2==true){	// filter sequences with gap count more than five
					Sequences2.add(element.Sequence); //saves sequence
				}

			}
			// make integer matricies dates ********************************************************************************			
			dateAvg = dateAvg/count;
			int[][] integer_matrix = new int[Sequences.size()][Sequences.get(0).length];
			String[] accessions = new String[Names.size()];
			// add sequences
			Iterator<int[]> It2 =  Sequences.iterator(); // converts object into a sequence matrix
			int l=0;
			while(It2.hasNext()){
				int[] tmp = It2.next();
				integer_matrix[l] = tmp ;
				l++;
			}
			Iterator<String> It3 =  Names.iterator(); //converts object into the accessions matrix
			l=0;
			while(It3.hasNext()){
				String tmp = It3.next();
				accessions[l] = tmp ;
				l++;
			}

			if(stop==-1){
				stop = integer_matrix[0].length;
			}

			// make ancestreal matricies dates ********************************************************************************	
			int[][] integer_ancestral = new int[Sequences2.size()][Sequences2.get(0).length];
			Iterator<int[]> It4 =  Sequences2.iterator(); // converts object into a sequence matrix
			l=0;
			while(It4.hasNext()){
				int[] tmp = It4.next();
				integer_ancestral[l] = tmp ;
				l++;
			}
			int[] ans = consensusArray(integer_ancestral);

			// adjust matrix size
			int[][] mat = new int[integer_matrix.length][stop-start];  //start=51, stop=1035 for ha1, ha2 is start 1036 - end
			int[] matans = new int[stop-start];

			for(int i=0;i<integer_matrix.length;i++){
				int p=0;
				for(int j=start;j<stop;j++){		
					mat[i][p] = integer_matrix[i][j];
					matans[p] = ans[j];
					p++;
				}
			}
			Williamson3bin ws = new Williamson3bin(mat,matans);


			ws.williamson3bin_method(nr,low, mid, high);

			if(which==1){
				System.out.print(dateAvg  + "\t" + ws.integer_matrix.length + "\t" + ws.low_R+"\t"+ws.low_S + "\t" + ws.Neut_R+"\t"+ws.Neut_S + "\t" + ws.High_R+"\t"+ws.High_S+"\n");
			} else {
				//		System.out.print(dateAvg + "\t" +  ws.Adapt + "\n" );
				System.out.print(dateAvg + "\t" +  ws.LowA+ "\t" +  ws.MidA+ "\t" +  ws.Adapt + "\n" );

			}



		}

	}


	public int[] consensusArray(int[][] integer_matrix){
		Methods preprocess = new Methods();
		int[] consensus = new int[integer_matrix[0].length];
		double[] counter = new double[5];
		for(int site=0;site<integer_matrix[0].length;site++){
			// count numbers
			counter[0] = preprocess.num_of_base(integer_matrix, 1, site);
			counter[1] = preprocess.num_of_base(integer_matrix, 2, site);
			counter[2] = preprocess.num_of_base(integer_matrix, 3, site);
			counter[3] = preprocess.num_of_base(integer_matrix, 4, site);
			counter[4] = preprocess.num_of_base(integer_matrix, 5, site);
			int length = counter.length;
			double max = -1;
			int position = 0;
			for (int i = 0; i < length; i++) {
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

	public void Extract(String Location, String Table, String Gene,String[] type,String Dates,String outname){
		DataSet d = new DataSet(Location,Table,Gene);
		Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();
		ArrayList<SequenceInfo> out = new ArrayList<SequenceInfo>();
		while(It.hasNext()){		
			SequenceInfo element = It.next();
			boolean flag = false;  // first flag checks through genotypes
			for(int f=0;f<type.length;f++){
				if(element.Genotype.equals(type[f])){ // fiter genotypes
					flag=true;
				}
			}
			String delimiter = "-";
			String[] temp = Dates.split(delimiter);
			boolean flag2 = false; // second flag checks through dates
			for(int n=0;n<temp.length;n++){
				if(element.Year==Double.valueOf(temp[n])){
					flag2=true;
				}
			}
			// adds the sequences we want according to date,genotype and gap uncertainty
			if(element.Gap<5.0 && flag==true && flag2==true){	// filter sequences with gap count more than five
				out.add(element);
			}		
		}		
		d.exportFASTA(outname, out);
	}

	public void print(WilliamsonOutputClass woc,int type,String outp){
		if(type==1){  // output adaptation and percentiles
			for(int i=0;i<woc.BS.size();i++){
				System.out.println(woc.dateavg.get(i)  + "\t" + woc.LC.get(i) + "\t" + woc.TotalA.get(i) + "\t" + woc.HC.get(i) );
			}
		}else if(type==2){ // out put silent and replace stats
			for(int i=0;i<woc.BS.size();i++){
				System.out.print(woc.dateavg.get(i)  + "\t" + woc.numSamples.get(i) + "\t" + woc.Low_R.get(i)+"\t"+woc.Low_S.get(i) + "\t" + woc.Mid_R.get(i)+"\t"+woc.Mid_S.get(i) + "\t" + woc.High_R.get(i)+"\t"+woc.High_S.get(i)+"\n");
			}

		}else if(type==3){ // output all adaptation from all ranges
			for(int i=0;i<woc.BS.size();i++){
				System.out.print(woc.dateavg.get(i) + "\t" +  woc.LowA.get(i)+ "\t" +  woc.MidA.get(i)+ "\t" +  woc.TotalA.get(i) + "\n" );
			}
		}else if(type==4){ // output just top adaptations
			for(int i=0;i<woc.BS.size();i++){
				System.out.print(woc.dateavg.get(i) + "\t" +  woc.TotalA.get(i) + "\n" );
			}
		} else if(type==5) {  // output silent and replacement from mid neutral bin
			for(int i=0;i<woc.BS.size();i++){
				System.out.print(woc.dateavg.get(i) + "\t" + woc.Mid_R.get(i)+"\t"+woc.Mid_S.get(i) + "\n" );
			}
		}else if(type==6) { // output bootstraps
			for(int j=0;j<woc.BS.get(0).length;j++){
				for(int i=0;i<woc.BS.size();i++){
					System.out.print(woc.BS.get(i)[j]+"\t");
				}
				System.out.println();
			}
		}else if(type==7) {  // write bootstraps
			try{
				FileWriter fstream = new FileWriter(outp+".txt");
				BufferedWriter out = new BufferedWriter(fstream);
				for(int j=0;j<woc.BS.get(0).length;j++){
					for(int i=0;i<woc.BS.size();i++){
						out.write(woc.BS.get(i)[j]+"\t");
					}
					out.write("\n");
				}
				out.close();

			}	catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}			

		} else if(type==8){  // stats
			Predict P = new Predict(woc.dateavg,woc.TotalA);
			double[] E = P.FitPoly(1);
			System.out.print(woc.neut + "\t"+ E[0]+"\t"+(E[0]/(woc.L/3.0)));
		}else if(type==9){  // stats
			try{
				FileWriter fstream = new FileWriter(outp+".txt");
				BufferedWriter out = new BufferedWriter(fstream);
				for(int i=0;i<woc.BS.size();i++){
					out.write(r(woc.dateavg.get(i)) + "\t\t " +  r(woc.TotalA.get(i)) + "\n" );
				}
				out.close();
			} 	catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}


		}else if(type==10){  // big output
			double[] E = Summarise(woc);
			double[] P = ModelSelection(woc);
			double[] ests = parameterEstimates(woc);
			double[] mse = MeanSquareErrors(woc);
			double[] prop = Proportion(woc);
			Date todaysDate = new java.util.Date();
			SimpleDateFormat formatter = new SimpleDateFormat("EEE, dd-MMM-yyyy HH:mm:ss");
			SimpleDateFormat formatter2 = new SimpleDateFormat("dd.MM.yy");
			String formattedDate = formatter.format(todaysDate);
			String formattedDate2 = formatter2.format(todaysDate);

			try{
				FileWriter fstream = new FileWriter(outp+"["+formattedDate2+"]"+".txt");
				BufferedWriter out = new BufferedWriter(fstream);
				// print basic data
				out.write("Results Date Stamp: "+formattedDate+"\n");
				out.write("The Sequence Length is "+woc.L+"\n");
				out.write("The number of bootstraps is: "+woc.BS.get(0).length+"\n");
				out.write("Neutral ratio: "+r(woc.neut) +"\n\n");
				out.write("The estimated rate for a 1st order linear model is: "+r(E[0])+"\n");
				out.write("The per codon estimated rate for a 1st order linear model is: "+(E[0]/(woc.L/3.0))+"\n");
				out.write("The Standard error is: "+r(E[3])+"\n");
				out.write("The 0.95 hi percentile is "+r(E[4])+"\n");
				out.write("The 0.05 low percentile is "+r(E[5])+"\n");
				out.write("The Per Codon 0.95 hi percentile is "+(E[4]/(woc.L/3.0))+"\n");
				out.write("The Per Codon  0.05 low percentile is "+(E[5]/(woc.L/3.0))+"\n");
				out.write("The Bias is: "+r(E[2])+"\n");
				out.write("The Mean is: "+r(E[1])+"\n");
				out.write("the proportion of adaptive fixations is: "+ String.valueOf(prop[0])+"\n");
				out.write("number of points used in proportion calculation: "+ String.valueOf(prop[1])+"\n");
				out.write("\n\n");
				// print the model data

				out.write("The mean square error for linear: bX model is: " +r(mse[0])+"\n");
				out.write("The mean square error for Quadratic: b1X^2 + bX model is: " +r(mse[1])+"\n");
				out.write("The mean square error for Cubic: b2X^3 + b1X^2 + bX model is: " +r(mse[2])+"\n\n");
				String lin = "linear: bX";
				String quad = "Quadratic: b1X^2 + bX";
				String cube = "Cubic: b2X^3 + b1X^2 + bX";
				String which="";
				if(P[0]==1){which=lin;}if(P[0]==2){which=quad;}if(P[0]==3){which=cube;}

				out.write("The best model for the actual data is: "+which+"\n");
				out.write("The parameter estimates are: \n");
				String[] coeff= {"b","b1","b2"};
				for(int i=0;i<ests.length;i++){
					out.write(coeff[i]+"="+r(ests[i])+"\n");
				}
				out.write("From the bootstraps the proportion which have bX as the best model is: "+ P[1]+"\n");
				out.write("From the bootstraps the proportion which have b1X^2 + bX as the best model is: "+ P[2]+"\n");
				out.write("From the bootstraps the proportion which have b2X^3 + b1X^2 + bX as the best model is: "+ P[3]+"\n");

				out.write("\n\n");
				out.write("----------------Data calculated Values---------------");
				out.write("\n\n");
				out.write("Adaptive substitutions per year \n");
				out.write("data\t\tNumber of adaptive fixations\n");
				for(int i=0;i<woc.BS.size();i++){
					out.write(r(woc.dateavg.get(i)) + "\t\t " +  r(woc.TotalA.get(i)) + "\n" );
				}
				out.write("\n\n");
				out.write("date\t\tNr\t\tNs\n");
				for(int i=0;i<woc.BS.size();i++){
					out.write(r(woc.dateavg.get(i)) + "\t\t" + r(woc.Mid_R.get(i))+"\t\t"+r(woc.Mid_S.get(i)) + "\n" );
				}
				out.write("\n\n");
				out.write("data\t\tNumber of Samples\n");
				for(int i=0;i<woc.BS.size();i++){
					out.write(r(woc.dateavg.get(i)) + "\t\t" + woc.numSamples.get(i)+"\n");
				}

				Predict Pest = new Predict(woc.dateavg,woc.TotalA,woc.BS);
				double [] estimates = Pest.FitLinear();
				out.write("\n\n");
		//		out.write("-------------These are the bootstrap regression values for a linear model--------------\n");
		//		for(double ez:estimates){
		//			out.write(ez+"\n");
		//		}

				out.close();
			} 	catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}	

		} else if(type==11){
			try{
				FileWriter fstream = new FileWriter(outp+".txt");
				BufferedWriter out = new BufferedWriter(fstream);
				// print actual results
				for(int i=0;i<woc.BS.size();i++){
					out.write(woc.dateavg.get(i) + "\t");
				}
				out.write("\n");
				for(int i=0;i<woc.BS.size();i++){
					out.write(woc.TotalA.get(i) + "\t" );
				}
				out.write("\n");
				// pirnt bootstraps
				for(int j=0;j<woc.BS.get(0).length;j++){
					for(int i=0;i<woc.BS.size();i++){
						out.write(woc.BS.get(i)[j]+"\t");
					}
					out.write("\n");
				}
				out.close();
			}	catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}	


		} else if(type==12){
			try{
				FileWriter fstream = new FileWriter(outp+".txt");
				BufferedWriter out = new BufferedWriter(fstream);
				// print actual results
				Predict P = new Predict(woc.dateavg,woc.TotalA,woc.BS);
				double[][] mses = P.bestModelDist();

				for(int i=0;i<mses.length;i++){
					out.write(mses[i][0] +"\t" + mses[i][1] +"\t" +mses[i][2] +"\n");
				}

				out.close();
			}	catch (Exception e){//Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}	


		}
		System.out.println();
	}
	// output rounding function
	public String r(double val){
		DecimalFormat df = new DecimalFormat("#.####");
		return df.format(val);


	}

	public double CalcNr(WilliamsonOutputClass woc){
		double S=0; double R=0;

		for(int i=1;i<woc.Mid_S.size();i++){
			S+=woc.Mid_S.get(i);
			R+=woc.Mid_R.get(i);
		}
		double neutralratio=R/S;
		return neutralratio;
	}

	public double[] Statistics(double actual,double[] data){

		// sd is sqrt of sum of (values-mean) squared divided by n - 1
		// Calculate the mean
		double mean = 0;
		final int n = data.length;
		for ( int i=0; i<n; i++ )
		{
			mean += data[i];
		}
		mean /= n;
		double bias = mean-actual;
		// calculate the sum of squares
		double sum = 0;
		for ( int i=0; i<n; i++ )
		{
			final double v = data[i] - mean;
			sum += v * v;
		}
		// Change to ( n - 1 ) to n if you have complete data instead of a sample.
		double std = Math.sqrt( sum / ( n - 1 ) );
		Arrays.sort(data);

		double h = 0.95*data.length; double l = 0.05*data.length; double m = 0.5*data.length;
		double hi = data[(int) h+1];
		double low = data[(int)l+1 ];
		double median = data[(int)m+1 ];

		double ests[] = {actual,mean,bias,std,hi,low,median};
		return ests;	
	}

	public double[] Summarise(WilliamsonOutputClass woc){

		Predict P = new Predict(woc.dateavg,woc.TotalA,woc.BS);
		double [] estimates = P.FitLinear();
		double [] actual = P.FitPoly(1);
		double[] stats = Statistics(actual[0],estimates);
		return stats;	

	}

	public double[] Proportion(WilliamsonOutputClass woc){
		double Avg = 0;
		double count=0;
		for(int i=1;i<woc.TotalA.size();i++){
			if(woc.High_R.get(i)+woc.Fix_R.get(i) != 0)	{
				Avg+=(woc.HighA.get(i)+woc.FixA.get(i))/(woc.High_R.get(i)+woc.Fix_R.get(i));
				count++;
			}		
		}
		Avg = Avg/count;
		double[] result = {Avg, count};
		return result;
	}

	public double[] ModelSelection(WilliamsonOutputClass woc){
		Predict P = new Predict(woc.dateavg,woc.TotalA,woc.BS);
		double [] p = P.bestModel();
		double  a = P.WhichModel();
		double[] vals = {a,p[0],p[1],p[2]};
		return vals;			
	}

	public double[] parameterEstimates(WilliamsonOutputClass woc){
		Predict P = new Predict(woc.dateavg,woc.TotalA,woc.BS);
		double w = P.WhichModel();
		double[] ests = P.FitPoly((int) w);
		return ests;
	}

	public double[] MeanSquareErrors(WilliamsonOutputClass woc){
		Predict P = new Predict(woc.dateavg,woc.TotalA,woc.BS);
		double[] mse = P.MSEval();
		return mse;
	}

	public void ExtractVal(String[] List,int line, int word,String loc){
		for(String x:List){
			if(new File(loc+x).exists()){
				try{ //read one text file
					FileReader fr = new FileReader (loc+x);
					BufferedReader br = new BufferedReader(fr);
					Scanner s = new Scanner(br); 
					int k=0;
					while (s.hasNextLine()) {
						String tmp = s.nextLine();
						if(k==line){
							String[] temp;
							String delimiter = " ";
							temp = tmp.split(delimiter);
							System.out.print(temp[word] + "\t" );					 
						}
						k++;

					}

					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			} else {
				System.out.print("\t");
			}
		}

		System.out.println();	
	}

}
