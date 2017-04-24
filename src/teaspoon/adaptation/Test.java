package teaspoon.adaptation;

import java.util.ArrayList;
import java.util.Iterator;

//import weka.core.ContingencyTables;


//java -classpath /Users/sam/Desktop/ECLIPSE\ FILES/Adaptive\ Substitition/ teaspoon.adaptation.Test
//XMLOutputter serializer = new XMLOutputter(Format.getRawFormat().setIndent(" ").setLineSeparator("\n"));
//import java.io.File;
//System.out.println(System.getProperty("java.class.path"));

public class Test {

	public static void main(String[] args)  {



		//		System.out.println(System.getProperty("java.class.path"));
		/*    System.out.println("Starting test");
	    long time = System.currentTimeMillis();
//		File input = new File("/Users/sam/Desktop/VirusDataSets/HCVseqs/seq/");		// Input 
//		File input = new File("/Users/sam/Desktop/CleanedUPFINAL/");		// Input 
			File input = new File("/Users/sam/Desktop/VirusDataSets/TestData/50/");		// Input 
		File[] list = input.listFiles();
		ArrayList<String> data = new ArrayList<String>();
		for(int i=0;i<list.length;i++){
			if(list[i].isHidden() == false && list[i].isDirectory() == false)	{
				data.add(list[i].getAbsolutePath());
			}
		}	
		Collections.sort(data);
		double[][] matrix = new double[data.size()][10];
		for(int x=0;x<data.size();x++){
	//		teaspoon.adaptation.Read_main re = new teaspoon.adaptation.Read_main("/Users/sam/Desktop/VirusDataSets/HCVseqs/ans/ans.phy");
			teaspoon.adaptation.Read_main re = new teaspoon.adaptation.Read_main("/Users/sam/Desktop/VirusDataSets/TestData/ans50.nex");
			teaspoon.adaptation.Read_main ra = new teaspoon.adaptation.Read_main(data.get(x));

			int[][] seq = ra.readNEXUS();
			int[][] tmp = re.readNEXUS();



			int[] ans = re.consensusArray(tmp);	
			teaspoon.adaptation.SiteEstMulti sm = new teaspoon.adaptation.SiteEstMulti(seq,ans);
			double[] L = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
			double[] H = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
//			double[] L = {0.0};
//			double[] H = {1.0};
			boolean[] Nvec = {false,false,false,false,true,true,false,false,false,false};

			double[][] bins = new double[2][L.length];
			for(int i=0;i<L.length;i++){
				bins[0][i]=L[i];
				bins[1][i]=H[i];
			}

		//	System.out.println(seq.length);
		//	System.out.println(data.get(x).toString());
		//	System.out.println(sm.howManyInvariant());


			double[] prior = {1,1,1,1};
			double[][] val = sm.Value_Matrix(Nvec, bins,prior,true);
		//	System.out.println(sm.howManySingPoly());

	 //		System.out.println(data.get(x).toString());
	//		System.out.print("	");

	//			System.out.print(sm.Distance());
	//			System.out.print("	");
			for(int i=0;i<val[2].length;i++){
	//			matrix[x][i]=val[3][i];
				System.out.print(val[2][i]);
				System.out.print("	");
			}
			System.out.println();

			//	for(int i=0;i<val.length;i++){
			//		for(int j=0;j<val[0].length;j++){
			//			System.out.print(val[i][j]);
			//		System.out.print("	");
			//		}
		//			System.out.println("");
		//		}


		} 
			File file = new File("/Users/sam/Desktop/HCV.txt");
	      try{
	    	    // Create file 
	    	    FileWriter fstream = new FileWriter(file);
	    	        BufferedWriter out = new BufferedWriter(fstream);
	    			for(int i=0;i<matrix.length;i++){
	    				for(int j=0;j<matrix[0].length;j++){
			    						out.write(new Double(matrix[i][j]).toString());
	    						out.write("	");
	    					}
	    					out.write("\n");
	    			}	    			

	    	    //Close the output stream
	    	    out.close();
	    	    }catch (Exception e){//Catch exception if any
	    	      System.err.println("Error: " + e.getMessage());
	    	    }
	    time = System.currentTimeMillis() - time;
	    System.out.println(" The test took " + time + " milliseconds");*/

		//		*********************
		//		*********************
		//		*********************		
		//		String sss = new String("/Users/sam/Desktop/NeutSeq/NeutLow/");

		// START ****************************************************************************************************************************************
		/*			double[] store = new double[3];
		double lowbin=0;
		double midbin=0;
		double highbin=0;

		double[] Tstore = new double[3];
		double Tlowbin=0;
		double Tmidbin=0;
		double Thighbin=0;
					String sss = new String("/Users/sam/Desktop/VirusDataSets/96VD/");
	//	String sss = new String("/Users/sam/Desktop/tmp/");
//		File input = new File(sss+"seqNOSTOP");		// Input 
		File input = new File(sss+"eseq");		// Input 
		File[] list = input.listFiles();
		ArrayList<String> data = new ArrayList<String>();
		ArrayList<String> names = new ArrayList<String>();
		for(int i=0;i<list.length;i++){
			if(list[i].isHidden() == false && list[i].isDirectory() == false)	{
				data.add(list[i].getAbsolutePath());
				names.add(list[i].getName());
			}
		}	
		Collections.sort(data);
		Collections.sort(names);

		// highlight for HCV data
//		File input2 = new File(sss+"ansNOSTOP");		// Input 
		File input2 = new File(sss+"eans");		// Input 
		File[] list2 = input2.listFiles();
		ArrayList<String> data2 = new ArrayList<String>();
		for(int i=0;i<list2.length;i++){
			if(list2[i].isHidden() == false && list2[i].isDirectory() == false)	{
				data2.add(list2[i].getAbsolutePath());
			}
		}	
		Collections.sort(data2);

		for(int x=0;x<data.size();x++){

			teaspoon.adaptation.Read_main ra = new teaspoon.adaptation.Read_main(data.get(x));
			//	System.out.println(data.get(x));
			//	System.out.println(data2.get(x));
			teaspoon.adaptation.Read_main re = new teaspoon.adaptation.Read_main(data2.get(x));


			int[][] tmp = re.read();
			int[] ans = re.consensusArray(tmp);	
			int[][] seq = ra.read();

			teaspoon.adaptation.Methods m = new teaspoon.adaptation.Methods();

			//		double[] L = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
			//		double[] H = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
					double[] L = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
					double[] H = {0.1,0.2,0.3,0.4,0.5,0.6,0.7};
			boolean[] Nvec = {false,false,false,false,true,true,false,false,false,false};

	//		double[] L = {0.0,0.33,0.66};
	//		double[] H = {0.33,0.66,1.0};

	//		double[] L = {0.0,0.2};
	//		double[] H = {0.2,1.0};

			double[][] bins = new double[2][L.length];
			for(int i=0;i<L.length;i++){
				bins[0][i]=L[i];
				bins[1][i]=H[i];
			}
			double[] prior = {1,1,1,1};


			teaspoon.adaptation.BhattMethod bm = new teaspoon.adaptation.BhattMethod(seq,ans);
		//	double[][] bc = sm.BaseComposition();

		//	for(int i=0;i<bc.length;i++){
		//		for(int j=0;j<bc[0].length;j++){
		//			System.out.print(bc[i][j]);System.out.print(" ");
		//		}
		//		System.out.println();
		//	}

		//	double[][] val = sm.Value_Matrix(Nvec, bins,prior,true);

//			for(int i=0;i<val[2].length;i++){
//				System.out.print(val[3][i]);
//				System.out.print("	");
//			}

//			System.out.println();
			bm.Method(bins,prior,true);
			bm.print(bm.ReplacementSilentRatio);

	//		lowbin+=val[3][0];
	//		midbin+=val[3][1];
	//		highbin+=val[3][2];

	//		Tlowbin+=val[2][0];
	//		Tmidbin+=val[2][1];
	//		Thighbin+=val[2][2];


		} 
		//end ************************************************************************************************************************
		 */
		/*	store[0]=lowbin/data.size();
		store[1]=midbin/data.size();
		store[2]=highbin/data.size();

		Tstore[0]=Tlowbin/data.size();
		Tstore[1]=Tmidbin/data.size();
		Tstore[2]=Thighbin/data.size();
		try{
			BufferedWriter writer = new BufferedWriter(	new FileWriter("/Users/sam/Desktop/OutPut.txt",true)) ;
			//			writer.newLine();
			writer.write(String.valueOf(store[0])) ;
			writer.write("\t") ;
			writer.write(String.valueOf(store[1])) ;
			writer.write("\t") ;
			writer.write(String.valueOf(store[2])) ;
			writer.write("\r") ;
			writer.close() ;
		} catch (IOException ex) {
			ex.printStackTrace();
		}
		try{
			BufferedWriter writer = new BufferedWriter(	new FileWriter("/Users/sam/Desktop/OutPutT.txt",true)) ;
			//			writer.newLine();
			writer.write(String.valueOf(Tstore[0])) ;
			writer.write("\t") ;
			writer.write(String.valueOf(Tstore[1])) ;
			writer.write("\t") ;
			writer.write(String.valueOf(Tstore[2])) ;
			writer.write("\r") ;
			writer.close() ;
		} catch (IOException ex) {
			ex.printStackTrace();
		}*/


		/*			File file = new File("/Users/sam/Desktop/sni2.txt");
	      try{
	    	    // Create file 
	    	    FileWriter fstream = new FileWriter(file);
	    	        BufferedWriter out = new BufferedWriter(fstream);
	    			for(int i=0;i<matrix.length;i++){
	    				for(int j=0;j<matrix[0].length;j++){
			    						out.write(new Double(matrix[i][j]).toString());
	    						out.write("	");
	    					}
	    					out.write("\n");
	    			}	    			

	    	    //Close the output stream
	    	    out.close();
	    	    }catch (Exception e){//Catch exception if any
	    	      System.err.println("Error: " + e.getMessage());
	    	    }*/



		//		**********************************************************************************************************************	

		//		*********************
		//		*********************
		//		*********************

		//		System.out.println(sm.totalNoAdapt(Nvec, bins, 0.131));

		//	teaspoon.adaptation.Williamson3bin ww = new teaspoon.adaptation.Williamson3bin(seq,ans);
		//	System.out.printlhhhn(ww.williamson3bin_method(0.131));

		/*	int N=1000;
		teaspoon.adaptation.FluAnalysis ff = new teaspoon.adaptation.FluAnalysis();
for(int x=0;x<25;x++){
	//	teaspoon.adaptation.Value[][] v = ff.clumpanalysish3n2BS(N);
		teaspoon.adaptation.Value[][] v = ff.clumpanalysish3n2(x);
	//	teaspoon.adaptation.Value[][] v = ff.Swine();
		//	teaspoon.adaptation.Value[][] v = ff.HAanalysisH3N2(N);
		//	System.out.println("---------------------------------");
		//	System.out.println(" 					Results");
		//	System.out.println("---------------------------------");
		for(int i=0;i<v.length;i++){
			System.out.print(v[i][0].row);
			System.out.print("\t");
			for(int j=0;j<v[0].length;j++){

				System.out.print(v[i][j].Nadapt);

				System.out.print("	");
			}
			System.out.println("");
		}
}*/

		/*		File file = new File("/Users/sam/Desktop/30.txt");
		// Send formatted output to the file.
	      try{
	    	    // Create file 
	    	    FileWriter fstream = new FileWriter(file);
	    	        BufferedWriter out = new BufferedWriter(fstream);
	    			for(int x=0;x<N;x++){
	    				for(int i=0;i<v.length;i++){
	    					for(int j=0;j<v[0].length;j++){
	    						out.write(new Double(v[i][j].Bstrap[x]).toString());

	    						out.write("	");
	    					}
	    					out.write("\n");
	    				}
	    			}

	    	    //Close the output stream
	    	    out.close();
	    	    }catch (Exception e){//Catch exception if any
	    	      System.err.println("Error: " + e.getMessage());
	    	    }*/




		/*		teaspoon.adaptation.Read_main ra = new teaspoon.adaptation.Read_main("/Users/sam/Desktop/1977_1979.NA.H1N1.nex");
		teaspoon.adaptation.Read_main re = new teaspoon.adaptation.Read_main("/Users/sam/Desktop/2007_2009.NA.H1N1.nex") ;
		int[][] integer_matrix = ra.readNEXUS();
		int[][] tmp = re.readNEXUS();
		int[] integer_ancestral = re.consensusArray(tmp);	
		teaspoon.adaptation.Methods m = new teaspoon.adaptation.Methods();
		teaspoon.adaptation.Neuraminidase Nd = new teaspoon.adaptation.Neuraminidase();
		int[][] sub = m.Subsetter(integer_matrix, Nd.badlisth1n1, 2);
		System.out.print(sub[0].length);
		for(int i=0;i<sub.length;i++){
			for(int j=0;j<sub[0].length;j++){
	//			System.out.print(sub[i][j]);
			}
		}*/
		/*				teaspoon.adaptation.Read_main ra = new teaspoon.adaptation.Read_main("/Users/sam/Desktop/ans");
		teaspoon.adaptation.Read_main re = new teaspoon.adaptation.Read_main("/Users/sam/Desktop/seq");

		int[][] seq = re.readNEXUS();
		int[][] tmp = ra.readNEXUS();
		int[] ans = re.consensusArray(tmp);	

		teaspoon.adaptation.SiteEstMulti sm = new teaspoon.adaptation.SiteEstMulti(seq,ans);
		boolean[] Nvec = {false,false,false,false,true,true,false,false,false,false};
		double[] L = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
		double[] H = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
		double[][] bins = new double[2][L.length];
		for(int i=0;i<L.length;i++){
			bins[0][i]=L[i];
			bins[1][i]=H[i];
		}

		double[] prior = {1,1,1,1};
	    System.out.println("Starting test");
	    long time = System.currentTimeMillis();
	//	double[] id = sm.SiteFreq(0.0, 1.0, prior, true);

		double[][] val = sm.Value_Matrix(Nvec, bins,prior,true);



		for(int i=0;i<val[2].length;i++){
					System.out.print(val[3][i]);
					System.out.print("	");
		}
			System.out.println();*/

		/*		teaspoon.adaptation.Read_main ra = new teaspoon.adaptation.Read_main("/Users/sam/Desktop/a");
		teaspoon.adaptation.Read_main re = new teaspoon.adaptation.Read_main("/Users/sam/Desktop/s");

		int[][] seq = re.readNEXUS();
		int[][] tmp = ra.readNEXUS();
		int[] ans = re.consensusArray(tmp);	

		teaspoon.adaptation.SiteEstMulti sm = new teaspoon.adaptation.SiteEstMulti(seq,ans);

		double[] cf = sm.EmpiricalCodonFreq();

		for(int i=0;i<cf.length;i++){
			System.out.println(cf[i]);
		}*/






		/*	// extrat dates from ever year for ESPript
		teaspoon.adaptation.DataSet d = new teaspoon.adaptation.DataSet("/Users/sam/Desktop/sw_ref_H1.fasta","/Users/sam/Desktop/Table.csv","HA");
	//	Iterator<teaspoon.adaptation.SequenceInfo> It =  d.DataMainFrame.iterator();
		ArrayList<teaspoon.adaptation.SequenceInfo> ss = new ArrayList<teaspoon.adaptation.SequenceInfo>();

	int flag = 0;
	for(int i=1930;i<2012;i++){
		flag=0;
		Iterator<teaspoon.adaptation.SequenceInfo> It =  d.DataMainFrame.iterator();

		while(It.hasNext() && flag==0){
			teaspoon.adaptation.SequenceInfo element = It.next();
			if(element.Year==Double.valueOf(i) && element.Genotype.equals("#4C1")){
			//	System.out.println(element.Year);
				ss.add(element);
				flag=1;
				break;
			}
		}
	}
	d.exportFASTA("/Users/sam/Desktop/al.fasta", ss);	
		 */	
		// count number of dates ********************************************************************************
		//// start 
		/*		teaspoon.adaptation.DataSet d = new teaspoon.adaptation.DataSet("/Users/sam/Desktop/sw_ref_H1.fasta","/Users/sam/Desktop/Table_AvianEAans.csv","HA");
		Iterator<teaspoon.adaptation.SequenceInfo> It =  d.DataMainFrame.iterator();
		double x=0;
		int[] dates = new int[2011-1930];
		while(It.hasNext()){
			teaspoon.adaptation.SequenceInfo element = It.next();
			x=element.Year-1931;
			if(element.Genotype.equals("#4C1")){
				dates[(int)x]++;
			}
		}

		x=1931;
		for(int i=0;i<dates.length;i++){
			//		System.out.print(x+"\t");
			//		System.out.print(dates[i]);
			//		System.out.println();
			x++;
		}
		// define bounds ************************

		final double[] low = {0.0,0.15};
		final double[] mid = {0.15,0.75};
		final double[] high = {0.75,1.0};	


		// define set of Sequences with too many gaps ********************
	//	String[] badSeqs = {"AB573421","AB574412","AB574404","U47306","U47304","U45452"};  /// bad Sequences for GENOTYPE C1
		String[] badSeqs = {"AJ344013","GYSWB0444"};  /// bad Sequences for GENOTYPE B1
		String genotype = "#4B1";
		// Specify ANCESTRAL SEQUENCE********************************************************************************
	//	String[] clumpdates = {"1966","1066","1066","1975","1976","1066","1977","1066","1066","1978","1066","1066","1979","1066","1066","1980","1981","1066","1986","1987","1988","1990","1991","1066","1993","1066","1066","1994","1066","1066","2000","1066","1066","2001","1066","1066","2002","2003","1066","2004","1066","1066","2005","2006","1066","2008","2009","2010"};
		String[] clumpdates = {"1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","1066","2005","1066","1066","2006","1066","1066","2007","2008","1066","2009","2010","1066"};

		It =  d.DataMainFrame.iterator();
		int[] ans = new int[d.DataMainFrame.get(0).Sequence.length];
		ArrayList<teaspoon.adaptation.SequenceInfo> outans = new ArrayList<teaspoon.adaptation.SequenceInfo>();
		while(It.hasNext()){		
			teaspoon.adaptation.SequenceInfo element = It.next();
			if(element.Genotype.equals(genotype)){
			//	if(element.Year==1966){  // C1
				if(element.Accession.equals("AF091313")){ //B1 avian ancestral
					outans.add(element);
					ans = element.Sequence;
				}
			}
		}

		// Specify MAIN SEQUENCE********************************************************************************
		double dateAvg = 0; // calculate date average
		for(int Z=3;Z<clumpdates.length-2;Z=Z+3){  // BIG loop through dates ************************ change date in loops
			It =  d.DataMainFrame.iterator();
			double count=0;
			dateAvg = 0;
			ArrayList<int[]> Sequences = new ArrayList<int[]>(); 
			ArrayList<String> Names = new ArrayList<String>(); 
			ArrayList<teaspoon.adaptation.SequenceInfo> out = new ArrayList<teaspoon.adaptation.SequenceInfo>();
			while(It.hasNext()){		
				teaspoon.adaptation.SequenceInfo element = It.next();
				if(element.Genotype.equals(genotype)){ // filter genotypes
					if(element.Year==Double.valueOf(clumpdates[Z]) || element.Year==Double.valueOf(clumpdates[Z+1]) || element.Year==Double.valueOf(clumpdates[Z+2])){ // group dates
						boolean flag=true;
						for(int y=0;y<badSeqs.length;y++){ //loop to check for bad sequences
							if(element.Accession.equals(badSeqs[y])){flag=false;}							 					
						}
						if(flag){
							Sequences.add(element.Sequence); //saves sequence
							Names.add(element.Accession); //saves name
							out.add(element); //for output if needed
							dateAvg+=element.DecimalDate; //caluclates average date
							count++;
						}
					}
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

			// make HA1 ********************************************************************************
			int[][] HA1mat = new int[integer_matrix.length][984];
			int[] HA1matans = new int[984];

			for(int i=0;i<integer_matrix.length;i++){
				int p=0;
				for(int j=51;j<1035;j++){		
					HA1mat[i][p] = integer_matrix[i][j];
					HA1matans[p] = ans[j];
					p++;
				}
			}

			// make HA2 ********************************************************************************
			int[][] HA2mat = new int[integer_matrix.length][663];
			int[] HA2matans = new int[663];

			for(int i=0;i<integer_matrix.length;i++){
				int p=0;
				for(int j=1035;j<integer_matrix[0].length;j++){
					HA2mat[i][p] = integer_matrix[i][j];
					HA2matans[p] = ans[j];
					p++;
				}
			}


			double nr = 0.238035264;
			//run williamson  ********************************************************************************
	//		teaspoon.adaptation.Williamson3bin ws = new teaspoon.adaptation.Williamson3bin(HA2mat,HA2matans);

	//		ws.williamson3bin_method(nr, low, mid, high);

	//		teaspoon.adaptation.Williamson3bin ws2 = new teaspoon.adaptation.Williamson3bin(HA1mat,HA1matans);

	//		ws2.williamson3bin_method(nr, low, mid, high);

	//		teaspoon.adaptation.Williamson3bin ws3 = new teaspoon.adaptation.Williamson3bin(integer_matrix,ans);

	//		ws3.williamson3bin_method(nr, low, mid, high);


		//	double[] L = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
		//	double[] H = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

		double[] L=	{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95};
		double[] H ={0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1};

//	boolean[] Nvec = {false,true,true,true,true,true,true,false,false,false};
	boolean[] Nvec = {false,true,true,true,true,true,true,false,false,false,false,true,true,true,true,true,true,false,false,false};
	double[][] bins = new double[2][L.length];
	for(int i=0;i<L.length;i++){
		bins[0][i]=L[i];
		bins[1][i]=H[i];
	}
	double[] prior = {1,1,1,1};
	teaspoon.adaptation.BhattMethod bm = new teaspoon.adaptation.BhattMethod(integer_matrix,ans);
	teaspoon.adaptation.BhattMethod bm1 = new teaspoon.adaptation.BhattMethod(HA1mat,HA1matans);
	teaspoon.adaptation.BhattMethod bm2 = new teaspoon.adaptation.BhattMethod(HA2mat,HA2matans);

//	bm.Method(bins,prior,true,Nvec);
//	bm.print(bm.ReplacementSilentRatio);
//	System.out.println(bm.neutralratio + "\t"+ bm.DeleteriousLoad +"\t"+ bm.Adaptation);

	bm1.Method(bins,prior,true,Nvec,nr);
	bm1.print(bm1.NonNeutralSubstitutions);
//	System.out.println(bm1.neutralratio + "\t"+ bm1.DeleteriousLoad +"\t"+ bm1.Adaptation);

//	bm2.Method(bins,prior,true,Nvec);
//	bm2.print(bm2.ReplacementSilentRatio);
//	System.out.println(bm2.neutralratio + "\t"+ bm2.DeleteriousLoad +"\t"+ bm2.Adaptation);

//	System.out.println();

			//check gapinfo********************************************************************************
		//	ws.gapInfo();
		//	ws2.gapInfo();
		//	ws3.gapInfo();
		//	for(int i=0;i<ws.gapcount.length;i++){
		//		System.out.println(accessions[i]+ "\t" + dateAvg + "\t" + ws.gapcount[i] + "\t" + ws2.gapcount[i] + "\t" + ws3.gapcount[i] + "\t");
		//	}

			//output  ********************************************************************************

		//		System.out.print(ws.integer_matrix.length + "\t" + ws.mid_R+"\t"+ws.mid_S + "\t" + ws2.mid_R+"\t"+ws2.mid_S + "\t" + ws3.mid_R+"\t"+ws3.mid_S+"\n");
		//		System.out.print(dateAvg + "\t" + ws.Adapt + "\t" + ws2.Adapt + "\t" + ws3.Adapt + "\n" );	


		}
		 */

		//FluAnalysis2 F = new FluAnalysis2();

		//********************/********************/********************/********************/********************/********************/********************/********************
		/*		// this bit of code gives the gap mainAnalysis
		String[] loc = {"sw_ref_PB2.fasta","sw_ref_PB1.fasta","sw_ref_PA.fasta","sw_ref_H1.fasta","sw_ref_NP.fasta","sw_ref_N1.fasta","sw_ref_MP.fasta","sw_ref_NS.fasta"};
		String[] Gene = {"PB2","PB1","PA","HA","NP","NA","MP","NS"};
		double[][] mat = new double[8][ts];
		String[][] mat2 = new String[8][ts];
		// this prints out the gap percentage
//		for(int i=0;i<loc.length;i++){
//			double[] tmp =F.GAPcounter("/Users/sam/Desktop/FluData/"+loc[i],"/Users/sam/Desktop/FluData/Table.csv",Gene[i], ts) ;
//			mat[i] = tmp;
//		}

		// prints the accession numbers 
		for(int i=0;i<loc.length;i++){
			String[] tmp =F.Acc("/Users/sam/Desktop/FluData/"+loc[i],"/Users/sam/Desktop/FluData/Table.csv",Gene[i], ts) ;
			mat2[i] = tmp;
		}

		for(int i=0;i<mat2[0].length;i++){
			for(int j=0;j<mat2.length;j++){
				System.out.print(mat2[j][i]+"\t");
			}
			System.out.println();
		}*/
		//********************/********************/********************/********************/********************/********************/********************/********************

		//********************/********************/********************/********************/********************/********************/********************/********************
		// this bit of code produces the date counts
		//	String[] Types = {"#6C1","#6C1a","#6C1b","#6C1c","#6C1d","#6C1e"};
		//	all
	//		String[] loc = {"sw_ref_PB2.fasta","sw_ref_PB1.fasta","sw_ref_PA.fasta","sw_ref_H1.fasta","sw_ref_NP.fasta","sw_ref_N1.fasta","sw_ref_M1.fasta","sw_ref_M2.fasta","sw_ref_NS1.fasta","sw_ref_NS2.fasta"};
	//	String[] Gene = {"PB2","PB1","PA","HA","NP","NA","MP","MP","NS","NS"};
	/*	String[] loc = {"sw_ref_H1.fasta"};
		String[] Gene = {"HA"};
	
		
		int[][] mat = new int[10][2011-1930];
	//	String[] Types = {"#1B1","#2B1","#3B1","#4B1","#5B1","#6B1","#7B1","#7B1","#8B1","#8B1"};
	//	String[] Types = {"#1T1","#2T1","#3T1","#4T1","#5T1","#6T1","#7T1","#7T1","#8T1","#8T1"};
	//	String[] Types = {"#1C1","#2C1","#3C1","#4C1","#5C1","#6C1","#7C1","#7C1","#8C1","#8C1"};
		String[] Types = {"#4B1"};
		for(int i=0;i<loc.length;i++){
			mat[i]= F.Datecounter2("/Users/sam/Desktop/FluData/"+loc[i],"/Users/sam/Desktop/FluData/TableT.csv",Gene[i],Types[i]);
		}

		for(int i=0;i<mat[0].length;i++){
			for(int j=0;j<mat.length;j++){
				System.out.print(mat[j][i]+"\t");
			}
			System.out.println();
		}
*/

		/*	String Gene = "NS1";
		int[] mat = new int[2011-1930];

		int L=8;
		String[] typeout = new String[18];
		String[] subs = {"","a","b","c","d","e","f","g","h",};
		String[] type={"C1","T1"};
		int k=0;
		for(int j=0;j<type.length;j++){
			for(int i=0;i<subs.length;i++){
				typeout[k] = "#"+String.valueOf(L)+type[j]+subs[i];
				k++;
			}
		}


			mat = F.Datecounter("/Users/sam/Desktop/FluData/sw_ref_"+Gene+".fasta","/Users/sam/Desktop/FluData/Table.csv",Gene,typeout);



			for(int j=0;j<mat.length;j++){
				System.out.println(mat[j]);
			}

		 */


		/* Exporter ***************************************************************************************************/
							int[] array = {0,1,2,3,4,5,6,6,7,7};
		for(int x=0;x<array.length;x++){
			int l=array[x];
			String[] typeout = new String[18];
			//	String[] typeout = new String[9];
		//
			String[] subs = {"","a","b","c","d","e","f","g","h",};
			String[] type={"T1"};
			int k=0;
			for(int j=0;j<type.length;j++){
				for(int i=0;i<9;i++){
					typeout[k] = "#"+String.valueOf(l+1)+type[j]+subs[i];
					k++;
				}
			}


			String[] NO = {"PB2","PB1","PA","H1","NP","N1","M1","M2","NS1","NS2"};
			DataSet d = new DataSet("/Users/Sam/Desktop/FluData/sw_ref_"+NO[x]+".fasta", "/Users/Sam/Desktop/FluData/Table.csv", NO[x]);
			ArrayList<SequenceInfo> output = new ArrayList<SequenceInfo>();
			Iterator<SequenceInfo> It =  d.DataMainFrame.iterator();
			while(It.hasNext()){
				SequenceInfo elementInfo = It.next();
				for(int i=0;i<typeout.length;i++)
					if(elementInfo.Genotype.equals(typeout[i])){
						output.add(elementInfo);
					}
			}
			d.exportFASTA("/Users/sam/Desktop/out2/"+NO[x]+"_"+type[0]+".fa", output);
		}

		/* continent count ***************************************************************************************************/
	/*	int[] array = {0,1,2,3,4,5,6,6,7,7};
		double[][] rel = new double[array.length][4];
		for(int x=0;x<array.length;x++){
			int l=array[x];
			String[] typeout = new String[18];
			//	String[] typeout = new String[9];
			//
			String[] subs = {"","a","b","c","d","e","f","g","h",};
			String[] type={"B1"};
			int k=0;
			for(int j=0;j<type.length;j++){
				for(int i=0;i<9;i++){
					typeout[k] = "#"+String.valueOf(l+1)+type[j]+subs[i];
					k++;
				}
			}
			double[] continent = new double[4];

			String[] NO = {"PB2","PB1","PA","H1","NP","N1","M1","M2","NS1","NS2"};
			teaspoon.adaptation.DataSet d = new teaspoon.adaptation.DataSet("/Users/sam/Desktop/FluData/sw_ref_"+NO[x]+".fasta", "/Users/sam/Desktop/FluData/Table.csv", NO[x]);
			Iterator<teaspoon.adaptation.SequenceInfo> It =  d.DataMainFrame.iterator();
			while(It.hasNext()){
				teaspoon.adaptation.SequenceInfo elementInfo = It.next();
				for(int i=0;i<typeout.length;i++){
					if(elementInfo.Genotype.equals(typeout[i])){
						if(elementInfo.Continent.equals("AS")){
							continent[0]++;
						}
						else if(elementInfo.Continent.equals("EU")){
							continent[1]++;
						}
						else if(elementInfo.Continent.equals("NA")){
							continent[2]++;
						}
						else {
							continent[3]++;
						}
					}
				}
			}

			rel[x][0]=continent[0];
			rel[x][1]=continent[1];
			rel[x][2]=continent[2];
			rel[x][3]=continent[3];

		}


		System.out.println("Asia:  "+rel[3][0]);
		System.out.println("Europe:  "+rel[3][1]);
		System.out.println("North America:  "+rel[3][2]);
		System.out.println("Other:  "+rel[3][3]);
*/

		/************************************ HUMAN SEGMENT*/		
		/*						int N=1000;
		teaspoon.adaptation.FluAnalysis ff = new teaspoon.adaptation.FluAnalysis();

		teaspoon.adaptation.Value[][] v = ff.clumpanalysish3n2BS(N);
	//	teaspoon.adaptation.Value[][] v = ff.clumpanalysish1n1();
	//	teaspoon.adaptation.Value[][] v = ff.Swine();
		//	teaspoon.adaptation.Value[][] v = ff.HAanalysisH3N2(N);
		//	System.out.println("---------------------------------");
		//	System.out.println(" 					Results");
		//	System.out.println("---------------------------------");
		for(int i=0;i<v.length;i++){
			System.out.print(v[i][0].row);
			System.out.print("\t");
			for(int j=0;j<v[0].length;j++){

				System.out.print(v[i][j].Nadapt);

				System.out.print("	");
			}
			System.out.println("");
		}


//				File file = new File("/Users/sam/Desktop/boots/HAHU2.txt");
//				File file1 = new File("/Users/sam/Desktop/boots/HA1HU2.txt");
//				File file2 = new File("/Users/sam/Desktop/boots/HA2HU2.txt");
//				File file3 = new File("/Users/sam/Desktop/boots/NAHU2.txt");
//				File file4 = new File("/Users/sam/Desktop/boots/NAHU2s.txt");
//				File file5 = new File("/Users/sam/Desktop/boots/NAHU2i.txt");


				File file = new File("/Users/sam/Desktop/boots/H3HU.txt");
				File file1 = new File("/Users/sam/Desktop/boots/H3HA1HU.txt");
				File file2 = new File("/Users/sam/Desktop/boots/H3HA2HU.txt");
				File file3 = new File("/Users/sam/Desktop/boots/N2HU.txt");
				File file4 = new File("/Users/sam/Desktop/boots/N2HUs.txt");
				File file5 = new File("/Users/sam/Desktop/boots/N2HUi.txt");
		// Send formatted output to the file.
	      try{
	    	    // Create file 
	    	    FileWriter fstream = new FileWriter(file);
	    	    FileWriter fstream1 = new FileWriter(file1);
	    	    FileWriter fstream2 = new FileWriter(file2);
	    	    FileWriter fstream3 = new FileWriter(file3);
	    	    FileWriter fstream4 = new FileWriter(file4);
	    	    FileWriter fstream5 = new FileWriter(file5);
	    	        BufferedWriter out = new BufferedWriter(fstream);
	    	        BufferedWriter out1 = new BufferedWriter(fstream1);
	    	        BufferedWriter out2 = new BufferedWriter(fstream2);
	    	        BufferedWriter out3 = new BufferedWriter(fstream3);
	    	        BufferedWriter out4 = new BufferedWriter(fstream4);
	    	        BufferedWriter out5 = new BufferedWriter(fstream5);
	    			for(int x=0;x<N;x++){
	    				for(int i=0;i<v.length;i++){
	    						out.write(new Double(v[i][3].Bstrap[x]).toString());
	    						out1.write(new Double(v[i][4].Bstrap[x]).toString());
	    						out2.write(new Double(v[i][5].Bstrap[x]).toString());
	    						out3.write(new Double(v[i][7].Bstrap[x]).toString());
	    						out4.write(new Double(v[i][8].Bstrap[x]).toString());
	    						out5.write(new Double(v[i][9].Bstrap[x]).toString());
	    						out.write("\t");
	    						out1.write("\t");
	    						out2.write("\t");
	    						out3.write("\t");
	    						out4.write("\t");
	    						out5.write("\t");


	    				}
	    				out.write("\n");
	    				out1.write("\n");
	    				out2.write("\n");
	    				out3.write("\n");
	    				out4.write("\n");
	    				out5.write("\n");
	    			}

	    	    //Close the output stream
	    	    out.close();
	    	    out1.close();
	    	    out2.close();
	    	    out3.close();
	    	    out4.close();
	    	    out5.close();
	    	    }catch (Exception e){//Catch exception if any
	    	      System.err.println("Error: " + e.getMessage());
	    	    }*/


		//	double[] xArray = {1977.50000000000,1979.65384615385,1984.63333333333,1987.50000000000,1991.30000000000,1993.50000000000,1994.50000000000,1995.50000000000,1996.50000000000,1997.50000000000,1998.50000000000,1999.50000000000,2000.50000000000,2001.50000000000,2002.50000000000,2003.50000000000,2004.50000000000,2005.50000000000,2006.50000000000,2007.50000000000,2008.50000000000,2009.50000000000};
		// observed y data array
		//	double[] yArray = {0,1,8.32673268200000,9.32673268200000,12.7673267465000,16.4752475452500,16.6608911089999,16.4214993135178,21.8906586484456,22.0922259097037,23.1668350736153,20.7771018855275,31.6761698336096,28.9727723010000,30.6664412075667,27.7804505041041,38.7893539771098,37.5273988058092,39.1435643870000,41.5863968071050,41.8075200689523,45.1683168645603};
		//	double[] yArray = {0,1,8.32673268200000,9.32673268200000,12.7673267465000,16.4752475452500,16.6608911089999,16.4214993135178,21.8906586484456,22.0922259097037,23.1668350736153,20.7771018855275,31.6761698336096,28.9727723010000,30.6664412075667,27.7804505041041,39.7893539771098,42.5273988058092,46.1435643870000,49.5863968071050,55.8075200689523,66.1683168645603};

		/*	double[] yArray2 ={0,
				0,
				0,
				3.1921,
				15.9421,
				17.0968,
				18.8438,
				22.0585,
				21.3725};

		double[] xArray2 ={1992.6202,		 
				1995.1038,		 
				1998.3791,		 
				2001.2723,		 
				2003.8679,		 
				2005.2476,		 
				2006.8279,		 
				2008.3257,		 
				2009.4689};	


		double[] yArray ={ 0,
				 14.7229,
				 8.5864,
				 15.9665,
				 10.5123,
				 13.4615,
				 13.4635,
				 11.1553,
				 14.7462,
				 22.9236,
				 15.5902,};


		double[] xArray ={1987.5,	
				1992.6958,
				1995.1038,
				1998.3791,
				2000.8308,
				2003.1072,
				2004.5236,
				2005.3081,
				2006.8345,
				2008.3089,
				2009.4282,
};	

		xArray=xArray;
		yArray=yArray;
		// estimates of the standard deviations of y
		double[] x = new double[xArray.length];
		double init = xArray[0];
		for(int i=0;i<xArray.length;i++){
			x[i] = xArray[i]-init;
		}


		Regression reg = new Regression(x, yArray);
		reg.supressYYplot();
		reg.bestPolynomialPlot(0.0);*/
		//	reg.fiveParameterLogisticPlot(0.0, 22.0);


		//	reg.polynomialPlot(4,0.0);


		//	reg2.polynomial(1,0);
		//	reg.polynomial(1,1975.5);
		//		reg.print();
		//	double[] est = reg.getBestEstimates();
		//	System.out.println(est[0]);
		//	double[] est2 = reg2.getBestEstimates();
		//	System.out.println(est2[0]);

		/*	System.out.println("*********Lets start the program This is a test run*********\n\n");
		GA G=new GA();
		teaspoon.adaptation.Hemagglutinin HD = new teaspoon.adaptation.Hemagglutinin();

		int[] OG = HD.HA;
		for(int i=0;i<OG.length;i++){
			if(OG[i]==2){
				OG[i]=1;
			}
		}



		int[] par = G.Maximise();
		for(int i=0;i<par.length;i++){
			System.out.print(par[i]);
		}*/
		//		System.out.println("\n\n");
		//		for(int i=0;i<G.partition.length;i++){
		//			System.out.print(OG[i]);
		//		}

		/*		FluAnalysis2 f = new FluAnalysis2();
		String[] EA = {"PB2_B1_[10.11.11].txt",
				"PB1_B1_[10.11.11].txt",
				"PA_B1_[10.11.11].txt",
				"H1_B1_HA_[10.11.11].txt",
				"H1_B1_HA1_[10.11.11].txt",
				"H1_B1_HA2_[10.11.11].txt",
				"H1_B1_HAsi_[10.11.11].txt",
				"H1_B1_HAs_[10.11.11].txt",
				"H1_B1_HAi_[10.11.11].txt",
				"NP_B1_[10.11.11].txt",
				"N1_B1_NA_[10.11.11].txt",
				"N1_B1_NAs_[10.11.11].txt",
				"N1_B1_NAi_[10.11.11].txt",
				"M1_B1_[10.11.11].txt",
				"M2_B1_[10.11.11].txt",
				"NS1_B1_[10.11.11].txt",
				"NS2_B1_[10.11.11].txt"};

		String[] CS = {"PB2_C1_[10.11.11].txt",
				"PB1_C1_[10.11.11].txt",
				"PA_C1_[10.11.11].txt",
				"H1_C1_HA_[10.11.11].txt",
				"H1_C1_HA1_[10.11.11].txt",
				"H1_C1_HA2_[10.11.11].txt",
				"H1_C1_HAsi_[10.11.11].txt",
				"H1_C1_HAs_[10.11.11].txt",
				"H1_C1_HAi_[10.11.11].txt",
				"NP_C1_[10.11.11].txt",
				"N1_C1_NA_[10.11.11].txt",
				"N1_C1_NAs_[10.11.11].txt",
				"N1_C1_NAi_[10.11.11].txt",
				"M1_C1_[10.11.11].txt",
				"M2_C1_[10.11.11].txt",
				"NS1_C1_[10.11.11].txt",
				"NS2_C1_[10.11.11].txt"};
		String[] TRIG = {"PB2_T1_[10.11.11].txt",
				"PB1_T1_[10.11.11].txt",
				"PA_T1_[10.11.11].txt",
				"H1_T1_HA_[10.11.11].txt",
				"H1_T1_HA1_[10.11.11].txt",
				"H1_T1_HA2_[10.11.11].txt",
				"H1_T1_HAsi_[10.11.11].txt",
				"H1_T1_HAs_[10.11.11].txt",
				"H1_T1_HAi_[10.11.11].txt",
				"NP_T1_[10.11.11].txt",
				"N1_T1_NA_[10.11.11].txt",
				"N1_T1_NAs_[10.11.11].txt",
				"N1_T1_NAi_[10.11.11].txt",
				"M1_T1_[10.11.11].txt",
				"M2_T1_[10.11.11].txt",
				"NS1_T1_[10.11.11].txt",
				"NS2_T1_[10.11.11].txt"};
		String[] H1N1 = {"PB2_H1N1_[04.11.11].txt",
				"PB1_H1N1_[04.11.11].txt",
				"PA_H1N1_[04.11.11].txt",
				"HA_H1N1_[04.11.11].txt",
				"HA1_H1N1_[04.11.11].txt",
				"HA2_H1N1_[04.11.11].txt",
				"HAsi_H1N1_[04.11.11].txt",
				"HAs_H1N1_[04.11.11].txt",
				"HAi_H1N1_[04.11.11].txt",
				"NP_H1N1_[04.11.11].txt",
				"NA_H1N1_[04.11.11].txt",
				"NAs_H1N1_[04.11.11].txt",
				"NAi_H1N1_[04.11.11].txt",
				"M1_H1N1_[04.11.11].txt",
				"M2_H1N1_[04.11.11].txt",
				"NS1_H1N1_[04.11.11].txt",
				"NS2_H1N1_[04.11.11].txt"};
		String[] H3N2 = {"PB2_H3N2_[04.11.11].txt",
				"PB1_H3N2_[04.11.11].txt",
				"PA_H3N2_[04.11.11].txt",
				"HA_H3N2_[04.11.11].txt",
				"HA1_H3N2_[04.11.11].txt",
				"HA2_H3N2_[04.11.11].txt",
				"HAsi_H3N2_[04.11.11].txt",
				"HAs_H3N2_[04.11.11].txt",
				"HAi_H3N2_[04.11.11].txt",
				"NP_H3N2_[04.11.11].txt",
				"NA_H3N2_[04.11.11].txt",
				"NAs_H3N2_[04.11.11].txt",
				"NAi_H3N2_[04.11.11].txt",
				"M1_H3N2_[04.11.11].txt",
				"M2_H3N2_[04.11.11].txt",
				"NS1_H3N2_[04.11.11].txt",
				"NS2_H3N2_[04.11.11].txt"};
		int[] l = {3,6,10,11,14,15,22};
		int[] word = {2,12,7,8,6,7,8};
		String loc = "/Users/sam/Desktop/out/Final8772/";
		for(int i=0;i<l.length;i++){
			f.ExtractVal(CS, l[i],word[i],loc);
			f.ExtractVal(EA, l[i],word[i],loc);
			f.ExtractVal(H1N1, l[i],word[i],loc);
			f.ExtractVal(H3N2,  l[i],word[i],loc);
			f.ExtractVal(TRIG,  l[i],word[i],loc);
			System.out.println("\n");
		}

		 */

	}
}












