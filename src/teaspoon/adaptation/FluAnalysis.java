package teaspoon.adaptation;

import java.io.File;
import java.util.Random;

public class FluAnalysis {

	public String folder;
	public String type;
	public String gene;
	public String location;
	public String date;
	public boolean isans=false;
	public String[] dates;
	public String[] types;
	public String[] genes;
	public double[] dateA;
	public int[] list;
	public int[] list2;
	public int[] list3;
	Random generator = new Random();
	final double[] low = {0.0,0.15};
	final double[] mid = {0.15,0.75};
	final double[] high = {0.75,1.0};	
	public FluAnalysis(){
		//	System.out.println("INTERVALS");
		//	System.out.println("low bound:  "+low[0]+"--->"+low[1]);
		//	System.out.println("mid bound:  "+mid[0]+"--->"+mid[1]);
		//	System.out.println("high bound:  "+high[0]+"--->"+high[1]);
	}



	public Value[][] clumpanalysish1n1(){		
		//		String[] datesh1n1 ={"1977","1978","1979","1980","1981","1982","1983","1984","1985","1986","1987","1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009"};
		//		String[] datesh1n1 = {"1977_1979","1983_1985","1995_1997","1998_2000","2001_2003","2004_2006","2007_2009"};
		String[] datesh1n1 = {"1977","1978-1979-1980","1982-1983-1984","1995","1996","1999-2000","2001","2002-2003","2004-2005","2006","2007","2008","2009"};
		//		String[] datesh1n1 = {"1977_1979","1982-1983-1984","1995","1996","1999-2000","2001","2002-2003","2004-2005","2006","2007","2008","2009"};

		String[] typesh1n1 = {"H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1"};
		String[] genesh1n1 = {"PB2","PB1","PA","HA","HA1","HA2","NP","NA","NA","NA","M1","M2","NS1"};

		//	String[] typesh1n1 = {"H1N1"};
		//	String[] genesh1n1 = {"HA"};
		Neuraminidase Nd = new Neuraminidase();
		double[] nr00={0.347958181,0.262833676,0.261582734,0.271944923,0.271944923,0.271944923,0.236294093,0.358890701,0.358890701,0.358890701,0.160869565,1.115241636,0.809906292};
		double[] nr05={0.189444501,0.15912031,0.136080332,0.125984252,0.125984252,0.125984252,0.120689655,0.152284264,0.152284264,0.152284264,0.07960199,1.130434783,0.573099415};
		double[] nr10={0.18256366,0.130558184,0.11713367,0.166089965,0.166089965,0.166089965,0.112903226,0.15862069,0.15862069,0.15862069,0.123076923,1.28,0.640776699};
		double[] nr15={0.175203447,0.09784324,0.10302409,0.19047619,0.19047619,0.19047619,0.125348189,0.153225806,0.153225806,0.153225806,0.105263158,1.16,0.846153846};
		double[] nr20={0.183238636,0.114858706,0.113293051,0.171122995,0.171122995,0.171122995,0.127659574,0.142857143,0.142857143,0.142857143,0.129032258,1.25,0.685185185};
		double[] nr25={0.196911197,0.126237624,0.107142857,0.154929577,0.154929577,0.154929577,0.152777778,0.15942029,0.15942029,0.15942029,0,1.25,0.608695652};
		double[] nr30={0.24,0.117437722,0.090778098,0.189473684,0.189473684,0.189473684,0.071428571,0.153846154,0.153846154,0.153846154,0,3,0.631578947};
		double[] nr35={0.266666667,0.117437722,0.090778098,0.105263158,0.105263158,0.105263158,0.076271186,0.10989011,0.10989011,0.10989011,0,3,0.631578947};
		double[] nr40={0.276497696,0.119771863,0.101286174,0.063157895,0.063157895,0.063157895,0.076271186,0.10989011,0.10989011,0.10989011,0,0,0.526315789};
		double[] nr45={0.262910798,0.125498008,0.10862069,0.063157895,0.063157895,0.063157895,0.076271186,0.120481928,0.120481928,0.120481928,0,0,0.526315789};
		//	double[] neutral_ratio = nr15;
		double[] nr = {0.207132915,0.133777952,0.123544093,0.211764706,0.211764706,0.211764706,0.17248062,0.151624549,0.151624549,0.151624549,0.179487179,1.469135802,0.594147093};// swine mainAnalysis
		//	double[] neutral_ratio = {0.211};
		String ansname = new String("");
		Value[][] value_matrix = new Value[datesh1n1.length][typesh1n1.length];
		for(int i=0;i<value_matrix.length;i++){
			for(int j=0;j<value_matrix[0].length;j++){
				value_matrix[i][j] = new Value(1000);
			}
		}	
		int f =0;
		for(int t=0,g=0,l=0;t<genesh1n1.length;t++, g++,l++){
			int count=0;
			for(int d=0;d<datesh1n1.length;d++){
				String whichfile = new StringBuffer().append(datesh1n1[d]).append(".").append(genesh1n1[g]).append(".").append(typesh1n1[t]).append(".").append("nex").toString();

				//	 	File fle = new File("/Users/sam/Desktop/FluData/Hu/"+ whichfile);
				File fle = new File("/Users/sam/Desktop/FluData/SplitSingle/"+ whichfile);
				// file reading ***********************************************************
				if(fle.exists()){
					// use first time point as ans
					count++;
					if(count==1){
						ansname = fle.toString();
					} 				

					Read_main ra = new Read_main(ansname);
					int[][] tmp = ra.readNEXUS();
					int[] ans = ra.consensusArray(tmp);
					//		int sequenceNum=num;  //max12
					//		for(int i=0;i<tmp[0].length;i++){
					//			ans[i] = tmp[sequenceNum][i];
					//		}

					Read_main r = new Read_main(fle.toString());
					int[][] seq = r.readNEXUS();	
					//na block *************
					if(t==8){
						Methods m = new Methods();
						int[][] seqNA = m.Subsetter(seq, Nd.badlisth1n1, 1);
						int[] ansNA = m.Subsetter(ans, Nd.badlisth1n1, 1);
						Williamson3bin w = new Williamson3bin(seqNA,ansNA);
						w.williamson3bin_method(nr[t],low,mid,high);
						value_matrix[d][t].Nadapt = w.Adapt;
						//		value_matrix[d][t].Nadapt = w.mid_S;
						value_matrix[d][t].row = datesh1n1[d];
						value_matrix[d][t].column = genesh1n1[t];
						value_matrix[d][t].dataset = "H1N1";
					} else if (t==9){
						Methods m = new Methods();
						int[][] seqNA = m.Subsetter(seq, Nd.badlisth1n1, 0);
						int[] ansNA = m.Subsetter(ans, Nd.badlisth1n1, 0);
						Williamson3bin w = new Williamson3bin(seqNA,ansNA);
						w.williamson3bin_method(nr[t],low,mid,high);
						value_matrix[d][t].Nadapt =  w.Adapt;
						//		value_matrix[d][t].Nadapt = w.mid_S;
						value_matrix[d][t].row = datesh1n1[d];
						value_matrix[d][t].column = genesh1n1[t];
						value_matrix[d][t].dataset = "H1N1";
					} else {
						Williamson3bin w = new Williamson3bin(seq,ans);
						w.williamson3bin_method(nr[t],low,mid,high);
						value_matrix[d][t].Nadapt =  w.Adapt;
						//		value_matrix[d][t].Nadapt = w.mid_S;
						value_matrix[d][t].row = datesh1n1[d];
						value_matrix[d][t].column = genesh1n1[t];
						value_matrix[d][t].dataset = "H1N1";
					}
				}

			}
		}
		return value_matrix;
	}

	public WilliamsonOutputClass[] HumanRunner(int N,double[] nr,String which){		
		if(which.equals("H1N1")){
			String[] datesh1n1 = {"1977","1978-1979-1980","1982-1983-1984","1995","1996","1999-2000","2001","2002-2003","2004-2005","2006","2007","2008","2009"};		
			double[] date =  {1977.5,1978.9375,1983.460784,1995.5,1996.5,2000.428571,2001.5,2003.269231,2005.447368,2006.5,2007.5,2008.5,2009.5};
			String[] typesh1n1 = {"H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1"};
			String[] genesh1n1 = {"PB2","PB1","PA","HA","HA","HA","HA","HA","HA","NP","NA","NA","NA","M1","M2","NS1"};
			dates=datesh1n1;
			types=typesh1n1;
			genes=genesh1n1;
			Neuraminidase Nd = new Neuraminidase();
			Hemagglutinin Hd = new Hemagglutinin();
			list3 = Hd.HAH1;
			list2 = Hd.H1;
			list=Nd.badlisth1n1;
			dateA=date;
		} else {
			String[] datesh3n2 ={"1977","1978-1979-1980","1983-1984-1985","1986-1987-1988","1990-1991-1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009"};  // new dates
			double[] date =  {1977.5,1979.653846,1984.633333,1987.5,1991.3,1993.5,1994.5,1995.5,1996.5,1997.5,1998.5,1999.5,2000.5,2001.5,2002.5,2003.5,2004.5,2005.5,2006.5,2007.5,2008.5,2009.5};	
			String[] typesh3n2 = {"H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2"};
			String[] genesh3n2 = {"PB2","PB1","PA","HA","HA","HA","HA","HA","HA","NP","NA","NA","NA","M1","M2","NS1"};
			dates=datesh3n2;
			types=typesh3n2;
			genes=genesh3n2;
			Hemagglutinin Hd = new Hemagglutinin();
			list3 = Hd.HAH3;
			list2 = Hd.H3;
			Neuraminidase Nd = new Neuraminidase();
			list=Nd.badlisth3n2;
			dateA=date;
		}



		// initialise williamson output classes
		String ansname = new String("");
		WilliamsonOutputClass[] woc = new WilliamsonOutputClass[genes.length];
		for(int u=0;u<woc.length;u++){
			woc[u]=new WilliamsonOutputClass(dates.length,N);
		}

		for(int run=0;run<N;run++){   // run through all boot straps
			for(int t=0,g=0,l=0;t<genes.length;t++, g++,l++){ // cycle through all genes
				int count=0;
				String whichfileans = new StringBuffer().append(dates[0]).append(".").append(genes[g]).append(".").append(types[t]).append(".").append("nex").toString();
				File fleans = new File("/Users/sam/Desktop/FluData/SplitSingle/"+ whichfileans);
				ansname = fleans.toString();
				Read_main ra = new Read_main(ansname); // read ans
				int[][] tmp = ra.readNEXUS();
				int[] ans = ra.consensusArray(tmp);			
				int[] sampler = new int[ans.length/3];
				Random generator = new Random();
				for(int x=0;x<sampler.length;x++){
					int randint = generator.nextInt(sampler.length-1);
					sampler[x] = randint;
				}
				if(t==4){ // subset for HA list surface
					int numsites=0;
					for(int j=0;j<list2.length;j++){if(list2[j]==1){numsites++;}}
					sampler = new int[(numsites)];
					Random generator2 = new Random();
					for(int x=0;x<sampler.length;x++){
						int randint = generator2.nextInt(sampler.length-1);
						sampler[x] = randint;
					}
				}
				if(t==5){ // subset for HA list internal
					int numsites=0;
					for(int j=0;j<list2.length;j++){if(list2[j]==0){numsites++;}}
					sampler = new int[(numsites)];
					Random generator2 = new Random();
					for(int x=0;x<sampler.length;x++){
						int randint = generator2.nextInt(sampler.length-1);
						sampler[x] = randint;
					}
				}

				if(t==7){ // subset for H1
					int numsites=0;
					for(int j=0;j<list3.length;j++){if(list3[j]==1){numsites++;}}
					sampler = new int[(numsites)];
					Random generator2 = new Random();
					for(int x=0;x<sampler.length;x++){
						int randint = generator2.nextInt(sampler.length-1);
						sampler[x] = randint;
					}
				}
				if(t==8){ // subset for H2
					int numsites=0;
					for(int j=0;j<list3.length;j++){if(list3[j]==0){numsites++;}}
					sampler = new int[(numsites)];
					Random generator2 = new Random();
					for(int x=0;x<sampler.length;x++){
						int randint = generator2.nextInt(sampler.length-1);
						sampler[x] = randint;
					}
				}
				

				if(t==11){ // subset for NA list surface
					int numsites=0;
					for(int j=0;j<list.length;j++){if(list[j]==1){numsites++;}}
					sampler = new int[(numsites)];
					Random generator2 = new Random();
					for(int x=0;x<sampler.length;x++){
						int randint = generator2.nextInt(sampler.length-1);
						sampler[x] = randint;
					}
				}
				if(t==12){ // subset for NA list internal   //this doesn't sample the internal codons randomly
					int numsites=0;
					for(int j=0;j<list.length;j++){
                        if(list[j]==0){
                            numsites++;
                        }
                    }
					sampler = new int[(numsites)];
					Random generator2 = new Random();
					for(int x=0;x<sampler.length;x++){
						int randint = generator2.nextInt(sampler.length-1);
						sampler[x] = randint;
					}
				}

				for(int d=0;d<dates.length;d++){	// cycle through all time points			

					String whichfile = new StringBuffer().append(dates[d]).append(".").append(genes[g]).append(".").append(types[t]).append(".").append("nex").toString();
					File fle = new File("/Users/sam/Desktop/FluData/SplitSingle/"+ whichfile);
					// file reading ***********************************************************
					if(fle.exists()){
						// use first time point as ans
						count++;
						Read_main r = new Read_main(fle.toString());
						int[][] seq = r.readNEXUS();	
						// ha block
//********************************************************************************************************************************************************							
						if(t==4){		// ha surface
							Methods m = new Methods();
							int[][] seqHA = m.Subsetter(seq, list2, 1);
							int[] ansHA = m.Subsetter(ans, list2, 1);
							Williamson3bin ww = new Williamson3bin(seqHA,ansHA); //******
							Store s = ww.CreateBlocks(3,seqHA[0].length,sampler); //******
							Williamson3bin ws = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							ws.williamson3bin_method(nr[t],low,mid,high);

							if(run==0 && d==0){
								woc[g].Low_S.add(Double.valueOf(0.0));
								woc[g].Low_R.add(Double.valueOf(0.0));
								woc[g].Mid_S.add(Double.valueOf(0.0));
								woc[g].Mid_R.add(Double.valueOf(0.0));
								woc[g].High_S.add(Double.valueOf(0.0));
								woc[g].High_R.add(Double.valueOf(0.0));
								woc[g].Fix_S.add(Double.valueOf(0.0));
								woc[g].Fix_R.add(Double.valueOf(0.0));
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(Double.valueOf(0.0));
								woc[g].MidA.add(Double.valueOf(0.0));
								woc[g].HighA.add(Double.valueOf(0.0));
								woc[g].FixA.add(Double.valueOf(0.0));
								woc[g].TotalA.add(Double.valueOf(0.0));
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansHA.length;
							}
							else if(run==0){
								ws = new Williamson3bin(seqHA,ansHA); //******
								ws.williamson3bin_method(nr[t],low,mid,high);							
								woc[g].Low_S.add(ws.low_S);
								woc[g].Low_R.add(ws.low_R);
								woc[g].Mid_S.add(ws.mid_S);
								woc[g].Mid_R.add(ws.mid_R);
								woc[g].High_S.add(ws.high_S);
								woc[g].High_R.add(ws.high_R);
								woc[g].Fix_S.add(ws.fix_S);
								woc[g].Fix_R.add(ws.fix_R);
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(ws.low_A);
								woc[g].MidA.add(ws.mid_A);
								woc[g].HighA.add(ws.high_A);
								woc[g].FixA.add(ws.fix_A);
								woc[g].TotalA.add(ws.Adapt);
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansHA.length;
							}
							woc[g].BS.get(d)[run]=ws.Adapt;
						} 
						else if(t==5){	// ha internal
							Methods m = new Methods();
							int[][] seqHA = m.Subsetter(seq, list2, 0);
							int[] ansHA = m.Subsetter(ans, list2, 0);
							Williamson3bin ww = new Williamson3bin(seqHA,ansHA); //******
							Store s = ww.CreateBlocks(3,seqHA[0].length,sampler); //******
							Williamson3bin ws = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							ws.williamson3bin_method(nr[t],low,mid,high);

							if(run==0 && d==0){
								woc[g].Low_S.add(Double.valueOf(0.0));
								woc[g].Low_R.add(Double.valueOf(0.0));
								woc[g].Mid_S.add(Double.valueOf(0.0));
								woc[g].Mid_R.add(Double.valueOf(0.0));
								woc[g].High_S.add(Double.valueOf(0.0));
								woc[g].High_R.add(Double.valueOf(0.0));
								woc[g].Fix_S.add(Double.valueOf(0.0));
								woc[g].Fix_R.add(Double.valueOf(0.0));
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(Double.valueOf(0.0));
								woc[g].MidA.add(Double.valueOf(0.0));
								woc[g].HighA.add(Double.valueOf(0.0));
								woc[g].FixA.add(Double.valueOf(0.0));
								woc[g].TotalA.add(Double.valueOf(0.0));
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansHA.length;
							}
							else if(run==0){
								ws = new Williamson3bin(seqHA,ansHA); //******
								ws.williamson3bin_method(nr[t],low,mid,high);							
								woc[g].Low_S.add(ws.low_S);
								woc[g].Low_R.add(ws.low_R);
								woc[g].Mid_S.add(ws.mid_S);
								woc[g].Mid_R.add(ws.mid_R);
								woc[g].High_S.add(ws.high_S);
								woc[g].High_R.add(ws.high_R);
								woc[g].Fix_S.add(ws.fix_S);
								woc[g].Fix_R.add(ws.fix_R);
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(ws.low_A);
								woc[g].MidA.add(ws.mid_A);
								woc[g].HighA.add(ws.high_A);
								woc[g].FixA.add(ws.fix_A);
								woc[g].TotalA.add(ws.Adapt);
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansHA.length;
							}
							woc[g].BS.get(d)[run]=ws.Adapt;
						} 
//********************************************************************************************************************************************************							
						else if(t==7){		// ha surface
							Methods m = new Methods();
							int[][] seqHA = m.Subsetter(seq, list3, 1);
							int[] ansHA = m.Subsetter(ans, list3, 1);
							Williamson3bin ww = new Williamson3bin(seqHA,ansHA); //******
							Store s = ww.CreateBlocks(3,seqHA[0].length,sampler); //******
							Williamson3bin ws = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							ws.williamson3bin_method(nr[t],low,mid,high);

							if(run==0 && d==0){
								woc[g].Low_S.add(Double.valueOf(0.0));
								woc[g].Low_R.add(Double.valueOf(0.0));
								woc[g].Mid_S.add(Double.valueOf(0.0));
								woc[g].Mid_R.add(Double.valueOf(0.0));
								woc[g].High_S.add(Double.valueOf(0.0));
								woc[g].High_R.add(Double.valueOf(0.0));
								woc[g].Fix_S.add(Double.valueOf(0.0));
								woc[g].Fix_R.add(Double.valueOf(0.0));
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(Double.valueOf(0.0));
								woc[g].MidA.add(Double.valueOf(0.0));
								woc[g].HighA.add(Double.valueOf(0.0));
								woc[g].FixA.add(Double.valueOf(0.0));
								woc[g].TotalA.add(Double.valueOf(0.0));
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansHA.length;
							}
							else if(run==0){
								ws = new Williamson3bin(seqHA,ansHA); //******
								ws.williamson3bin_method(nr[t],low,mid,high);							
								woc[g].Low_S.add(ws.low_S);
								woc[g].Low_R.add(ws.low_R);
								woc[g].Mid_S.add(ws.mid_S);
								woc[g].Mid_R.add(ws.mid_R);
								woc[g].High_S.add(ws.high_S);
								woc[g].High_R.add(ws.high_R);
								woc[g].Fix_S.add(ws.fix_S);
								woc[g].Fix_R.add(ws.fix_R);
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(ws.low_A);
								woc[g].MidA.add(ws.mid_A);
								woc[g].HighA.add(ws.high_A);
								woc[g].FixA.add(ws.fix_A);
								woc[g].TotalA.add(ws.Adapt);
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansHA.length;
							}
							woc[g].BS.get(d)[run]=ws.Adapt;
						} 
						else if(t==8){	// ha internal
							Methods m = new Methods();
							int[][] seqHA = m.Subsetter(seq, list3, 0);
							int[] ansHA = m.Subsetter(ans, list3, 0);
							Williamson3bin ww = new Williamson3bin(seqHA,ansHA); //******
							Store s = ww.CreateBlocks(3,seqHA[0].length,sampler); //******
							Williamson3bin ws = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							ws.williamson3bin_method(nr[t],low,mid,high);

							if(run==0 && d==0){
								woc[g].Low_S.add(Double.valueOf(0.0));
								woc[g].Low_R.add(Double.valueOf(0.0));
								woc[g].Mid_S.add(Double.valueOf(0.0));
								woc[g].Mid_R.add(Double.valueOf(0.0));
								woc[g].High_S.add(Double.valueOf(0.0));
								woc[g].High_R.add(Double.valueOf(0.0));
								woc[g].Fix_S.add(Double.valueOf(0.0));
								woc[g].Fix_R.add(Double.valueOf(0.0));
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(Double.valueOf(0.0));
								woc[g].MidA.add(Double.valueOf(0.0));
								woc[g].HighA.add(Double.valueOf(0.0));
								woc[g].FixA.add(Double.valueOf(0.0));
								woc[g].TotalA.add(Double.valueOf(0.0));
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansHA.length;
							}
							else if(run==0){
								ws = new Williamson3bin(seqHA,ansHA); //******
								ws.williamson3bin_method(nr[t],low,mid,high);							
								woc[g].Low_S.add(ws.low_S);
								woc[g].Low_R.add(ws.low_R);
								woc[g].Mid_S.add(ws.mid_S);
								woc[g].Mid_R.add(ws.mid_R);
								woc[g].High_S.add(ws.high_S);
								woc[g].High_R.add(ws.high_R);
								woc[g].Fix_S.add(ws.fix_S);
								woc[g].Fix_R.add(ws.fix_R);
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(ws.low_A);
								woc[g].MidA.add(ws.mid_A);
								woc[g].HighA.add(ws.high_A);
								woc[g].FixA.add(ws.fix_A);
								woc[g].TotalA.add(ws.Adapt);
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansHA.length;
							}
							woc[g].BS.get(d)[run]=ws.Adapt;
						} 
						
						
//********************************************************************************************************************************************************							
						//na block *************
						else if(t==11){	// na surface
							Methods m = new Methods();
							int[][] seqNA = m.Subsetter(seq, list, 1);
							int[] ansNA = m.Subsetter(ans, list, 1);
							Williamson3bin ww = new Williamson3bin(seqNA,ansNA); //******
							Store s = ww.CreateBlocks(3,seqNA[0].length,sampler); //******
							Williamson3bin ws = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							ws.williamson3bin_method(nr[t],low,mid,high);

							if(run==0 && d==0){
								woc[g].Low_S.add(Double.valueOf(0.0));
								woc[g].Low_R.add(Double.valueOf(0.0));
								woc[g].Mid_S.add(Double.valueOf(0.0));
								woc[g].Mid_R.add(Double.valueOf(0.0));
								woc[g].High_S.add(Double.valueOf(0.0));
								woc[g].High_R.add(Double.valueOf(0.0));
								woc[g].Fix_S.add(Double.valueOf(0.0));
								woc[g].Fix_R.add(Double.valueOf(0.0));
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(Double.valueOf(0.0));
								woc[g].MidA.add(Double.valueOf(0.0));
								woc[g].HighA.add(Double.valueOf(0.0));
								woc[g].FixA.add(Double.valueOf(0.0));
								woc[g].TotalA.add(Double.valueOf(0.0));
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansNA.length;
							}
							else if(run==0){

								ws = new Williamson3bin(seqNA,ansNA); //******
								ws.williamson3bin_method(nr[t],low,mid,high);							
								woc[g].Low_S.add(ws.low_S);
								woc[g].Low_R.add(ws.low_R);
								woc[g].Mid_S.add(ws.mid_S);
								woc[g].Mid_R.add(ws.mid_R);
								woc[g].High_S.add(ws.high_S);
								woc[g].High_R.add(ws.high_R);
								woc[g].Fix_S.add(ws.fix_S);
								woc[g].Fix_R.add(ws.fix_R);
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(ws.low_A);
								woc[g].MidA.add(ws.mid_A);
								woc[g].HighA.add(ws.high_A);
								woc[g].FixA.add(ws.fix_A);
								woc[g].TotalA.add(ws.Adapt);
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansNA.length;
							}
							woc[g].BS.get(d)[run]=ws.Adapt;
						} else if (t==12){
							Methods m = new Methods();
							int[][] seqNA = m.Subsetter(seq, list, 0);
							int[] ansNA = m.Subsetter(ans, list, 0);
							Williamson3bin ww = new Williamson3bin(seqNA,ansNA); //******
							Store s = ww.CreateBlocks(3,seqNA[0].length,sampler); //******
							Williamson3bin ws = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							ws.williamson3bin_method(nr[t],low,mid,high);

							if(run==0 && d==0){
								woc[g].Low_S.add(Double.valueOf(0.0));
								woc[g].Low_R.add(Double.valueOf(0.0));
								woc[g].Mid_S.add(Double.valueOf(0.0));
								woc[g].Mid_R.add(Double.valueOf(0.0));
								woc[g].High_S.add(Double.valueOf(0.0));
								woc[g].High_R.add(Double.valueOf(0.0));
								woc[g].Fix_S.add(Double.valueOf(0.0));
								woc[g].Fix_R.add(Double.valueOf(0.0));
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(Double.valueOf(0.0));
								woc[g].MidA.add(Double.valueOf(0.0));
								woc[g].HighA.add(Double.valueOf(0.0));
								woc[g].FixA.add(Double.valueOf(0.0));
								woc[g].TotalA.add(Double.valueOf(0.0));
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansNA.length;
							}
							else if(run==0){
								ws = new Williamson3bin(seqNA,ansNA); //******
								ws.williamson3bin_method(nr[t],low,mid,high);

								woc[g].Low_S.add(ws.low_S);
								woc[g].Low_R.add(ws.low_R);
								woc[g].Mid_S.add(ws.mid_S);
								woc[g].Mid_R.add(ws.mid_R);
								woc[g].High_S.add(ws.high_S);
								woc[g].High_R.add(ws.high_R);
								woc[g].Fix_S.add(ws.fix_S);
								woc[g].Fix_R.add(ws.fix_R);
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(ws.low_A);
								woc[g].MidA.add(ws.mid_A);
								woc[g].HighA.add(ws.high_A);
								woc[g].FixA.add(ws.fix_A);
								woc[g].TotalA.add(ws.Adapt);
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ansNA.length;
							}
							woc[g].BS.get(d)[run]=ws.Adapt;
//********************************************************************************************************************************************************							
							// all other genes
						} else {
							Williamson3bin ww = new Williamson3bin(seq,ans); //******
							Store s = ww.CreateBlocks(3,seq[0].length,sampler); //******
							Williamson3bin ws = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							ws.williamson3bin_method(nr[t],low,mid,high);
							if(run==0 && d==0){
								woc[g].Low_S.add(Double.valueOf(0.0));
								woc[g].Low_R.add(Double.valueOf(0.0));
								woc[g].Mid_S.add(Double.valueOf(0.0));
								woc[g].Mid_R.add(Double.valueOf(0.0));
								woc[g].High_S.add(Double.valueOf(0.0));
								woc[g].High_R.add(Double.valueOf(0.0));
								woc[g].Fix_S.add(Double.valueOf(0.0));
								woc[g].Fix_R.add(Double.valueOf(0.0));
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(Double.valueOf(0.0));
								woc[g].MidA.add(Double.valueOf(0.0));
								woc[g].HighA.add(Double.valueOf(0.0));
								woc[g].FixA.add(Double.valueOf(0.0));
								woc[g].TotalA.add(Double.valueOf(0.0));
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ans.length;
							}
							else if(run==0){
								ws = new Williamson3bin(seq,ans); //******
								ws.williamson3bin_method(nr[t],low,mid,high);
								woc[g].Low_S.add(ws.low_S);
								woc[g].Low_R.add(ws.low_R);
								woc[g].Mid_S.add(ws.mid_S);
								woc[g].Mid_R.add(ws.mid_R);
								woc[g].High_S.add(ws.high_S);
								woc[g].High_R.add(ws.high_R);
								woc[g].Fix_S.add(ws.fix_S);
								woc[g].Fix_R.add(ws.fix_R);
								woc[g].dateavg.add(dateA[d]);
								woc[g].LowA.add(ws.low_A);
								woc[g].MidA.add(ws.mid_A);
								woc[g].HighA.add(ws.high_A);
								woc[g].FixA.add(ws.fix_A);
								woc[g].TotalA.add(ws.Adapt);
								woc[g].numSamples.add(seq.length);
								woc[g].neut=nr[t];
								woc[g].L=ans.length;
							}
							woc[g].BS.get(d)[run]=ws.Adapt;

						}
					}

				}
			}
			System.out.println("I am on Run  "+run+"  of "+ N);
		}

		return woc;

	}


	public WilliamsonOutputClass clumpanalysish1n1BS(int N){		
		String[] datesh1n1 = {"1977","1978-1979-1980","1982-1983-1984","1995","1996","1999-2000","2001","2002-2003","2004-2005","2006","2007","2008","2009"};
		double[] dateA =  {1977.5,1978.9375,1983.460784,1995.5,1996.5,2000.428571,2001.5,2003.269231,2005.447368,2006.5,2007.5,2008.5,2009.5};

		//		String[] typesh1n1 = {"H1N1"};
		//		String[] genesh1n1 = {"HA"};

		//	String[] datesh1n1 = {"1977_1979","1983_1985","1995_1997","1998_2000","2001_2003","2004_2006","2007_2009"};
		String[] typesh1n1 = {"H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1","H1N1"};
		String[] genesh1n1 = {"PB2","PB1","PA","HA","HA1","HA2","NP","NA","NA","NA","M1","M2","NS1"};
		Neuraminidase Nd = new Neuraminidase();

		String ansname = new String("");
		WilliamsonOutputClass woc = new WilliamsonOutputClass(datesh1n1.length,N);
		//	double[] neutral_ratio = {0.185604529,0.138175376,0.115634474,0.163822526,0.163822526,0.163822526,0.1367713,0.154362416,0.154362416,0.154362416,0.153846154,1.28,0.640776699};
		//		double[] neutral_ratio={0.175203447,0.09784324,0.10302409,0.19047619,0.19047619,0.19047619,0.125348189,0.153225806,0.153225806,0.153225806,0.105263158,1.16,0.846153846};
		double[] nr = {0.207132915,0.133777952,0.123544093,0.211764706,0.211764706,0.211764706,0.17248062,0.151624549,0.151624549,0.151624549,0.179487179,1.469135802,0.594147093};// swine mainAnalysis

		Value[][] value_matrix = new Value[datesh1n1.length][typesh1n1.length];
		for(int i=0;i<value_matrix.length;i++){
			for(int j=0;j<value_matrix[0].length;j++){
				value_matrix[i][j] = new Value(N);
			}
		}	
		for(int run=0;run<N;run++){   // run through all boot straps
			for(int t=0,g=0,l=0;t<genesh1n1.length;t++, g++,l++){ // cycle through all genes
				int count=0;
				String whichfileans = new StringBuffer().append(datesh1n1[0]).append(".").append(genesh1n1[g]).append(".").append(typesh1n1[t]).append(".").append("nex").toString();
				File fleans = new File("/Users/sam/Desktop/FluData/SplitSingle/"+ whichfileans);
				ansname = fleans.toString();
				Read_main ra = new Read_main(ansname); // read ans
				int[][] tmp = ra.readNEXUS();
				int[] ans = ra.consensusArray(tmp);			
				int[] sampler = new int[ans.length/3];
				Random generator = new Random();
				for(int x=0;x<sampler.length;x++){
					int randint = generator.nextInt(sampler.length-1);			
					sampler[x] = randint;
				}
				if(t==8){ // subset for NA list surface
					int numsites=0;
					for(int j=0;j<Nd.badlisth1n1.length;j++){if(Nd.badlisth1n1[j]==1){numsites++;}}
					sampler = new int[(numsites)];
					Random generator2 = new Random();
					for(int x=0;x<sampler.length;x++){
						int randint = generator2.nextInt(sampler.length-1);			
						sampler[x] = randint;
					}
				}
				if(t==9){ // subset for NA list internal
					int numsites=0;
					for(int j=0;j<Nd.badlisth1n1.length;j++){if(Nd.badlisth1n1[j]==0){numsites++;}}
					sampler = new int[(numsites)];
					Random generator2 = new Random();
					for(int x=0;x<sampler.length;x++){
						int randint = generator2.nextInt(sampler.length-1);			
						sampler[x] = randint;
					}
				}

				for(int d=0;d<datesh1n1.length;d++){	// cycle through all time points			

					String whichfile = new StringBuffer().append(datesh1n1[d]).append(".").append(genesh1n1[g]).append(".").append(typesh1n1[t]).append(".").append("nex").toString();
					File fle = new File("/Users/sam/Desktop/FluData/SplitSingle/"+ whichfile);
					// file reading ***********************************************************
					if(fle.exists()){
						// use first time point as ans
						count++;
						Read_main r = new Read_main(fle.toString());
						int[][] seq = r.readNEXUS();	


						//na block *************
						if(t==8){
							Methods m = new Methods();
							int[][] seqNA = m.Subsetter(seq, Nd.badlisth1n1, 1);
							int[] ansNA = m.Subsetter(ans, Nd.badlisth1n1, 1);
							Williamson3bin ww = new Williamson3bin(seqNA,ansNA); //******
							Store s = ww.CreateBlocks(3,seqNA[0].length,sampler); //******
							Williamson3bin ws = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							ws.williamson3bin_method(nr[t],low,mid,high);

							if(run==0 && d==0){
								woc.Low_S.add(Double.valueOf(0.0));
								woc.Low_R.add(Double.valueOf(0.0));
								woc.Mid_S.add(Double.valueOf(0.0));
								woc.Mid_R.add(Double.valueOf(0.0));
								woc.High_S.add(Double.valueOf(0.0));
								woc.High_R.add(Double.valueOf(0.0));
								woc.dateavg.add(dateA[d]);
								woc.LowA.add(Double.valueOf(0.0));
								woc.MidA.add(Double.valueOf(0.0));
								woc.HighA.add(Double.valueOf(0.0));
								woc.FixA.add(Double.valueOf(0.0));
								woc.TotalA.add(Double.valueOf(0.0));
								woc.numSamples.add(seq.length);
								woc.neut=nr[t];
								woc.L=ans.length;
							}
							if(run==0){
								woc.Low_S.add(ws.low_S);
								woc.Low_R.add(ws.low_R);
								woc.Mid_S.add(ws.mid_S);
								woc.Mid_R.add(ws.mid_R);
								woc.High_S.add(ws.high_S);
								woc.High_R.add(ws.high_R);
								woc.dateavg.add(dateA[d]);
								woc.LowA.add(ws.low_A);
								woc.MidA.add(ws.mid_A);
								woc.HighA.add(ws.high_A);
								woc.FixA.add(ws.fix_A);
								woc.TotalA.add(ws.Adapt);
								woc.numSamples.add(seq.length);
								woc.neut=nr[t];
								woc.L=ans.length;
							}
							woc.BS.get(d)[run]=ws.Adapt;
						} else if (t==9){
							Methods m = new Methods();
							int[][] seqNA = m.Subsetter(seq, Nd.badlisth1n1, 0);
							int[] ansNA = m.Subsetter(ans, Nd.badlisth1n1, 0);
							Williamson3bin ww = new Williamson3bin(seqNA,ansNA); //******
							Store s = ww.CreateBlocks(3,seqNA[0].length,sampler); //******
							Williamson3bin ws = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							ws.williamson3bin_method(nr[t],low,mid,high);

							if(run==0 && d==0){
								woc.Low_S.add(ws.low_S);
								woc.Low_R.add(ws.low_R);
								woc.Mid_S.add(ws.mid_S);
								woc.Mid_R.add(ws.mid_R);
								woc.High_S.add(ws.high_S);
								woc.High_R.add(ws.high_R);
								woc.dateavg.add(dateA[d]);
								woc.LowA.add(ws.low_A);
								woc.MidA.add(ws.mid_A);
								woc.HighA.add(ws.high_A);
								woc.FixA.add(ws.fix_A);
								woc.TotalA.add(ws.Adapt);
								woc.numSamples.add(seq.length);
								woc.neut=nr[t];
								woc.L=ans.length;
							}
							if(run==0){
								woc.Low_S.add(ws.low_S);
								woc.Low_R.add(ws.low_R);
								woc.Mid_S.add(ws.mid_S);
								woc.Mid_R.add(ws.mid_R);
								woc.High_S.add(ws.high_S);
								woc.High_R.add(ws.high_R);
								woc.dateavg.add(dateA[d]);
								woc.LowA.add(ws.low_A);
								woc.MidA.add(ws.mid_A);
								woc.HighA.add(ws.high_A);
								woc.FixA.add(ws.fix_A);
								woc.TotalA.add(ws.Adapt);
								woc.numSamples.add(seq.length);
								woc.neut=nr[t];
								woc.L=ans.length;
							}
							woc.BS.get(d)[run]=ws.Adapt;
						} else {
							Williamson3bin ww = new Williamson3bin(seq,ans); //******
							Store s = ww.CreateBlocks(3,seq[0].length,sampler); //******
							Williamson3bin ws = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							ws.williamson3bin_method(nr[t],low,mid,high);
							if(run==0 && d==0){
								woc.Low_S.add(ws.low_S);
								woc.Low_R.add(ws.low_R);
								woc.Mid_S.add(ws.mid_S);
								woc.Mid_R.add(ws.mid_R);
								woc.High_S.add(ws.high_S);
								woc.High_R.add(ws.high_R);
								woc.dateavg.add(dateA[d]);
								woc.LowA.add(ws.low_A);
								woc.MidA.add(ws.mid_A);
								woc.HighA.add(ws.high_A);
								woc.FixA.add(ws.fix_A);
								woc.TotalA.add(ws.Adapt);
								woc.numSamples.add(seq.length);
								woc.neut=nr[t];
								woc.L=ans.length;
							}
							if(run==0){
								woc.Low_S.add(ws.low_S);
								woc.Low_R.add(ws.low_R);
								woc.Mid_S.add(ws.mid_S);
								woc.Mid_R.add(ws.mid_R);
								woc.High_S.add(ws.high_S);
								woc.High_R.add(ws.high_R);
								woc.dateavg.add(dateA[d]);
								woc.LowA.add(ws.low_A);
								woc.MidA.add(ws.mid_A);
								woc.HighA.add(ws.high_A);
								woc.FixA.add(ws.fix_A);
								woc.TotalA.add(ws.Adapt);
								woc.numSamples.add(seq.length);
								woc.neut=nr[t];
								woc.L=ans.length;
							}
							woc.BS.get(d)[run]=ws.Adapt;

						}
					}

				}
			}
			System.out.println("I am on Run  "+run+"  of 1000");
		}

		return woc;

	}


	public Value[][] clumpanalysish3n2(){		
		//	String[] datesh3n2 ={"1977","1978","1979","1980","1981","1982","1983","1984","1985","1986","1987","1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009"};
		//	String[] datesh3n2 =  {"1977_1979","1980_1982","1983_1985","1986_1988","1992_1994","1995_1997","1998_2000","2001_2003","2004_2006","2007_2009"};
		Neuraminidase Nd = new Neuraminidase();
		String[] datesh3n2 ={"1977","1978-1979-1980","1983-1984-1985","1986-1987-1988","1990-1991-1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009"};  // new dates

		String[] typesh3n2 = {"H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2"};
		String[] genesh3n2 = {"PB2","PB1","PA","HA","HA1","HA2","NP","NA","NA","NA","M1","M2","NS1"};
		double[] nr00={0.337528771,0.365306046,0.362118466,0.431014151,0.431014151,0.431014151,0.278018168,0.314742237,0.314742237,0.314742237,0.199342105,1.047104247,1.09495549};
		double[] nr05={0.140006516,0.126206897,0.155296629,0.21875,0.21875,0.21875,0.273537755,0.171428571,0.171428571,0.171428571,0.193353474,0.382022472,0.8};
		double[] nr10={0.102505695,0.100700848,0.16549559,0.183908046,0.183908046,0.183908046,0.343782929,0.136170213,0.136170213,0.136170213,0.201005025,0.192307692,0.962962963};
		double[] nr15={0.095502779,0.098823529,0.175540327,0.172248804,0.172248804,0.172248804,0.341470588,0.112149533,0.112149533,0.112149533,0.208092486,0.227272727,1};
		double[] nr20={0.113621262,0.08988764,0.170897833,0.225,0.225,0.225,0.294510908,0.135135135,0.135135135,0.135135135,0.171428571,0.357142857,0.722222222};
		double[] nr25={0.071672355,0.062593145,0.158385093,0.258064516,0.258064516,0.258064516,0.365313653,0.036363636,0.036363636,0.036363636,0.142857143,0.25,0.818181818};
		double[] nr30={0.066176471,0.073043478,0.153681964,0.259259259,0.259259259,0.259259259,0.33496572,0.042553191,0.042553191,0.042553191,0.166666667,0.375,0.777777778};
		double[] nr35={0.069230769,0.067542214,0.15982242,0.269230769,0.269230769,0.269230769,0.272409779,0.044444444,0.044444444,0.044444444,0.176470588,0.5,0.75};
		double[] nr40={0.098092643,0.056470588,0.183079057,0.318181818,0.318181818,0.318181818,0.266846361,0.066666667,0.066666667,0.066666667,0.176470588,0.75,2};
		double[] nr45={0.102739726,0.054830287,0.170616114,0.358974359,0.358974359,0.358974359,0.275766017,0.074074074,0.074074074,0.074074074,0.214285714,1,1.5};	


		double[] nr = {0.115915692,0.095906121,0.192656532,0.160427807,0.160427807,0.160427807,0.25374128,0.148514851,0.148514851,0.148514851,0.307692308,0.3,0.722891566};
		String ansname = new String("");

		Value[][] value_matrix = new Value[datesh3n2.length][typesh3n2.length];
		for(int i=0;i<value_matrix.length;i++){
			for(int j=0;j<value_matrix[0].length;j++){
				value_matrix[i][j] = new Value(1000);
			}
		}	
		for(int t=0,g=0,l=0;t<genesh3n2.length;t++, g++,l++){
			int count=0;
			for(int d=0;d<datesh3n2.length;d++){
				String whichfile = new StringBuffer().append(datesh3n2[d]).append(".").append(genesh3n2[g]).append(".").append(typesh3n2[t]).append(".").append("nex").toString();
				File fle = new File("/Users/sam/Desktop/FluData/SplitSingle2/"+ whichfile);
				// file reading ***********************************************************
				// use first time point as ans
				if(fle.exists()){
					count++;
					if(count==1){
						ansname = fle.toString();
					}			
					Read_main ra = new Read_main(ansname);
					int[][] tmp = ra.readNEXUS();
					int[] ans = ra.consensusArray(tmp);
					//	int sequenceNum=num;  //max12
					//	for(int i=0;i<tmp[0].length;i++){
					//		ans[i] = tmp[sequenceNum][i];
					//	}


					Read_main r = new Read_main(fle.toString());
					int[][] seq = r.readNEXUS();	
					//na block *************
					if(t==8){
						Methods m = new Methods();
						int[][] seqNA = m.Subsetter(seq, Nd.badlisth3n2, 1);
						int[] ansNA = m.Subsetter(ans, Nd.badlisth3n2, 1);
						Williamson3bin w = new Williamson3bin(seqNA,ansNA);
						w.williamson3bin_method(nr[t],low,mid,high);
						value_matrix[d][t].Nadapt = w.Adapt;
						//		value_matrix[d][t].Nadapt = w.mid_R;
						//		value_matrix[d][t].Nadapt = w.L;
						value_matrix[d][t].row = datesh3n2[d];
						value_matrix[d][t].column = genesh3n2[t];
						value_matrix[d][t].dataset = "H3N2";
					} else if (t==9){
						Methods m = new Methods();
						int[][] seqNA = m.Subsetter(seq, Nd.badlisth3n2, 0);
						int[] ansNA = m.Subsetter(ans, Nd.badlisth3n2, 0);
						Williamson3bin w = new Williamson3bin(seqNA,ansNA);
						w.williamson3bin_method(nr[t],low,mid,high);
						value_matrix[d][t].Nadapt =  w.Adapt;
						//		value_matrix[d][t].Nadapt = w.mid_R;
						//		value_matrix[d][t].Nadapt = w.L;
						value_matrix[d][t].row = datesh3n2[d];
						value_matrix[d][t].column = genesh3n2[t];
						value_matrix[d][t].dataset = "H3N2";
					} else {
						Williamson3bin w = new Williamson3bin(seq,ans);
						w.williamson3bin_method(nr[t],low,mid,high);
						value_matrix[d][t].Nadapt =  w.Adapt;
						//		value_matrix[d][t].Nadapt = w.mid_R;
						//		value_matrix[d][t].Nadapt = w.L;
						value_matrix[d][t].row = datesh3n2[d];
						value_matrix[d][t].column = genesh3n2[t];
						value_matrix[d][t].dataset = "H3N2";
					}
				}
			}
		}
		return value_matrix;
	}

	public Value[][] clumpanalysish3n2BS(int N){		
		//	String[] datesh3n2 =  {"1977_1979","1980_1982","1983_1985","1986_1988","1992_1994","1995_1997","1998_2000","2001_2003","2004_2006","2007_2009"};
		String[] datesh3n2 ={"1977","1978-1979-1980","1983-1984-1985","1986-1987-1988","1990-1991-1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009"};  // new dates

		String ansname = new String("");
		Neuraminidase Nd = new Neuraminidase();
		String[] typesh3n2 = {"H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2","H3N2"};
		String[] genesh3n2 = {"PB2","PB1","PA","HA","HA1","HA2","NP","NA","NA","NA","M1","M2","NS1"};
		double[] nr = {0.115915692,0.095906121,0.192656532,0.160427807,0.160427807,0.160427807,0.25374128,0.148514851,0.148514851,0.148514851,0.307692308,0.3,0.722891566};

		Value[][] value_matrix = new Value[datesh3n2.length][typesh3n2.length];
		for(int i=0;i<value_matrix.length;i++){
			for(int j=0;j<value_matrix[0].length;j++){
				value_matrix[i][j] = new Value(N);
			}
		}	
		for(int run=0;run<N;run++){

			for(int t=0,g=0,l=0;t<genesh3n2.length;t++, g++,l++){
				int count=0;
				String whichfileans = new StringBuffer().append(datesh3n2[0]).append(".").append(genesh3n2[g]).append(".").append(typesh3n2[t]).append(".").append("nex").toString();
				File fleans = new File("/Users/sam/Desktop/FluData/SplitSingle2/"+ whichfileans);
				ansname = fleans.toString();
				Read_main ra = new Read_main(ansname);
				int[][] tmp = ra.readNEXUS();
				int[] ans = ra.consensusArray(tmp);			
				int[] sampler = new int[ans.length/3];
				Random generator = new Random();
				for(int x=0;x<sampler.length;x++){
					int randint = generator.nextInt(sampler.length);//generator.nextInt(sampler.length-1);
					sampler[x] = randint;
				}
				if(t==8){
					int numsites=0;
					for(int j=0;j<Nd.badlisth3n2.length;j++){if(Nd.badlisth3n2[j]==1){numsites++;}}
					sampler = new int[(numsites)];
					Random generator2 = new Random();
					for(int x=0;x<sampler.length;x++){
						int randint = generator2.nextInt(sampler.length);//generator2.nextInt(sampler.length-1);
						sampler[x] = randint;
					}
				}
				if(t==9){
					int numsites=0;
					for(int j=0;j<Nd.badlisth3n2.length;j++){if(Nd.badlisth3n2[j]==0){numsites++;}}
					sampler = new int[(numsites)];
					Random generator2 = new Random();
					for(int x=0;x<sampler.length;x++){
						int randint = generator2.nextInt(sampler.length); //generator2.nextInt(sampler.length-1);
						sampler[x] = randint;
					}
				}


				for(int d=0;d<datesh3n2.length;d++){

					String whichfile = new StringBuffer().append(datesh3n2[d]).append(".").append(genesh3n2[g]).append(".").append(typesh3n2[t]).append(".").append("nex").toString();
					File fle = new File("/Users/sam/Desktop/FluData/SplitSingle2/"+ whichfile);
					// file reading ***********************************************************
					if(fle.exists()){
						// use first time point as ans
						count++;
						Read_main r = new Read_main(fle.toString());
						int[][] seq = r.readNEXUS();	


						//na block *************
						if(t==8){
							Methods m = new Methods();
							int[][] seqNA = m.Subsetter(seq, Nd.badlisth3n2, 1);
							int[] ansNA = m.Subsetter(ans, Nd.badlisth3n2, 1);
							Williamson3bin ww = new Williamson3bin(seqNA,ansNA); //******
							Store s = ww.CreateBlocks(3,seqNA[0].length,sampler); //******
							Williamson3bin w = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							w.williamson3bin_method(nr[t],low,mid,high);
							value_matrix[d][t].row = datesh3n2[d];
							value_matrix[d][t].column = genesh3n2[g];
							value_matrix[d][t].dataset = "H3N2";
							value_matrix[d][t].Bstrap[run] = w.Adapt;
						} else if (t==9){
							Methods m = new Methods();
							int[][] seqNA = m.Subsetter(seq, Nd.badlisth3n2, 0);
							int[] ansNA = m.Subsetter(ans, Nd.badlisth3n2, 0);
							Williamson3bin ww = new Williamson3bin(seqNA,ansNA); //******
							Store s = ww.CreateBlocks(3,seqNA[0].length,sampler); //******
							Williamson3bin w = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							w.williamson3bin_method(nr[t],low,mid,high);
							value_matrix[d][t].row = datesh3n2[d];
							value_matrix[d][t].column = genesh3n2[g];
							value_matrix[d][t].dataset = "H1N1";
							value_matrix[d][t].Bstrap[run] = w.Adapt;
						} else {
							Williamson3bin ww = new Williamson3bin(seq,ans); //******
							Store s = ww.CreateBlocks(3,seq[0].length,sampler); //******
							Williamson3bin w = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
							w.williamson3bin_method(nr[t],low,mid,high);
							value_matrix[d][t].row = datesh3n2[d];
							value_matrix[d][t].column = genesh3n2[g];
							value_matrix[d][t].dataset = "H1N1";
							value_matrix[d][t].Bstrap[run] = w.Adapt;	
						}

					}

				}
			}
			System.out.println("I am on Run  "+run+"  of 1000");
		}
		return value_matrix;
	}








}



