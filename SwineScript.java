package teaspoon;

import teaspoon.adaptation.FluAnalysis;
import teaspoon.adaptation.Hemagglutinin;
import teaspoon.adaptation.Neuraminidase;
import teaspoon.adaptation.WilliamsonOutputClass;

public class SwineScript {

	
	public static void main(String[] args) {
		FluAnalysis2 F = new FluAnalysis2();
		int boot=10;
		int output=1;
		boolean printseq=true;
		 System.out.println("Starting Program - Gene: "+args[1]+" and Genotype: "+args[0]+"\n");	 
		 long time = System.currentTimeMillis();

		if(args[1].equals("PB2") || args[1].equals("PB1") || args[1].equals("PA")){
			String Name = args[1];  // specify name
			String Num="#9";
			if(Name.equals("PB2")){
				Num = "#1"; //type
			} else if(Name.equals("PB1")){
				Num = "#2"; //type
			} else if(Name.equals("PA")){
				Num = "#3"; //type
			}
			String Geno = args[0]; //type
			Neuraminidase Nd = new Neuraminidase();
			double nr=0.1;
			// create types
			String[] type = {Num+Geno,Num+Geno+"a",Num+Geno+"b",Num+Geno+"c",Num+Geno+"d",Num+Geno+"e",Num+Geno+"f",Num+Geno+"g",Num+Geno+"h"};
			//		F.Extract("/Users/Sam/Desktop/FluData/sw_ref_PB2.fasta", "/Users/Sam/Desktop/FluData/Table.csv", "PB2",type, "1999-2000", "/Users/Sam/Desktop/out.fa");
			if(Geno.equals("C1")){
				// define dayes
//				String[] dates={"1975","1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1998-1999-2000","2001","2002-2003","2004-2005-2006","2008-2009"};			
				String[] dates={"1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1998-1999-2000","2001","2002-2003","2004-2005-2006","2008-2009"};			

				// run to calculate neutral ratio
				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0);
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0,Geno);
				}
				double neutralratio = F.CalcNr(woc);  //save neutral ratio 
				// run again to get restuls with this neutral ratio
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_");
			} else if(Geno.equals("B1")){
//				String[] dates={"1984","1991-1992-1993","1994-1995","1996-1997","1998-1999-2000","2001","2002","2003","2004","2005","2006","2007-2008","2009-2010"};
				String[] dates={"1987","1991-1992-1993","1994-1995","1996-1997","1998-1999-2000","2001","2002","2003","2004","2005","2006","2007-2008","2009-2010"};

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0);
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0,Geno);
				}
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_");
			} else if(Geno.equals("T1")){
				String[] dates={"1999","2000-2001-2003","2003","2004","2005","2006","2007","2008","2009-2010"};
				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0);
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0,Geno);
				}
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_");
			}
		}

		//********************/********************/********************/********************/********************/********************/********************/********************
		//HA HA1 HA2 teaspoon.adaptation.Split
		if(args[1].equals("H1")){

			String Name = args[1];  // specify name
			String Num = "#4"; //type
			String Geno = args[0]; //type
			int A=10;
			Hemagglutinin Hd = new Hemagglutinin();  /// important needed to change lists
			// create types
			String[] type = {Num+Geno,Num+Geno+"a",Num+Geno+"b",Num+Geno+"c",Num+Geno+"d",Num+Geno+"e",Num+Geno+"f",Num+Geno+"g",Num+Geno+"h"};

			if(args[0].equals("B1")){
				String[] dates ={"1987","1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001-2002","2003-2004","2005","2006-2007","2008","2009-2010"}; //ea
//				String[] dates ={"1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001-2002","2003-2004","2005","2006-2007","2008","2009-2010"}; //ea

				double nr = 0.1; // EA

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.HAH1,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Hd.HAH1,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HA2"+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.HAH1,0,"B1_HA2");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.HAH1,1,"B1_HA1");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Hd.HAH1,0,"B1_HA");

				}


				WilliamsonOutputClass woc3 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, neutralratio,boot,Hd.HAH1,0);
				F.print(woc3, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HA"+"_");


				WilliamsonOutputClass woc4 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Hd.HAH1,1);
				F.print(woc4, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HA1"+"_");

			//	String[] dateswa ={"1979","1980-1981-1982-1983","1984-1985-1986-1987","1990-1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001-2002","2003-2004","2005","2006-2007","2008","2009-2010"}; //eawa
			//	String[] dateswa ={"1980-1981-1982-1983","1984-1985-1986-1987","1990-1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001-2002","2003-2004","2005","2006-2007","2008","2009-2010"}; //eawa

			//	teaspoon.adaptation.WilliamsonOutputClass woc5 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+"_wa.fasta", "/Users/Sam/Desktop/FluData/Table_H1.csv",Name, type,dateswa,0, -2, neutralratio,boot,Hd.HAH1,1);
			//	F.print(woc5, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HA1_With Additions"+"_");
			} else if(args[0].equals("C1")){
		//		String[] dates ={"1975","1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1999","2000","2001","2002","2003","2004","2005","2006-2007","2008","2009-2010"}; //cs
				String[] dates ={"1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1999-2000","2001","2002","2003","2004","2005","2006-2007","2008","2009-2010"}; //cs

				double nr = 0.1; // EA

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.HAH1,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Hd.HAH1,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HA2"+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.HAH1,0,"C1_HA2");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.HAH1,1,"C1_HA1");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Hd.HAH1,0,"C1_HA");

				}

				WilliamsonOutputClass woc3 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, neutralratio,boot,Hd.HAH1,0);
				F.print(woc3, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HA"+"_");


				WilliamsonOutputClass woc4 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Hd.HAH1,1);
				F.print(woc4, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HA1"+"_");


			}



		}
		//********************/********************/********************/********************/********************/********************/********************/********************
		//HA surface internal split
		if(args[1].equals("H1si")){

			String Name = "H1";  // specify name
			String Num = "#4"; //type
			String Geno = args[0]; //type
			int A=10;
			Hemagglutinin Hd = new Hemagglutinin();  /// important needed to change lists
			// create types
			String[] type = {Num+Geno,Num+Geno+"a",Num+Geno+"b",Num+Geno+"c",Num+Geno+"d",Num+Geno+"e",Num+Geno+"f",Num+Geno+"g",Num+Geno+"h"};
	//		String[] type = {Num+Geno,Num+Geno+"a"};

			if(args[0].equals("B1")){
		//		String[] dates ={"1984","1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001-2002","2003-2004","2005","2006-2007","2008","2009-2010"}; //ea
				String[] dates ={"1987","1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001-2002","2003-2004","2005","2006-2007","2008","2009-2010"}; //ea

	//			String[] dates ={"1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001-2002","2003-2004","2005","2006-2007","2008","2009-2010"}; //ea

				double nr = 0.1; // EA

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.H1,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Hd.H1,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HAi"+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.H1,0,"B1_HAi");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.H1,1,"B1_HAs");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Hd.H1,0,"B1_HA");

				}

				WilliamsonOutputClass woc3 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, neutralratio,boot,Hd.H1,0);
				F.print(woc3, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HAsi"+"_");


				WilliamsonOutputClass woc4 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Hd.H1,1);
				F.print(woc4, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HAs"+"_");

			//	String[] dateswa ={"1979","1980-1981-1982-1983","1984-1985-1986-1987","1990-1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001-2002","2003-2004","2005","2006-2007","2008","2009-2010"}; //eawa
			//	String[] dateswa ={"1980-1981-1982-1983","1984-1985-1986-1987","1990-1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001-2002","2003-2004","2005","2006-2007","2008","2009-2010"}; //eawa

			//	teaspoon.adaptation.WilliamsonOutputClass woc5 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+"_wa.fasta", "/Users/Sam/Desktop/FluData/Table_H1.csv",Name, type,dateswa,0, -2, neutralratio,boot,Hd.H1,1);
			//	F.print(woc5, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HAs_With Additions"+"_");
			} else if(args[0].equals("C1")){
		//		String[] dates ={"1975","1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1999","2000","2001","2002","2003","2004","2005","2006-2007","2008","2009-2010"}; //cs
				String[] dates ={"1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1999","2000","2001","2002","2003","2004","2005","2006-2007","2008","2009-2010"}; //cs

				double nr = 0.1; // EA

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.H1,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Hd.H1,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HAi"+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.H1,0,"C1_HAi");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Hd.H1,1,"C1_HAs");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Hd.H1,0,"C1_HA");

				}

				WilliamsonOutputClass woc3 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, neutralratio,boot,Hd.H1,0);
				F.print(woc3, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HAsi"+"_");


				WilliamsonOutputClass woc4 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Hd.H1,1);
				F.print(woc4, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_HAs"+"_");


			}



		}
		//********************/********************/********************/********************/********************/********************/********************/********************

		
		
		if(args[1].equals("NP")){

			//NP
			String Name = args[1];  // specify name
			String Num = "#5"; //type
			String Geno = args[0]; //type
			Neuraminidase Nd = new Neuraminidase();
			
			// create types
			String[] type = {Num+Geno,Num+Geno+"a",Num+Geno+"b",Num+Geno+"c",Num+Geno+"d",Num+Geno+"e",Num+Geno+"f",Num+Geno+"g",Num+Geno+"h"};
			double nr=0.1;
			if(args[0].equals("B1")){
//				String[] dates ={"1984","1991-1992-1993","1994-1995-1996","1997-1998","1999","2000-2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};
				String[] dates ={"1987","1991-1992-1993","1994-1995-1996","1997-1998","1999","2000-2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_NP.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_NP.csv",Name, type,dates,0, -1, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0,Geno);
				}
			}else if(args[0].equals("C1")){
		//		String[] dates ={"1975","1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};
				String[] dates ={"1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_NP.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_NP.csv",Name, type,dates,0, -1, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0,Geno);
				}
			}
		}


		//********************/********************/********************/********************/********************/********************/********************/********************

		if(args[1].equals("N1")){
			String Name = args[1];  // specify name
			String Num = "#6"; //type
			String Geno = args[0]; //type
			Neuraminidase Nd = new Neuraminidase();
			// create types
			String[] type = {Num+Geno,Num+Geno+"a",Num+Geno+"b",Num+Geno+"c",Num+Geno+"d",Num+Geno+"e",Num+Geno+"f",Num+Geno+"g",Num+Geno+"h"};
			if(args[0].equals("B1")){
//				String[] dates = {"1984","1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001","2002-2003","2004","2005","2006-2007","2008","2009"};
				String[] dates = {"1987","1991-1992-1993","1994-1995-1996","1997-1998-1999","2000-2001","2002-2003","2004","2005","2006-2007","2008","2009"};

				double nr = 0.1; // EA

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Nd.badlisth1n1_NoStopCodon,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_NAi"+"_");


				WilliamsonOutputClass woc3 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,0);
				F.print(woc3, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_NA"+"_");


				WilliamsonOutputClass woc4 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,1);
				F.print(woc4, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_NAs"+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Nd.badlisth1n1_NoStopCodon,0,"B1_NAi");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Nd.badlisth1n1_NoStopCodon,1,"B1_NAs");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0,"B1_NA");

				}

				
			}else if(args[0].equals("C1")){
//				String[] dates = {"1975","1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1999-2000","2001","2002-2003","2004-2005","2006-2007","2008-2009"};
				String[] dates = {"1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1999-2000","2001","2002-2003","2004-2005","2006-2007","2008-2009"};

				double nr = 0.1; // EA

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Nd.badlisth1n1_NoStopCodon,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_NAi"+"_");


				WilliamsonOutputClass woc3 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,0);
				F.print(woc3, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_NA"+"_");


				WilliamsonOutputClass woc4 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, neutralratio,boot,Nd.badlisth1n1_NoStopCodon,1);
				F.print(woc4, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_NAs"+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Nd.badlisth1n1_NoStopCodon,0,"C1_NAi");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -2, nr,1,Nd.badlisth1n1_NoStopCodon,1,"C1_NAs");
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Nd.badlisth1n1_NoStopCodon,0,"C1_NA");

				}
				
			}
		}
		//********************/********************/********************/********************/********************/********************/********************/********************

		if(args[1].equals("M1") || args[1].equals("M2")  ){
			//M1 M2
			String Name = args[1];  // specify name
			String Num = "#7"; //type
			String Geno = args[0]; //type
			Hemagglutinin Hd = new Hemagglutinin();
			// create types
			String[] type = {Num+Geno,Num+Geno+"a",Num+Geno+"b",Num+Geno+"c",Num+Geno+"d",Num+Geno+"e",Num+Geno+"f",Num+Geno+"g",Num+Geno+"h"};
			double nr=0.1;
			if(args[0].equals("C1")){
//				String[] dates ={"1975","1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1999-2000","2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};
				String[] dates ={"1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1999-2000","2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};
			
				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_MP.csv",Name, type,dates,0, -1, nr,1,Hd.HA,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_MP.csv",Name, type,dates,0, -1, neutralratio,boot,Hd.HA,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Hd.HA,0,Geno);
				}

			}else 	if(args[0].equals("B1")){
//				String[] dates ={"1984","1991-1992-1993","1994-1995-1996","1997-1998","1999","2000-2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};
				String[] dates ={"1987","1991-1992-1993","1994-1995-1996","1997-1998","1999","2000-2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_MP.csv",Name, type,dates,0, -1, nr,1,Hd.HA,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_MP.csv",Name, type,dates,0, -1, neutralratio,boot,Hd.HA,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Hd.HA,0,Geno);
				}
			}
		}

		//********************/********************/********************/********************/********************/********************/********************/********************
		if(args[1].equals("NS1") || args[1].equals("NS2")  ){
			//M1 M2
			String Name = args[1];  // specify name
			String Num = "#8"; //type
			String Geno = args[0]; //type
			Hemagglutinin Hd = new Hemagglutinin();
			// create types
			String[] type = {Num+Geno,Num+Geno+"a",Num+Geno+"b",Num+Geno+"c",Num+Geno+"d",Num+Geno+"e",Num+Geno+"f",Num+Geno+"g",Num+Geno+"h"};
			double nr=0.1;
			if(args[0].equals("C1")){
		// 		String[] dates ={"1975","1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1998-1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};
		 		String[] dates ={"1976","1977","1978","1979","1980-1981","1985-1986-1987","1988-1990-1991","1993","1994","1998-1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_NS.csv",Name, type,dates,0, -1, nr,1,Hd.HA,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_NS.csv",Name, type,dates,0, -1, neutralratio,boot,Hd.HA,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Hd.HA,0,Geno);
				}

			}else 	if(args[0].equals("B1")){
	//	 		String[] dates ={"1984","1991-1992-1993","1994-1995-1996","1997-1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};
		 		String[] dates ={"1987","1991-1992-1993","1994-1995-1996","1997-1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009-2010"};

				WilliamsonOutputClass woc = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_NS.csv",Name, type,dates,0, -1, nr,1,Hd.HA,0);
				double neutralratio = F.CalcNr(woc);
				WilliamsonOutputClass woc2 = F.RunnerBS("/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table_NS.csv",Name, type,dates,0, -1, neutralratio,boot,Hd.HA,0);
				F.print(woc2, output,"/Users/Sam/Desktop/out/"+Name+"_"+Geno+"_");
				if(printseq){
					F.WriteSequencesGeno("/Users/Sam/Desktop/out2/Swine","/Users/Sam/Desktop/FluData/sw_ref_"+Name+".fasta", "/Users/Sam/Desktop/FluData/Table.csv",Name, type,dates,0, -1, nr,1,Hd.HA,0,Geno);
				}
			}
		}
		//********************/********************/********************/********************/********************/********************/********************/********************
		if(args[1].equals("H1N1") || args[1].equals("H3N2") ){
			FluAnalysis ff = new FluAnalysis();
			double[] nr = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
			String[] names =  {"PB2","PB1","PA","HAsi","HAs","HAi","HA","HA1","HA2","NP","NA","NAs","NAi","M1","M2","NS1"};
			WilliamsonOutputClass[] woc = ff.HumanRunner(1, nr,args[1]);
			int k=0;
			for(WilliamsonOutputClass i :woc){
				nr[k]=F.CalcNr(i);
				k++;
			}
			nr[3]=nr[5];nr[4]=nr[5];
			nr[6]=nr[8];nr[7]=nr[8];
			nr[10]=nr[12];nr[11]=nr[12];	
			WilliamsonOutputClass[] woc2 = ff.HumanRunner(boot, nr,args[1]);
			k=0;
			for(WilliamsonOutputClass i :woc2){
				F.print(i, output, "/Users/Sam/Desktop/out/"+names[k]+"_"+args[1]+"_");
				k++;
				
			}
		}
		
		
		

		 time = System.currentTimeMillis() - time;
		 System.out.println("Program - Gene: "+args[1]+" and Genotype: "+args[0]+" took "+ time/1000 + " seconds");
	}

}
