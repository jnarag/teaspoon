package teaspoon;

import teaspoon.adaptation.Hemagglutinin;
import teaspoon.adaptation.Methods;
import teaspoon.adaptation.Neuraminidase;
import teaspoon.adaptation.WilliamsonOutputClass;

public class FluB {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Methods m = new Methods();
		m.SubsetAlignment("/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/cleaned_LL1404_B_FullHA_Victoria.fas", "/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/External_Victoria.fas", "/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/Internal_Victoria.fas", new Hemagglutinin().FluB );
		m.SubsetAlignment("/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/cleaned_LL1404_B_FullHA_Yamagata.fas", "/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/External_Yamagata.fas", "/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/Internal_Yamagata.fas", new Hemagglutinin().FluB );
		m.SubsetAlignment("/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/pure_Victoria_FluB_NA.fas", "/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/External_Victoria_NA.fas", "/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/Internal_Victoria_NA.fas", new Neuraminidase().FluB );
		m.SubsetAlignment("/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/pure_Yamagata_FluB_NA.fas", "/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/External_Yamagata_NA.fas", "/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/Internal_Yamagata_NA.fas", new Neuraminidase().FluB );

		/*		int start=1971;
		int end=2014;
		String[] files = {"cleaned_LL535_B_HA2_combined.fas",
				"cleaned_LL535_B_HA2_Victoria.fas",
				"cleaned_LL535_B_HA2_Yamagata.fas",
				"cleaned_LL868_B_HA1_combined.fas",
				"cleaned_LL868_B_HA1_Victoria.fas",
				"cleaned_LL868_B_HA1_Yamagata.fas",
				"cleaned_LL1404_B_FullHA_combined.fas",
				"cleaned_LL1404_B_FullHA_Victoria.fas",
				"cleaned_LL1404_B_FullHA_Yamagata.fas",
				"pure_Victoria_FluB_NA.fas",
				"pure_Yamagata_FluB_NA.fas"};
		int size = end-start;
		int[][] store = new int[size][files.length];
		for(int x=0;x<files.length;x++){
		    teaspoon.adaptation.DataSet data = new teaspoon.adaptation.DataSet("/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/"+files[x]);
			int[][] counter = data.dateCounter(start, end);

			for(int i=0;i<counter.length;i++){
				store[i][x]=counter[i][1];
			}		

		}

			for(int j=0;j<size;j++){
				for(int i=0;i<files.length;i++){
				System.out.print(store[j][i]+"\t");
			}
			System.out.println();	

		}*/



		//victoria HA1
		String[] DatesVHA1={"1975_1976_1977_1978_1979_1980","1982_1983_1984_1985","1986_1987","1988_1989","1990","1991_1992_1993_1994","1995_1996_1997","1998_1999_2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012"};
		//victoria HA1
		String[] DatesVHA2={"1985_1986","1987_1988_1989_1990","1993_1994_1995_1996_1997","1998_1999_2000_2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011_2012"};
		//Yamagataha1
		String[] DatesYHA1={"1975_1976_1977_1978_1979_1980_1981","1982_1983_1984_1985","1987_1988_1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012"};
		//Yamagata ha2
		String[] DatesYHA2={"1979_1980_1981_1982_1983","1988_1989_1990_1991","1992_1993","1994_1995","1996","1997","1998","1999","2000","2001","2002_2003","2004","2005","2006","2007","2008","2009","2010","2011_2012"};

		String[] DatesVHA = {"1985_1986","1987_1988_1989_1990","1993_1994_1995_1996_1997_1998","2001_2002","2003","2004","2005","2006","2007","2008","2009","2010","2011_2012"};

		String[] DatesYHA = {"1979_1980_1981_1982_1983","1988_1989_1990_1991","1992_1993","1994_1995","1996","1997","1998","1999","2000","2001","2002_2003","2004","2005","2006","2007","2008","2009","2010","2011_2012"};

		String[] DatesVNA = {"2002","2003","2004","2005","2006","2007","2008","2009","2010","2011"};

		String[] DatesYNA = {"1988_1989_1990","1991_1992_1993","1994_1995","1996","1997","1998","1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011_2012"};
		String gene="NA";

		if(gene=="HA"){
			String choice="V";
			int indexHA1=-9;
			int indexHA2=-9;
			int indexHA=-9;
			String[] DatesHA1=DatesVHA1;
			String[] DatesHA2=DatesVHA2;
			String[] DatesHA=DatesVHA;
			String[] files = {
					"Internal_Victoria.fas",
					"Internal_Yamagata.fas",
					"External_Victoria.fas",
			"External_Yamagata.fas"};

			if(choice=="V"){
				indexHA1=2;
				indexHA2=0;
				DatesHA1=DatesVHA;
				DatesHA2=DatesVHA;
			}
			if(choice=="Y"){
				indexHA1=3;
				indexHA2=1;
				DatesHA1=DatesYHA;
				DatesHA2=DatesYHA;

			}


			FluAnalysis2 F = new FluAnalysis2();
			Hemagglutinin H = new Hemagglutinin();

			WilliamsonOutputClass woc = F.GenericRunner("/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/"+files[indexHA2], DatesHA2, 0.1,2);
			double nr = F.CalcNr(woc);
			System.out.println(nr);

			woc = F.GenericRunner("/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/"+files[indexHA1], DatesHA1, nr, 100);
			F.print(woc,10,"/Users/Sam/Desktop/VictoriaHA1.txt");


			for(int i=0;i<woc.TotalA.size();i++){
				System.out.println(woc.dateavg.get(i)+"\t"+woc.TotalA.get(i)+"\t"+woc.Mid_R.get(i)+"\t"+woc.Mid_S.get(i));
			}
		} else if(gene=="NA"){
			
			String choice="Y";
			int indexI=-9;
			int indexE=-9;
			String[] DatesNA=DatesVNA;
			String[] files = {
					"Internal_Victoria_NA.fas",
					"Internal_Yamagata_NA.fas",
					"External_Victoria_NA.fas",
					"External_Yamagata_NA.fas"};

			if(choice=="V"){
				indexE=2;
				indexI=0;
				DatesNA=DatesVNA;
				DatesNA=DatesVNA;
			}
			if(choice=="Y"){
				indexE=3;
				indexI=1;
				DatesNA=DatesYNA;
				DatesNA=DatesYNA;

			}


			FluAnalysis2 F = new FluAnalysis2();


			WilliamsonOutputClass woc = F.GenericRunner("/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/"+files[indexI], DatesNA, 0.1,1);
			
			double nr = F.CalcNr(woc);
			System.out.println(nr);

			woc = F.GenericRunner("/Users/Sam/Desktop/PAPERS/FLUB/B_ready_HA/"+files[indexE], DatesNA, nr, 1);
	//		F.print(woc,10,"/Users/Sam/Desktop/VictoriaHA1.txt");


			for(int i=0;i<woc.TotalA.size();i++){
				System.out.println(woc.dateavg.get(i)+"\t"+woc.TotalA.get(i)+"\t"+woc.Mid_R.get(i)+"\t"+woc.Mid_S.get(i));
			}
			
			
		}




	}
}



