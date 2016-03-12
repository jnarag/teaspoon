package teaspoon.adaptation;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;


public class PBtest {

	public final String[] AA =	{"K","N","K","N","T","T","T","T","R","S","R","S","I","I","M","I","Q","H","Q","H","P","P","P","P",
			"R","R","R","R","L","L","L","L","E","D","E","D","A","A","A","A","G","G","G","G","V","V","V","V",
			"X","Y","X","Y","S","S","S","S","X","C","W","C","L","F","L","F","?","-","?" };

	public PBtest(){

	}


	public static void main(String[] args) {


		// START ****************************************************************************************************************************************


		String sss = new String("/Users/sam/Desktop/tmp/");
		File input = new File(sss+"seqHCV");		// Input 
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
		File input2 = new File(sss+"ansHCV");		// Input 
		File[] list2 = input2.listFiles();
		ArrayList<String> data2 = new ArrayList<String>();
		for(int i=0;i<list2.length;i++){
			if(list2[i].isHidden() == false && list2[i].isDirectory() == false)	{
				data2.add(list2[i].getAbsolutePath());
			}
		}	
		Collections.sort(data2);

		// tracking code
	//	teaspoon.adaptation.Read_main re1 = new teaspoon.adaptation.Read_main(data2.get(1));
	//	int[][] tmp1 = re1.readFASTA();
	//	int[] ans1 = re1.consensusArray(tmp1);	

	//	ArrayList<teaspoon.adaptation.Mutation> Fstore = new ArrayList<teaspoon.adaptation.Mutation>();

	//	for(int i=0;i<ans1.length;i++){
	//		for(int k=0;k<4;k++){	
	//			if(k+1!=ans1[i]){
	//				teaspoon.adaptation.Mutation M = new teaspoon.adaptation.Mutation();
	//				M.base=k+1;
	//				M.site=i;
	//				Fstore.add(M);
	//			}
	//		}
	//	}

		// diag
		//	int y=0;
		/*	for(int i=0;i<ans1.length;i++){
			System.out.print(ans1[i]+"\t");
			System.out.print(Fstore.get(y).base+"\t");
			System.out.print(Fstore.get(y+1).base+"\t");
			System.out.print(Fstore.get(y+2).base+"\t");
			System.out.println();
			y=y+3;
		}*/


		for(int x=0;x<data.size();x++){ /// loop

			Read_main ra = new Read_main(data.get(x));
			Read_main re = new Read_main(data2.get(x));


	//		int[][] tmp = re.readFASTA();
	//		int[] ans = re.consensusArray(tmp);	
	//		int[][] seq = ra.readFASTA();

				int[][] tmp = re.read();
				int[] ans = re.consensusArray(tmp);	
				int[][] seq = ra.read();



	//		teaspoon.adaptation.BhattMethod bm = new teaspoon.adaptation.BhattMethod(seq,ans);
	//		double[] prior = {1,1,1,1};

	//		Fstore = bm.Tracking(Fstore, prior, true);
	//		System.out.println(Fstore.size());

					Methods m = new Methods();

			double[] L = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
			double[] H = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

			//		double[] L = {0.0};
			//		double[] H = {1.0};


			boolean[] Nvec = {false,true,true,true,false,false,false,false,false,false};


			double[][] bins = new double[2][L.length];
			for(int i=0;i<L.length;i++){
				bins[0][i]=L[i];
				bins[1][i]=H[i];
			}
			double[] prior = {1,1,1,1};


			BhattMethod bm = new BhattMethod(seq,ans);
		//	bm.MethodNG(bins,prior,true,Nvec);
				bm.Method(bins,prior,true,Nvec);
	//		bm.print(bm.NonNeutralSubstitutions);
	//		bm.print(bm.ReplacementSilentRatio);
	//		bm.print(bm.TotalCountArray);
			bm.print(bm.TotalCountArrayNoInvariant);
	//		bm.print(bm.ReplacementCountArray);	
				
		//	System.out.println(bm.neutralratio + "\t"+ bm.DeleteriousLoad +"\t"+ bm.Adaptation);

			 


		} 
		//  END

	/*	boolean[] flag = new boolean[Fstore.size()];
		
		
		for(int i=0;i<Fstore.size();i++){
			for(int h=1;h<data.size();h++){	
				if(Fstore.get(i).timeProb.get(h)>Fstore.get(i).H.get(h-1) || Fstore.get(i).timeProb.get(h)<Fstore.get(i).L.get(h-1)){
					flag[i]=true;
				} else {
					flag[i]=false;
				}
			}
			
		}
		
		for(int i=0;i<Fstore.size();i++){
			for(int h=1;h<data.size();h++){	
				if(Fstore.get(i).timeProb.get(h)>0.2){
					flag[i]=true;
				} else {
					flag[i]=false;
				}
			}
			
		}
		
		for(int i=0;i<Fstore.size();i++){
			if(flag[i]){
				System.out.print(Fstore.get(i).site+"\t");
			}
			for(int h=0;h<data.size();h++){	
				if(flag[i]){
				System.out.print(Fstore.get(i).timeProb.get(h)+"\t");
				}
			}
			if(flag[i]){
			System.out.println();
			}
		}
*/

	}




}
