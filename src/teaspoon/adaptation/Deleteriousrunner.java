package teaspoon.adaptation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;

public class Deleteriousrunner {

	public static void main(String[] args) {

		String sss = new String("/Users/sam/Desktop/96VD/");
		File input = new File(sss+"names");		// Input 
		File[] list = input.listFiles();
		ArrayList<String> names = new ArrayList<String>();
		for(int i=0;i<list.length;i++){
			if(list[i].isHidden() == false && list[i].isDirectory() == false)	{
				names.add(list[i].getName());
			}
		}	
		Collections.sort(names);
//		************************************************** Outside LOOP	*************************************************	
		for(int x=0;x<names.size();x++){
			double[][] matrix = new double[300][14]; //Size of Each XLS file

//			************************************************** Inside LOOP	*************************************************	 
			for(int z=0;z<300;z++){

				Read_main ra = new Read_main(sss+"tseq/"+String.valueOf(z+1)+"."+names.get(x).toString());
				Read_main re = new Read_main(sss+"tans/"+String.valueOf(z+1)+"."+names.get(x).toString());
				int[][] tmp = re.read();
				int[] ans = re.consensusArray(tmp);	
				int[][] seq = ra.read();
				Methods m = new Methods();
				double[] L = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
				double[] H = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7};	
				boolean[] Nvec = {false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true};

				double[][] bins = new double[2][L.length];
				for(int i=0;i<L.length;i++){
					bins[0][i]=L[i];
					bins[1][i]=H[i];
				}
				double[] prior = {1,1,1,1};
				SiteEstMulti sm = new SiteEstMulti(seq,ans);
				double[][] val = sm.Value_Matrix(Nvec, bins,prior,true);	

				for(int i=0;i<val[0].length;i++){
					matrix[z][i] = val[3][i];
				}

			}
//			**********************************************************************************************************************		
			File file = new File("/Users/sam/Desktop/Sto/"+names.get(x).toString()+".xls");
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

			System.out.println("************* "+x+" ****************");
		}
	}

}
