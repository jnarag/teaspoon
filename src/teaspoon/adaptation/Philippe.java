package teaspoon.adaptation;

public class Philippe {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String [] list = {"A_2008-08-04.fasta","A_2008-03-03.fasta","A_2006-03-01.fasta","A_2005-10-17.fasta","A_2000-09-25.fasta"}; // patient A
//		String [] list = {"B_1990-05-07.fasta","C_1994-01-10.fasta","D_1999-05-26.fasta","E_2001-02-21.fasta","F_2002-04-17.fasta","G_2002-04-16.fasta","H_2000-05-04.fasta","I_1999-10-07.fasta","K_2004-09-30.fasta","L_2006-03-24.fasta"};
	//	double nr[] = {0.430795691,0.02662666,0.311831706,0.295596847,0.288161474,5.102194927,0.886829134,0.659794598,0.016315855,0.683172985}; //bhatt
	//	double nr[] ={0.380206544,0.380355224,0.524164209,0.386927997,0.629391129,0.35875859,0.522249225,0.616629582,0.373219096,0.445210172}; //williamson
	double nr[] = {0.449358329};
		for(int run=0;run<list.length;run++){
			// TODO Auto-generated method stub
			DataSet a = new DataSet("/Users/sam/Desktop/HIV/AncestralSeq.fasta");
			DataSet d = new DataSet("/Users/sam/Desktop/HIV/"+list[run]);



		//	final double[] low = {0.0,0.15};
		//	final double[] mid = {0.15,0.75};
		//	final double[] high = {0.75,1.0};	
			final double[] low = {0.0,0.0};
			final double[] mid = {0.0,0.5};
			final double[] high = {0.5,1.0};	

			

			int[] ans = a.integer_matrix[0];
			boolean[] bad = new boolean[d.integer_matrix[0].length];
			Williamson3bin ww = new Williamson3bin(d.integer_matrix,ans,bad);

			ww.badsites_IncludeMainGaps();
			ww.williamson3bin_method_IncludeMainGaps(nr[0], low, mid, high);
			System.out.println(list[run]+"\t"+ww.NumSample+"\t"+ww.low_R + "\t" +ww.low_S+"\t"+ww.mid_R + "\t" +ww.mid_S +"\t"+ww.high_R + "\t" +ww.high_S +"\t"+ww.Adapt);
		}
	}

}
