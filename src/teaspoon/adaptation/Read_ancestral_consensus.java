package teaspoon.adaptation;

import java.io.File;


public class Read_ancestral_consensus {
	File input;
	Read_main Reader; 
	public Read_ancestral_consensus(String Name){
		input = new File(Name);  // The file object
		Reader = new Read_main(Name);
	}

	

//	**********************************************************************
//	Read in sequences from txt file to ArrayList object	
	public int[] read(){
		int[][] int_matrix = Reader.read();	
		int[] int_array = consensusArray(int_matrix);
		return int_array;
	}
	
	public int[] readFASTA(){
		int[][] int_matrix = Reader.readFASTA();				// Convert to integer matrix
		int[] int_array = consensusArray(int_matrix);
		return int_array;
	}
	
	public int[] readNEXUS(){
		int[][] int_matrix = Reader.readNEXUS();				// Convert to integer matrix
		int[] int_array = consensusArray(int_matrix);
		return int_array;
	}

//	**********************************************************************	
//	teaspoon.adaptation.Methods
//	**********************************************************************
	

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

}
