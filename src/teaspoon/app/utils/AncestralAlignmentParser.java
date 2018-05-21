package teaspoon.app.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;


public class AncestralAlignmentParser {
	File input;

	public AncestralAlignmentParser(String fileName){
			input = new File(fileName);  // The file object
	}

	
//	**********************************************************************
//	Read in sequences from txt file to ArrayList object	
	public int[] read(){
		ArrayList data = new ArrayList();	

		data = readfile();		
		char[] array = make_array(data);						// Convert from array list to matrix
		int[] int_array = convert2int(array);					// Convert to integer matrix
		
		return int_array;
	}

//	**********************************************************************	
//	teaspoon.adaptation.Methods
//	**********************************************************************
	
	public ArrayList readfile()  {
		Scanner s = null;
		ArrayList<String> data = new ArrayList<String>();
		try {
			s = new Scanner(new BufferedReader(new FileReader(input)));

			while (s.hasNext()) {
				data.add(s.next());
			}
		}
		catch (IOException e) {
			System.err.println("Caught IOException: " +  e.getMessage());
		}
		finally {
			if (s != null) {
				s.close();
			}
		}
		return data;
	}
	
	public char[] make_array(ArrayList Readsequence) {
		int start=0;
		// find when tages end and sequences start, by searning for long strings (bases)
		for(int x=0;x<Readsequence.size();x++){
			char[] temp = Readsequence.get(x).toString().toCharArray();	
			if(temp.length>40){
				start = x;
				break;
			}
			
		}
		char[] array = new char[Readsequence.get(start).toString().length()];	
		array = Readsequence.get(start).toString().toCharArray();			

		double iswhole = (array.length)/3;					// Check and see there are no incomplete numCodons
		double part = Math.floor(iswhole);
		if (iswhole/part != 1){
			throw new RuntimeException("incomplete numCodons found. please use whole numCodons");
		}
		
		return array;
	}
	
	public int[] convert2int(char[] array){
		int[] integer_array = new int[array.length];
		for (int i=0; i< array.length; i++ ){
				if (array[i] == 'A'){					// If A, then use integer 1
					integer_array[i] = 1;
				} else if (array[i] == 'C'){			// If T, then use integer 2
					integer_array[i] = 2;
				} else if (array[i] == 'G'){			// If G, then use integer 3
					integer_array[i] = 3;
				} else if (array[i] == 'T'){			// If C, then use integer 4
					integer_array[i] = 4;
				} else if (array[i] == '-'){			// If -, then use integer 0
					integer_array[i] = 5;
				} else {
					integer_array[i] = 5;
				}
		}
		return integer_array;
	}

}



