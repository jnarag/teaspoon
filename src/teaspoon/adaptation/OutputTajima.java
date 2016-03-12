package teaspoon.adaptation;

import java.io.File;
import java.io.FileNotFoundException;

import java.util.Formatter;

public class OutputTajima {

	public File file;
	public File[] data;
	public String name;
	
//	Constructers		
	public OutputTajima(){
		// default no-arg constructor
		throw new RuntimeException("please imput seqeunce file location");
	}
//	Constructor	
//	public teaspoon.adaptation.OutputTajima(String inputsequenceslocation){
//	input = new File(inputsequenceslocation);												// Input 
//	data = input.listFiles();
//	file = new File("files.txt");
//	}

	public OutputTajima(String inputsequenceslocation){
		name = inputsequenceslocation;
		file = new File("files.txt");
	}

	public void output(){
		// Formatter Steps
		Formatter formatter = null;
		// Send formatted output to the file.
		try {formatter = new Formatter (file);}
		catch  (FileNotFoundException e) {
			// File not found exception thrown since this is a new file name. However, Formatter will create the new file.
		}

		Read_main readseq = new Read_main(name);
		int[][] m = readseq.read();
		DiversityStats TD = new DiversityStats(m);
		String[] headers = TD.headers();
		double[] values = TD.TajimasDValues();
		formatter.format("%2s",name);
		formatter.format("%n");
		for(int index=0;index<headers.length;index++){
			formatter.format("%-20s", headers[index]);
			formatter.format("%-13s", values[index]);
			formatter.format("%n");
		}
		formatter.format("%n");

		formatter.format("%n");
		formatter.format("%n");

		formatter.flush ();
		formatter.close ();

	}


	/*	public void outputmany(){
		// Formatter Steps
		Formatter formatter = null;
		// Send formatted output to the file.
		try {formatter = new Formatter (file);}
		catch  (FileNotFoundException e) {
			// File not found exception thrown since this is a new file name. However, Formatter will create the new file.
		}

		for(int i=1;i<data.length;i++){
			teaspoon.adaptation.Read_main readseq = new teaspoon.adaptation.Read_main(data[i].getAbsoluteFile().toString());
			int[][] m = readseq.read();
			teaspoon.adaptation.DiversityStats TD = new teaspoon.adaptation.DiversityStats(m);
			String[] headers = TD.headers();
			double[] values = TD.TajimasDValues();
			formatter.format("%2s",data[i]);
			formatter.format("%n");
			for(int index=0;index<headers.length;index++){
				formatter.format("%-20s", headers[index]);
				formatter.format("%-13s", values[index]);
				formatter.format("%n");
			}
			formatter.format("%n");

			formatter.format("%n");
			formatter.format("%n");

		}
		formatter.flush ();
		formatter.close ();

	} */


}
