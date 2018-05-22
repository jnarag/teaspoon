package teaspoon.adaptation;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;


public class DataSet {
	ArrayList<SequenceInfo> DataMainFrame = new ArrayList<SequenceInfo>();
	String DataSetName = new String("");
	private int[][] integer_matrix;
	private String[] taxon_matrix;
	String[][] InfoTable;
	File input;
	File table;
	String WhichGene;
	String tableStringName;
	int tablelength;
	public DataSet(){ // null constructor for methods

	}

	public DataSet(String Name,String HashTableName,String Gene){
		this.DataSetName = Name;
		this.input = new File(Name);  // The file object
		this.table = new File(HashTableName);
		this.WhichGene = Gene;
		this.tableStringName = HashTableName;
		this.tablelength = count(HashTableName);
		readFASTA();
		LoadHashTable(tablelength+2,24);
		CreateDataFrame(WhichGene);
		//	OutputStats();
	}

	public DataSet(String Name){
		this.DataSetName = Name;
		this.input = new File(Name);  // The file object
		readFASTA();
	}

	public void OutputStats(){
		System.out.println("I have read  "+getIntegerMatrix().length+" sequences");
		Print(InfoTable);
	}

	// print matrix
	public void Print(int[][] mat){
		for(int i=0;i<mat.length;i++){
			for(int j=0;j<mat[0].length;j++){
				System.out.print(mat[i][j]);
			}
			System.out.println();
		}
	}
	// overloaded print array
	public void Print(int[] mat){
		for(int i=0;i<mat.length;i++){
			System.out.print(mat[i]);		
		}
		System.out.println();
	}
	public void Print(double[][] mat){
		for(int i=0;i<mat.length;i++){
			for(int j=0;j<mat[0].length;j++){
				System.out.print(mat[i][j]);
			}
			System.out.println();
		}
	}
	public void Print(String[][] mat){
		for(int i=0;i<mat.length;i++){
			for(int j=0;j<mat[0].length;j++){
				System.out.print(mat[i][j]+"\t");
			}
			System.out.println();
		}
	}

	public void CreateDataFrame(){
		for(int i=0;i<getTaxonMatrix().length;i++){
			String line = getTaxonMatrix()[i];
			int[] seq= getIntegerMatrix()[i];
			SequenceInfo temp = new SequenceInfo();
			temp.setTaxon(line);
			temp.setSequence(seq);
			DataMainFrame.add(temp);
		}
	}


	public void CreateDataFrame(String Gene){
		int index = 0;
		for(int i=0;i<getTaxonMatrix().length;i++){
			String line = getTaxonMatrix()[i];
			int[] seq= getIntegerMatrix()[i];
			if(Gene.equals("PB2")){	// matches gene
				SequenceInfo temp = new SequenceInfo();
				for(int x=1;x<InfoTable.length;x++){
					if(line.matches("(?i).*"+InfoTable[x][10]+".*")){  // matches the accession number
						index=x;	
						temp.Accession=InfoTable[index][10];
						temp.GenomeAccession=InfoTable[index][0];
						temp.Name=InfoTable[index][1];
						temp.Genotype=InfoTable[index][2];
						temp.Species=InfoTable[index][18];
						temp.Strain=InfoTable[index][19];
						temp.Country=InfoTable[index][20];
						temp.Continent=InfoTable[index][21];
						temp.Year=Double.valueOf(InfoTable[index][22]);
						temp.DecimalDate=Double.valueOf(InfoTable[index][23]);
						temp.setSequence(seq);
						temp.setTaxon(line);
						temp.Gap = gapInfo(seq);
						DataMainFrame.add(temp);
						break;
					}
				}
			}
			if(Gene.equals("PB1")){
				SequenceInfo temp = new SequenceInfo();
				for(int x=1;x<InfoTable.length;x++){
					if(line.matches("(?i).*"+InfoTable[x][10+1]+".*")){
						index=x;
						temp.Accession=InfoTable[index][10+1];
						temp.GenomeAccession=InfoTable[index][0];
						temp.Name=InfoTable[index][1];
						temp.Genotype=InfoTable[index][2+1];
						temp.Species=InfoTable[index][18];
						temp.Strain=InfoTable[index][19];
						temp.Country=InfoTable[index][20];
						temp.Continent=InfoTable[index][21];
						temp.Year=Double.valueOf(InfoTable[index][22]);
						temp.DecimalDate=Double.valueOf(InfoTable[index][23]);
						temp.setSequence(seq);
						temp.setTaxon(line);
						temp.Gap = gapInfo(seq);
						DataMainFrame.add(temp);
						break;
					}
				}
			}
			if(Gene.equals("PA")){
				SequenceInfo temp = new SequenceInfo();
				for(int x=1;x<InfoTable.length;x++){
					if(line.matches("(?i).*"+InfoTable[x][10+2]+".*")){
						index=x;
						temp.Accession=InfoTable[index][10+2];
						temp.GenomeAccession=InfoTable[index][0];
						temp.Name=InfoTable[index][1];
						temp.Genotype=InfoTable[index][2+2];
						temp.Species=InfoTable[index][18];
						temp.Strain=InfoTable[index][19];
						temp.Country=InfoTable[index][20];
						temp.Continent=InfoTable[index][21];
						temp.Year=Double.valueOf(InfoTable[index][22]);
						temp.DecimalDate=Double.valueOf(InfoTable[index][23]);
						temp.setSequence(seq);
						temp.setTaxon(line);
						temp.Gap = gapInfo(seq);
						DataMainFrame.add(temp);
						break;
					}
				}
			}
			if(Gene.equals("HA") || Gene.equals("H1") || Gene.equals("H3")){
				SequenceInfo temp = new SequenceInfo();
				for(int x=1;x<InfoTable.length;x++){
					if(line.matches("(?i).*"+InfoTable[x][10+3]+".*")){
						index=x;
						temp.Accession=InfoTable[index][10+3];
						temp.GenomeAccession=InfoTable[index][0];
						temp.Name=InfoTable[index][1];
						temp.Genotype=InfoTable[index][2+3];
						temp.Species=InfoTable[index][18];
						temp.Strain=InfoTable[index][19];
						temp.Country=InfoTable[index][20];
						temp.Continent=InfoTable[index][21];
						temp.Year=Double.valueOf(InfoTable[index][22]);
						temp.DecimalDate=Double.valueOf(InfoTable[index][23]);
						temp.setSequence(seq);
						temp.setTaxon(line);
						temp.Gap = gapInfo(seq);
						DataMainFrame.add(temp);
						break;
					}
				}
			}
			if(Gene.equals("NP")){
				SequenceInfo temp = new SequenceInfo();
				for(int x=1;x<InfoTable.length;x++){
					if(line.matches("(?i).*"+InfoTable[x][10+4]+".*")){
						index=x;
						temp.Accession=InfoTable[index][10+4];
						temp.GenomeAccession=InfoTable[index][0];
						temp.Name=InfoTable[index][1];
						temp.Genotype=InfoTable[index][2+4];
						temp.Species=InfoTable[index][18];
						temp.Strain=InfoTable[index][19];
						temp.Country=InfoTable[index][20];
						temp.Continent=InfoTable[index][21];
						temp.Year=Double.valueOf(InfoTable[index][22]);
						temp.DecimalDate=Double.valueOf(InfoTable[index][23]);
						temp.setSequence(seq);
						temp.setTaxon(line);
						temp.Gap = gapInfo(seq);
						DataMainFrame.add(temp);
						break;
					}
				}
			}
			if(Gene.equals("NA") || Gene.equals("N1") || Gene.equals("N2")){
				SequenceInfo temp = new SequenceInfo();
				for(int x=1;x<InfoTable.length;x++){
					if(line.matches("(?i).*"+InfoTable[x][10+5]+".*")){
						index=x;
						temp.Accession=InfoTable[index][10+5];
						temp.GenomeAccession=InfoTable[index][0];
						temp.Name=InfoTable[index][1];
						temp.Genotype=InfoTable[index][2+5];
						temp.Species=InfoTable[index][18];
						temp.Strain=InfoTable[index][19];
						temp.Country=InfoTable[index][20];
						temp.Continent=InfoTable[index][21];
						temp.Year=Double.valueOf(InfoTable[index][22]);
						temp.DecimalDate=Double.valueOf(InfoTable[index][23]);
						temp.setSequence(seq);
						temp.setTaxon(line);
						temp.Gap = gapInfo(seq);
						DataMainFrame.add(temp);
						break;
					}
				}
			}
			if(Gene.equals("MP") || Gene.equals("M1") || Gene.equals("M2")){
				SequenceInfo temp = new SequenceInfo();
				for(int x=1;x<InfoTable.length;x++){
					if(line.matches("(?i).*"+InfoTable[x][10+6]+".*")){
						index=x;
						temp.Accession=InfoTable[index][10+6];
						temp.GenomeAccession=InfoTable[index][0];
						temp.Name=InfoTable[index][1];
						temp.Genotype=InfoTable[index][2+6];
						temp.Species=InfoTable[index][18];
						temp.Strain=InfoTable[index][19];
						temp.Country=InfoTable[index][20];
						temp.Continent=InfoTable[index][21];
						temp.Year=Double.valueOf(InfoTable[index][22]);
						temp.DecimalDate=Double.valueOf(InfoTable[index][23]);
						temp.setSequence(seq);
						temp.setTaxon(line);
						temp.Gap = gapInfo(seq);
						DataMainFrame.add(temp);
						break;
					}
				}
			}
			if(Gene.equals("NS") || Gene.equals("NS1") || Gene.equals("NS2")){
				SequenceInfo temp = new SequenceInfo();
				for(int x=1;x<InfoTable.length;x++){
					if(line.matches("(?i).*"+InfoTable[x][10+7]+".*")){
						index=x;
						temp.Accession=InfoTable[index][10+7];
						temp.GenomeAccession=InfoTable[index][0];
						temp.Name=InfoTable[index][1];
						temp.Genotype=InfoTable[index][2+7];
						temp.Species=InfoTable[index][18];
						temp.Strain=InfoTable[index][19];
						temp.Country=InfoTable[index][20];
						temp.Continent=InfoTable[index][21];
						temp.Year=Double.valueOf(InfoTable[index][22]);
						temp.DecimalDate=Double.valueOf(InfoTable[index][23]);
						temp.setSequence(seq);
						temp.setTaxon(line);
						temp.Gap = gapInfo(seq);
						DataMainFrame.add(temp);
						break;
					}
				}
			}
		}
	}


	public void LoadHashTable(int r,int c) {

		String[][] numbers = new String[r][c];

		try {

			BufferedReader bufRdr = new BufferedReader(new FileReader(table));
			String line = null;
			int row = 0;
			int col = 0;

			//read each line of text file
			while((line = bufRdr.readLine()) != null && row < r)
			{
				StringTokenizer st = new StringTokenizer(line,",");
				while (st.hasMoreTokens())
				{
					//get next token and store it in the array
					numbers[row][col] = st.nextToken();
					col++;
				}
				col = 0;
				row++;
			}

		}	catch (IOException e) {
			System.err.println("Caught IOException: " +  e.getMessage());
		}
		InfoTable= numbers;

	}

	public void exportFASTA(String OutputName, ArrayList<SequenceInfo> ExportData){
		try{
			// Create file 
			System.out.println("-->Exporting Sequence alignment to: "+ OutputName + ": Number of Sequences: " + ExportData.size());
			FileWriter fstream = new FileWriter(OutputName);
			BufferedWriter out = new BufferedWriter(fstream);
			Iterator<SequenceInfo> It =  ExportData.iterator();
			while(It.hasNext()){
				SequenceInfo element = It.next();
				out.write(element.getTaxon());
				out.write("\n");
				out.write(new String(convert2char(element.getSequence())));
				out.write("\n");
			}
			//Close the output stream
			out.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}		
	}

	public void readFASTA() {

		ArrayList<String> data = new ArrayList<String>();
		try {

			BufferedReader br = new BufferedReader(
					new InputStreamReader(
							new FileInputStream(input)));

			String line;
			while((line = br.readLine()) != null)

				if(line.length() != 0){
					data.add(line);
				}

			br.close();		
		}
		catch (IOException e) {
			System.err.println("Caught IOException: " +  e.getMessage());
		}
		data.add(">ForParsing");
		ArrayList<String> sequence = new ArrayList<String>();
		ArrayList<String> taxons = new ArrayList<String>();

		String line = "";
		String Taxon = "";
		String Sequence = "";
		Iterator<String> It =  data.iterator();
		int flag = 0;
		while(It.hasNext()){
			line = It.next();
			if(line.matches("(?i).*>.*")){
				if(flag==1){
					sequence.add(Sequence);			
				}
				flag = 1;
				Taxon=line;
				if(line!=">ForParsing"){
					taxons.add(Taxon);
				}
				Sequence = "";
			} else {
				Sequence += line;
			}
		}

		char[][] matrix = new char[sequence.size()][sequence.get(0).toString().length()];	
		It =  sequence.iterator();	
		int i=0;
		while(It.hasNext()){
			matrix[i] = It.next().toString().toCharArray();
			i++;
		}

		String[] namematrix = new String[taxons.size()];	
		It =  taxons.iterator();	
		i=0;
		while(It.hasNext()){
			namematrix[i] = It.next().toString();
			i++;
		}


		this.setIntegerMatrix(convert2int(matrix));
		this.setTaxonMatrix(namematrix);
	}

	public int[][] readFASTAfast(){
		ArrayList<String> data = new ArrayList<String>();

		try {

			BufferedReader br = new BufferedReader(
					new InputStreamReader(
							new FileInputStream(input)));

			String line;
			while((line = br.readLine()) != null)

				if(line.length() != 0){
					data.add(line);
				}

			br.close();		
		}
		catch (IOException e) {
			System.err.println("Caught IOException: " +  e.getMessage());
		}

		ArrayList<String> sequence = new ArrayList<String>();


		String line = "";
		for(int i=1;i<data.size();i++){
			char[] name = data.get(i).toString().toCharArray();	
			if(new Character(name[0]).equals('>')){
				sequence.add(line);
				line = "";		
			} else {
				line += data.get(i);
			}

		}
		sequence.add(line);

		char[][] matrix = new char[sequence.size()][sequence.get(0).toString().length()];	

		for(int i=0; i< sequence.size();i++){
			matrix[i] = sequence.get(i).toString().toCharArray();			
		}

		int[][] integer_matrix = convert2int(matrix);
		return integer_matrix;
	}

	public int[][] convert2int(char[][] matrix){
		int[][] integer_matrix = new int[matrix.length][matrix[0].length];
		for (int i=0; i< matrix.length; i++ ){
			for(int j=0; j< matrix[0].length; j++){
				if (matrix[i][j] == 'A'){					// If A, then use integer 1
					integer_matrix[i][j] = 1;
				} else if (matrix[i][j] == 'C'){			// If C, then use integer 2
					integer_matrix[i][j] = 2;
				} else if (matrix[i][j] == 'G'){			// If G, then use integer 3
					integer_matrix[i][j] = 3;
				} else if (matrix[i][j] == 'T'){			// If T, then use integer 4
					integer_matrix[i][j] = 4; 
				} else {
					integer_matrix[i][j] = 5;
				}
			}
		}
		return integer_matrix;
	}


	public char[] convert2char(int[] matrix){
		char[] char_matrix = new char[matrix.length];
		for (int i=0; i< matrix.length; i++ ){
			if (matrix[i] == 1){					// If A, then use integer 1
				char_matrix[i] = 'A';
			} else if (matrix[i] == 2){			// If C, then use integer 2
				char_matrix[i] = 'C';
			} else if (matrix[i] == 3){			// If G, then use integer 3
				char_matrix[i] = 'G';
			} else if (matrix[i] == 4){			// If T, then use integer 4
				char_matrix[i] = 'T'; 
			} else {
				char_matrix[i] = '-';
			}

		}
		return char_matrix;
	}

	public double gapInfo(int[] seq){
		double invalidcount = 0.0;
		double count=0;
		for(int j=0;j<seq.length;j++){
			if(seq[j]==5){
				count++;
			}
		}
		invalidcount = 100.0*(count/(double) seq.length);
		return invalidcount;
	}

	public int count(String filename) {
		int count = 0;
		try {
			InputStream is = new BufferedInputStream(new FileInputStream(filename));
			byte[] c = new byte[1024];

			int readChars = 0;
			while ((readChars = is.read(c)) != -1) {
				for (int i = 0; i < readChars; ++i) {
					if (c[i] == '\n' || c[i] == '\r')
						++count;
				}
			}
			is.close();
		}catch (IOException e) {
			System.err.println("Caught IOException: " +  e.getMessage());
		}
		return count+1;
	}

	public int[][] dateCounter(int start, int end){
		int[] dates = new int[end-start];
		for(int i=0;i<dates.length;i++){
			dates[i]=i+start;
		}
		int[] counter = new int[dates.length];
		
		for(int i=0;i<getTaxonMatrix().length;i++){
			String line = getTaxonMatrix()[i];
			for(int j=0;j<dates.length;j++){
				if(line.matches("(?i).*_"+String.valueOf(dates[j])+".*")){
					counter[j]++;
				}			
				
			}
		}
		int[][] vec = new int[counter.length][2];
		for(int i=0;i<counter.length;i++){
			vec[i][1]=counter[i];
			vec[i][0]=dates[i];
			
		}
		return vec;
		
	}

	/**
	 * @return the integerMatrix
	 */
	public int[][] getIntegerMatrix() {
		return integer_matrix;
	}

	/**
	 * @param integerMatrix the integerMatrix to set
	 */
	public void setIntegerMatrix(int[][] integer_matrix) {
		this.integer_matrix = integer_matrix;
	}

	/**
	 * @return the taxon_matrix
	 */
	public String[] getTaxonMatrix() {
		return taxon_matrix;
	}

	/**
	 * @param taxon_matrix the taxon_matrix to set
	 */
	public void setTaxonMatrix(String[] taxon_matrix) {
		this.taxon_matrix = taxon_matrix;
	}


}
