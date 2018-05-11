package teaspoon.adaptation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;

import jebl.evolution.taxa.Taxon;


public class Read_main {
    public File input;
    public List<Double> sampleTimes;
    public List<String> sampleNames;
    public List<Taxon> taxa;
    public boolean getSampleTimes;
    int[][] sequenceMatrix;

    public Read_main(String FileName){
        input = new File(FileName);  // The file object
        getSampleTimes = false;


    }

    public Read_main(String FileName, boolean getSampleTimes){
        input = new File(FileName);  // The file object
        sampleTimes = new ArrayList<Double>();
        sampleNames = new ArrayList<String>();
        this.getSampleTimes = getSampleTimes;


    }

    public Read_main(){
        // The file object
    }

    //	**********************************************************************
//	Read in sequences from txt file to ArrayList object	 general read file
    public int[][] read(){
        ArrayList data = new ArrayList();

        data = readfile();
        char[][] matrix = make_matrixNEXUS(data);						// Convert from array list to matrix
        int[][] int_matrix = convert2int(matrix);					// Convert to integer matrix

        return int_matrix;
    }

    //	**********************************************************************
//	teaspoon.adaptation.Methods
//	**********************************************************************
    public int[][] readNEXUS(){
        ArrayList readSequence = new ArrayList();
        readSequence = readfile();
        ArrayList<String> store = new ArrayList<String>();
        // find when tages end and sequences start, by searching for long strings (bases)
        for(int x=0;x<readSequence.size();x++){
            char[] temp = readSequence.get(x).toString().toCharArray();
            if(temp.length>100){
                store.add(readSequence.get(x).toString());
            }
        }
        char[][] matrix = new char[store.size()][store.get(0).toString().length()];
        for(int i=0;i<store.size();i++){
            matrix[i] = store.get(i).toString().toCharArray();
        }
        int[][] integer_matrix = convert2int(matrix);
        return integer_matrix;

    }

    public int[][] readFASTA(){

        ArrayList<String> data = new ArrayList<String>();

        try {

            BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(input)));
            String line;
            while((line = br.readLine()) != null)

                if(line.length() != 0){
                    data.add(line);
                }

            br.close();
        }
        catch (IOException e) {
            System.err.println("Caught IOException: " +  e.getMessage());
            System.out.println(input);
        }

        ArrayList<String> sequence = new ArrayList<String>();

        String line = "";
        for(int i=0;i< data.size();i++){
            line = data.get(i);
            char[] name = line.toString().toCharArray();


            if(!new Character(name[0]).equals('>')){

                sequence.add(line);
            }
            else{
                //System.out.println(line);
                if(getSampleTimes) {
                    String [] parts = line.split("_");
                    double date = Double.parseDouble(parts[parts.length - 1]);
                    sampleTimes.add(date);
                    String seqname = line.replace(">", "");
                    sampleNames.add(seqname);
                    Taxon t = Taxon.getTaxon(seqname);
                    t.setAttribute("time",date);

                }
            }

        }
        //sequence.add(line);
        //System.out.println(sequence.get(0).toString().length());
        if(sequence.size() == 0) {
            return null;
        }
        char[][] matrix = new char[sequence.size()][sequence.get(0).toString().length()];


        //I could limit this so the sequences are read from non-zero start to elsewhere// helpful for whole genome mainAnalysis
        for(int i=0; i< sequence.size();i++){
            matrix[i] = sequence.get(i).toString().toCharArray();
        }

        sequenceMatrix = convert2int(matrix);

        return sequenceMatrix;
    }

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

    public char[][] make_matrix(ArrayList Readsequence) {
        int start=0;
        // find when tages end and sequences start, by searching for long strings (bases)
        for(int x=0;x<Readsequence.size();x++){
            char[] temp = Readsequence.get(x).toString().toCharArray();
            if(temp.length>30){
                start = x;
                break;
            }

        }
        int nrow = ((Readsequence.size() - start)+1)/2;

        char[][] matrix = new char[nrow][Readsequence.get(start).toString().length()];

        for(int i=start, k = 0; i< Readsequence.size();k++, i=i+2){
            matrix[k] = Readsequence.get(i).toString().toCharArray();
        }
        double iswhole = (double)(matrix[0].length)/3;					// Check and see there are no incomplete codons
        double part = Math.floor(iswhole);
        if (iswhole/part != 1){
            throw new RuntimeException("incomplete codons found. please use whole codons");
        }

        return matrix;
    }

    public char[][] make_matrixNEXUS(ArrayList Readsequence) {
        ArrayList<String> store = new ArrayList<String>();
        // find when tages end and sequences start, by searching for long strings (bases)
        for(int x=0;x<Readsequence.size();x++){
            char[] temp = Readsequence.get(x).toString().toCharArray();
            if(temp.length>100){
                store.add(Readsequence.get(x).toString());
            }
        }
        char[][] matrix = new char[store.size()][store.get(0).toString().length()];
        for(int i=0;i<store.size();i++){
            matrix[i] = store.get(i).toString().toCharArray();
        }
        return matrix;

    }



    public int[][] convert2int(char[][] matrix){
        //System.out.println(matrix[0].length+","+matrix.length);
        int[][] integer_matrix = new int[matrix.length][matrix[0].length];
        for (int i=0; i< matrix.length; i++ ){

            if(matrix[i].length != matrix[0].length) {
//
//                System.out.println(matrix[i].length);

                break;
            }
            //System.out.println(sampleNames.get(i));
            for(int j=0; j< matrix[0].length; j++){

                //System.out.println(matrix[0].length+","+matrix.length);
                if (matrix[i][j] == 'A'){					// If A, then use integer 1
                    integer_matrix[i][j] = 1;
                } else if (matrix[i][j] == 'C'){			// If C, then use integer 2
                    integer_matrix[i][j] = 2;
                } else if (matrix[i][j] == 'G'){			// If G, then use integer 3
                    integer_matrix[i][j] = 3;
                } else if (matrix[i][j] == 'T'){			// If T, then use integer 4
                    integer_matrix[i][j] = 4;
                } else {
                    integer_matrix[i][j] = 5;               // If gap or N, use integer 5
                }
            }
        }
        return integer_matrix;
    }

    public int[][] subMatrix(int start, int end) {


        int window_length = end-start;
        List<int[]> sub_matrix = new ArrayList<int[]>();

        int limit = (int) Math.round(window_length*0.05);     // tolerates 5% of total sites with gaps

        int n = 0;
        List<String> cnames = new ArrayList<>();
        for (int[] aSequenceMatrix : sequenceMatrix) {



            int[] sub_sub_mat = Arrays.copyOfRange(aSequenceMatrix, start, end);

            if (gapCheck(limit, sub_sub_mat)) {
                sub_matrix.add(sub_sub_mat);
                //cnames.add(sampleNames.get(n));
            }


        }

        int[][] sub_arrayOfarray = new int[sub_matrix.size()][window_length];

        if(sub_matrix.size()>100000) {

            List<Integer> indices = new ArrayList<Integer>();
            sub_arrayOfarray = new int[100000][window_length];
            List<Double> new_sub_sampleTimes = new ArrayList<Double>();
            for(int i = 0; i< sub_matrix.size(); i++) {

                indices.add(i);
            }

            Collections.shuffle(indices);

            for(int i = 0; i < 100000; i++) {

                sub_arrayOfarray[i] = sub_matrix.get(indices.get(i));
                //new_sub_sampleTimes.add(sub_sampleTimes.get(indices.get(i)));

            }
            //sub_sampleTimes.clear();
            //sub_sampleTimes.addAll(new_sub_sampleTimes);
        }
        else{

            sub_arrayOfarray = new int[sub_matrix.size()][window_length];


            for(int i=0; i < sub_matrix.size(); i++) {

                sub_arrayOfarray[i] = sub_matrix.get(i);
            }
        }


        return sub_arrayOfarray;
    }

    public int[][] subMatrix(int[][] sequenceMatrix, int start, int end) {


        int window_length = end-start;
        List<int[]> sub_matrix = new ArrayList<int[]>();

        int limit = (int) Math.round(window_length*0.05);     // tolerates 5% of total sites with gaps

        int n = 0;
        List<String> cnames = new ArrayList<>();
        for (int[] aSequenceMatrix : sequenceMatrix) {



            int[] sub_sub_mat = Arrays.copyOfRange(aSequenceMatrix, start, end);

            if (gapCheck(limit, sub_sub_mat)) {
                sub_matrix.add(sub_sub_mat);
                cnames.add(sampleNames.get(n));
            }

            n++;
        }

        int[][] sub_arrayOfarray = new int[sub_matrix.size()][window_length];

        for(int i=0; i < sub_matrix.size(); i++) {

            sub_arrayOfarray[i] = sub_matrix.get(i);
        }


        return sub_arrayOfarray;
    }

    public int[][] subMatrix(int[][] sequenceMatrix, List<Double> sampleTimes, int start, int end) {

        // do i need to create an object that comprises submatrix and list of sampletimes?

        int window_length = end-start;
        List<int[]> sub_matrix = new ArrayList<int[]>();
        List<Double> sub_sampleTimes = new ArrayList<Double>();
        int limit = (int) Math.round(window_length*0.05);     // tolerates 5% of total sites with gaps

        for (int i = 0; i < sequenceMatrix.length; i++) {

            int[] sub_sub_mat = Arrays.copyOfRange(sequenceMatrix[i], start, end);

            if (gapCheck(limit, sub_sub_mat)) {
                sub_matrix.add(sub_sub_mat);
                sub_sampleTimes.add(sampleTimes.get(i));
            }


        }
        int [][] sub_array_of_array;
        if(sub_matrix.size()>3000) {

            List<Integer> indices = new ArrayList<Integer>();
            sub_array_of_array = new int[3000][window_length];
            List<Double> new_sub_sampleTimes = new ArrayList<Double>();
            for(int i = 0; i< sub_matrix.size(); i++) {

                indices.add(i);
            }

            Collections.shuffle(indices);

            for(int i = 0; i < 3000; i++) {

                sub_array_of_array[i] = sub_matrix.get(indices.get(i));
                new_sub_sampleTimes.add(sub_sampleTimes.get(indices.get(i)));

            }
            sub_sampleTimes.clear();
            sub_sampleTimes.addAll(new_sub_sampleTimes);
        }
        else{

            sub_array_of_array = new int[sub_matrix.size()][window_length];


            for(int i=0; i < sub_matrix.size(); i++) {

                sub_array_of_array[i] = sub_matrix.get(i);
            }
        }

        sampleTimes.removeAll(sampleTimes);
        sampleTimes.addAll(sub_sampleTimes);

        return sub_array_of_array;

    }

    public boolean gapCheck(int limit, int[] codonArray) {

        boolean gapLessThanLimit = true;
        int count = 0;
        for (int aCodonArray : codonArray) {

            if (aCodonArray == 5) { // 5 is 'gap'
                count++;
            }

            //only one gap allowed
            if (count > limit) {
                gapLessThanLimit = false;
                break;
            }

        }
        return gapLessThanLimit;
    }

    public int[] consensusArray(int[][] integer_matrix){
        Methods preprocess = new Methods();
        int[] consensus = new int[integer_matrix[0].length];
        //ignore gaps - so counter will be array of 4, not 5
        double[] counter = new double[4];
        for(int site=0;site<integer_matrix[0].length;site++){
            // count numbers
            counter[0] = preprocess.num_of_base(integer_matrix, 1, site);
            counter[1] = preprocess.num_of_base(integer_matrix, 2, site);
            counter[2] = preprocess.num_of_base(integer_matrix, 3, site);
            counter[3] = preprocess.num_of_base(integer_matrix, 4, site);
            //counter[4] = preprocess.num_of_base(integer_matrix, 5, site);
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


