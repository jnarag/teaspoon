package teaspoon.SUPGMA;

import jebl.evolution.alignments.BasicAlignment;
import jebl.evolution.distances.DistanceMatrix;
import jebl.evolution.distances.HKYDistanceMatrix;
import jebl.evolution.distances.JukesCantorDistanceMatrix;
import jebl.evolution.distances.TamuraNeiDistanceMatrix;
import jebl.evolution.graphs.Node;
import jebl.evolution.io.FastaImporter;
import jebl.evolution.io.ImportException;
import jebl.evolution.sequences.BasicSequence;
import jebl.evolution.sequences.Sequence;
import jebl.evolution.sequences.SequenceType;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;
import jebl.util.ProgressListener;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import java.io.*;
import java.util.*;


/**
 * Created by jayna on 14/02/2014.
 */
public class SUPGMA {


    List<Double> samplingTimes;
    List<String> sampleNames;
    String [] genes;
    List<Double> samplingTimesSet;
    List<Taxon> taxa;
    List<Sequence> sequences;
    int date_position_in_name;

    String sequenceFileName;
    public SUPGMA(String filename, int date_position_in_name) {

        this.sequenceFileName = filename;
        samplingTimes = new ArrayList<Double>();
//        samplingTimesSet = null;
//        sampleNames = null;
//        taxa = null;
        genes = new String [] {"gag","pol","env","nef"};
        this.date_position_in_name = date_position_in_name;


    }

    public DistanceMatrix getDistances(char model) {
        //need to convert distance matrix into single array

        //associate appropriate sampling interval

        //decide whether to estimate theta per sample timepoint or constant?

        DistanceMatrix pairwiseDistances = null;

        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(sequenceFileName));
            FastaImporter importer = new FastaImporter(reader, SequenceType.NUCLEOTIDE);
            this.sequences = importer.importSequences();
            BasicAlignment alignment = new BasicAlignment(sequences);
            reader.close();
            switch(model) {
                case 'J':
                    pairwiseDistances = new JukesCantorDistanceMatrix(alignment, ProgressListener.EMPTY );
                    break;
                case 'H':
                    pairwiseDistances = new HKYDistanceMatrix(alignment, ProgressListener.EMPTY);
                    break;
                case 'T':
                    pairwiseDistances = new TamuraNeiDistanceMatrix(alignment, ProgressListener.EMPTY);

            }

            assert pairwiseDistances != null;
            List<Taxon> taxa = pairwiseDistances.getTaxa();

            samplingTimes.clear();
            for(Taxon t: taxa) {

                String name = t.getName();
                String [] parts =  name.split("_");
                Double date = Double.parseDouble(parts[date_position_in_name-1].trim());

                samplingTimes.add(date);
                //System.out.println(date);


            }



        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return pairwiseDistances;
    }

    public DistanceMatrix getDistances(char model, List<Sequence> sequences) {
        //need to convert distance matrix into single array

        //associate appropriate sampling interval

        //decide whether to estimate theta per sample timepoint or constant?

        DistanceMatrix pairwiseDistances = null;

        BasicAlignment alignment = new BasicAlignment(sequences);
        switch(model) {
            case 'J':
                pairwiseDistances = new JukesCantorDistanceMatrix(alignment, ProgressListener.EMPTY );
                break;
            case 'H':
                pairwiseDistances = new HKYDistanceMatrix(alignment, ProgressListener.EMPTY);
                break;
            case 'T':
                pairwiseDistances = new TamuraNeiDistanceMatrix(alignment, ProgressListener.EMPTY);

        }

        assert pairwiseDistances != null;
        List<Taxon> taxa = pairwiseDistances.getTaxa();

        samplingTimes.clear();
        for(Taxon t: taxa) {

            String name = t.getName();
            String [] parts =  name.split("_");
            Double date = Double.parseDouble(parts[date_position_in_name-1].trim());

            samplingTimes.add(date);
            //System.out.println(date);


        }



        return pairwiseDistances;
    }

    private List<Double> getSamplingTimesSet(List<Double> samplingTimes) {



        List<Double> sampleTimesSet_list = new ArrayList<Double>();
        Set<Double> sampleTimesSet = new HashSet<Double>();
        sampleTimesSet.addAll(samplingTimes);
        sampleTimesSet_list.addAll(sampleTimesSet);
        Collections.sort(sampleTimesSet_list);

        return sampleTimesSet_list;
    }

    private double[] estimateMultiTheta(DistanceMatrix distanceMatrix) {

        double[][] pairwiseDistances = distanceMatrix.getDistances();


        int dimensions = (int)Math.floor((Math.pow(pairwiseDistances.length, 2)-pairwiseDistances.length)*0.5);

        int c = 0;
        int p = 0;
        boolean multipleTheta = false;

        Set<Double> uniqueSamplingTimes = new HashSet<Double>();
        uniqueSamplingTimes.addAll(samplingTimes);
        List<Double> timepoints = new ArrayList<Double>();
        timepoints.addAll(uniqueSamplingTimes);
        Collections.sort(timepoints);

        int no_of_timepoints = uniqueSamplingTimes.size();

        double [] y = new double[dimensions];
        double [][] linearModelMatrix = new double[dimensions][no_of_timepoints+1];

        for(int i=0; i < pairwiseDistances.length-1; i++) {

            for(int j=(i+1); j < pairwiseDistances.length; j++) {

                y[c] = pairwiseDistances[i][j];

                if(samplingTimes.get(i)-samplingTimes.get(j)==0) {

                    linearModelMatrix[c][timepoints.indexOf(samplingTimes.get(i))] = 1;
                }

                linearModelMatrix[c][no_of_timepoints] = Math.abs(samplingTimes.get(i)-samplingTimes.get(j));

                c++;



            }
        }

        //System.out.println(linearModelMatrix.length + " c "+ c);

        OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
        regression.setNoIntercept(true);
        regression.newSampleData(y, linearModelMatrix);


        double[] beta = regression.estimateRegressionParameters();

        for(double d: beta) {
            System.out.println(d);
        }
        return beta;
    }

    private double[] estimateMultiTheta(DistanceMatrix distanceMatrix, List<Double> samplingTimes) {


        double[][] pairwiseDistances = distanceMatrix.getDistances();

        int dimensions = (int)Math.floor((Math.pow(pairwiseDistances.length, 2)-pairwiseDistances.length)*0.5);

        int c = 0;
        int p = 0;


        Set<Double> uniqueSamplingTimes = new HashSet<Double>();
        uniqueSamplingTimes.addAll(samplingTimes);
        samplingTimesSet = new ArrayList<Double>();
        samplingTimesSet.addAll(uniqueSamplingTimes);
        Collections.sort(samplingTimesSet);

        int no_of_timepoints = uniqueSamplingTimes.size();

        double [] y = new double[dimensions];
        double [][] linearModelMatrix = new double[dimensions][no_of_timepoints+1];

        //System.out.println("dimensions: "+dimensions);
        for(int i=0; i < pairwiseDistances.length-1; i++) {

            for(int j=(i+1); j < pairwiseDistances.length; j++) {

                y[c] = pairwiseDistances[i][j];

                if(samplingTimes.get(i)-samplingTimes.get(j)==0) {

                    linearModelMatrix[c][samplingTimesSet.indexOf(samplingTimes.get(i))] = 1;
                }
                //System.out.println(i+","+j+" :"+samplingTimes.get(i)+" and "+samplingTimes.get(j)+ " pop "+linearModelMatrix[c][timepoints.indexOf(samplingTimes.get(i))]);

                linearModelMatrix[c][no_of_timepoints] = Math.abs(samplingTimes.get(i)-samplingTimes.get(j));

                c++;



            }
        }
        OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
        regression.setNoIntercept(true);
        regression.newSampleData(y, linearModelMatrix);

        try{
            double[] beta = regression.estimateRegressionParameters();
            int d = 0;
            while (d < beta.length-1) {
                System.out.println(samplingTimesSet.get(d)+":"+beta[d]);
                d++;
            }
            System.out.println("rate: "+beta[beta.length-1]);


            System.out.println();



            return beta;
        }
        catch(SingularMatrixException e) {

            System.out.println("Exception thrown :"+ e);
            return null;

        }
//        Array2DRowRealMatrix mat = new Array2DRowRealMatrix(regression.estimateRegressionParameters());
//        QRDecomposition e = new QRDecomposition(mat);
//        boolean determinant = e.getSolver().isNonSingular();
//        System.out.println(determinant);


        //System.out.print(y.length+","+linearModelMatrix.length+","+linearModelMatrix[0].length);




    }

    private double[][] updateDistanceMatrix(double[]beta, DistanceMatrix distanceMatrix) {

        double[][] pairwiseDistances = distanceMatrix.getDistances();

        double[][] updatedDistanceMatrix = new double[pairwiseDistances.length][pairwiseDistances.length];

        double mostRecentSamplingTime = Collections.max(samplingTimesSet);
        for(int i=0; i < pairwiseDistances.length; i++) {

            String s = "";

            for (int j = 0; j < pairwiseDistances.length; j++) {


                int index1 = samplingTimesSet.indexOf(samplingTimes.get(i));
                int index2 = samplingTimesSet.indexOf(samplingTimes.get(j));


                int min_index1 = Math.min(index1, index2);
                int max_index2 = Math.max(index1, index2);


                double samplingInterval = (mostRecentSamplingTime-samplingTimes.get(i))+(mostRecentSamplingTime-samplingTimes.get(j))     ;


                updatedDistanceMatrix[i][j] = pairwiseDistances[i][j]+beta[beta.length-1]*samplingInterval;

//                double total_dist = 0.0;
//
//                if((max_index2>0) && i!=j) {
//
//                    double dist_min = 0.0;
//                    double dist_max = 0.0;
//
//
//                    if(min_index1>0) {
//
//
//                        dist_min = beta[min_index1-1];
//
//                        for(int a=0; a<min_index1; a++) {
//
//                            dist_min+= beta[a];
//                            //System.out.println(a);
//                        }
//                    }
//
//                    for(int b=0; b<max_index2; b++) {
//
//                        dist_max+= beta[b];
//                        //System.out.println(b);
//
//                    }
//
//                    total_dist = dist_min+dist_max;
//
//                    updatedDistanceMatrix[i][j] +=  total_dist;
                //System.out.println(i+","+j+","+pairwiseDistances[i][j]+","+updatedDistanceMatrix[i][j]+","+index1+","+index2+","+samplingTimes.get(j));
                //pairwiseDistances[i][j] += total_dist;


                //System.out.println(total_dist);
                //              }


            }
        }

//        for(int i = 0; i < updatedDistanceMatrix.length; i++) {
//
//            String s = "";
//            for(int j=0; j < updatedDistanceMatrix.length-1; j++) {
//
//                s+=","+updatedDistanceMatrix[i][j];
//            }
//            System.out.println(s);
//        }
        return updatedDistanceMatrix;
    }

    private double[] estimateTheta() {

        DistanceMatrix distanceMatrix = getDistances('H');
        double[][] pairwiseDistances = distanceMatrix.getDistances();

        int dimensions = (int)Math.floor((Math.pow(pairwiseDistances.length, 2)-pairwiseDistances.length)*0.5);

        int c = 0;

        double [] y = new double[dimensions];
        double [][] linearModelMatrix = new double[dimensions][2];

        for(int i=0; i < pairwiseDistances.length-1; i++) {

            for(int j=(i+1); j < pairwiseDistances.length; j++) {


                y[c] = pairwiseDistances[i][j];

                linearModelMatrix[c][0] = 1;
                //System.out.println(i+","+j+" :"+samplingTimes.get(i)+" and "+samplingTimes.get(j)+ " pop "+linearModelMatrix[c][0]);

                linearModelMatrix[c][1] = (samplingTimes.get(i) - samplingTimes.get(j));


                c++;
            }


        }

        //System.out.println(linearModelMatrix.length + " c "+ c);

        OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
        regression.setNoIntercept(true);
        regression.newSampleData(y, linearModelMatrix);


        double[] beta = regression.estimateRegressionParameters();

        for(double d: beta) {
            System.out.println(d);
        }
        System.out.println();
        return beta;
    }

    private double[] estimateTheta(DistanceMatrix distanceMatrix, List<Double> samplingTimes) {


        double [][] pairwiseDistances = distanceMatrix.getDistances();
        int dimensions = (int)Math.floor((Math.pow(pairwiseDistances.length, 2)-pairwiseDistances.length)*0.5);

        Set<Double> uniqueSamplingTimes = new HashSet<Double>();
        uniqueSamplingTimes.addAll(samplingTimes);
        samplingTimesSet = new ArrayList<Double>();
        samplingTimesSet.addAll(uniqueSamplingTimes);
        Collections.sort(samplingTimesSet);

        System.out.println("dimensions: "+dimensions);
        int c = 0;

        double [] y = new double[dimensions];
        double [][] linearModelMatrix = new double[dimensions][2];

        //assuming single/constant population
        for(int i=0; i < pairwiseDistances.length-1; i++) {

            for(int j=(i+1); j < pairwiseDistances.length; j++) {

                y[c] = pairwiseDistances[i][j];

                linearModelMatrix[c][0] = 1;
                //System.out.println(i+","+j+" :"+samplingTimes.get(i)+" and "+samplingTimes.get(j)+ " pop "+linearModelMatrix[c][0]);

                linearModelMatrix[c][1] = Math.abs(samplingTimes.get(i) - samplingTimes.get(j));

                c++;


            }
        }


        //System.out.println(linearModelMatrix.length + " c "+ c);

        OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
        regression.setNoIntercept(true);
        regression.newSampleData(y, linearModelMatrix);


        double[] beta = regression.estimateRegressionParameters();

        for(double d: beta) {
            System.out.println(d);
        }
        return beta;
    }

    private List<Sequence> subAlignment(int start, int end) {

        List<Sequence> subAlignment = new ArrayList<>();
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(sequenceFileName));
            FastaImporter importer = new FastaImporter(reader, SequenceType.NUCLEOTIDE);
            this.sequences = importer.importSequences();
            BasicAlignment alignment = new BasicAlignment(sequences);
            reader.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        for(Sequence s: sequences) {

            String seq = s.getString();
            if((start-1) < seq.length() && (end < seq.length())) {
                String subSeq = seq.substring(start - 1, end);
                int count = StringUtils.countMatches(subSeq, '-');
                if (count <= subSeq.length() * 0.05) {
                    BasicSequence newSeq = new BasicSequence(SequenceType.NUCLEOTIDE, s.getTaxon(), subSeq);
                    subAlignment.add(newSeq);
                }
            }
        }

        return subAlignment;
    }

    private void estimateThetaPerWindow(int geneStart, int geneEnd, int window_length) {

        int geneLength = geneEnd - geneStart;
        int no_of_windows = (int)Math.floor(geneLength/window_length);

        int window_start = geneStart;

        List<double[]> results = new ArrayList<>();

        for(int i = 0; i < no_of_windows; i++) {

            int window_end = window_start+window_length;
            if(window_end > geneEnd) {

                window_end = geneEnd;
            }
            List<Sequence> subSequences = subAlignment(window_start,window_end);
            DistanceMatrix distanceMatrix = getDistances('H', subSequences);

            List<Double> subSampleTime = new ArrayList<Double>();
            subSampleTime.addAll(samplingTimes);

            System.out.println("window "+ i+ " at pos: "+window_start+" to "+window_end);

            double[] beta = estimateTheta(distanceMatrix, subSampleTime);

            window_start+= window_length;
            //System.out.println(beta[0]+","+beta[1]);

        }


        //return array of beta...

    }

    private void estimateMultiThetaPerWindow(int geneStart, int geneEnd, int window_length) {

        int geneLength = geneEnd - geneStart;
        int no_of_windows = (int) Math.round(geneLength / (double) window_length);

        System.out.println("windows "+no_of_windows);
        int window_start = geneStart;

        List<double[]> results = new ArrayList<>();
        List<List<Double>> times_list = new ArrayList<>();
        for (int i = 0; i < no_of_windows; i++) {

            System.out.println("w "+(i+1));
            int window_end = window_start + window_length;
            if (window_end > geneEnd) {

                window_end = geneEnd;
            }
            List<Sequence> subSequences = subAlignment(window_start, window_end);
            DistanceMatrix distanceMatrix = getDistances('H', subSequences);

            List<Double> subSampleTime = new ArrayList<Double>();
            subSampleTime.addAll(samplingTimes);


            double[] beta = null;

            if(subSequences.size()>100) {
                beta = estimateMultiTheta(distanceMatrix, subSampleTime);
            }
            System.out.println("window " + i + " at pos: " + window_start + " to " + window_end+","+samplingTimesSet.size());


            times_list.add(samplingTimesSet);

            if (beta == null) {

                beta = new double[samplingTimesSet.size() + 1];
                Arrays.fill(beta, Double.NaN);

            }


            //System.out.println(beta.length+","+samplingTimesSet.size());
            results.add(beta);
//            else{
//
//            }
            window_start += window_length;
            //System.out.println(beta[0]+","+beta[1]);

        }


        writeResults(results, times_list, String.valueOf(geneStart) + "_" + String.valueOf(geneEnd) + "_SUPGMA_results.csv");

//        String r = "";
//        for (int d=0; d<(samplingTimesSet.size()+1); d++) {
//
//
//            if(d<samplingTimesSet.size()) {
//                r += String.valueOf(samplingTimesSet.get(d));
//            }
//            else{
//                r += "rate";
//            }
//            for (double[] beta : results) {
//
//                //if(d<samplingTimesSet.size()) {
//                    r += ","+(String.valueOf(beta[d]));
////                }
////                else{
////                    r += ","+String.valueOf(beta[beta.length-1]);
////                }
//
//            }
//            r+="\n";
//
//
//        }
//        System.out.println(r);



        //return array of beta...

    }

    private void writeResults(List<double[]> results, List<List<Double>> time_list, String filename) {

        filename = sequenceFileName.replace(".fasta","_"+filename);
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(filename)));


            //String r = "Time";


            Map<Double, String> times_string_map = new HashMap<>();

            times_string_map.put(3000.0,"");
            for(int i=0; i<time_list.size(); i++) {

                List<Double> times = time_list.get(i);
                double beta[] = results.get(i);

                Map<Double, Double> times_beta_map = new HashMap<>();


                for (int t = 0; t < times.size(); t++) {

                    times_beta_map.put(times.get(t), beta[t]);
                    if (!times_string_map.containsKey(times.get(t))) {

                        times_string_map.put(times.get(t), "");
                    }

                }

                for (Double t : times_beta_map.keySet()) {

                    if (times_string_map.containsKey(t)) {

                        times_string_map.put(t, times_string_map.get(t) +","+ times_beta_map.get(t));
                    } else {
                        times_string_map.put(t, times_string_map.get(t) +","+Double.NaN);
                    }
                }
                times_string_map.put(3000.0, times_string_map.get(3000.0) + "," + beta[beta.length-1]);
            }


            List<Double> times = new ArrayList<>(times_string_map.keySet());
            Collections.sort(times);
            String s = "Time";

            for(int r=0; r<results.size(); r++){
                s+=",Window"+(r+1);
            }
            s+="\n";

            for(Double t: times) {

                if(t==3000.0) {
                    s+="\nrate"+times_string_map.get(t)+"\n";
                    //System.out.println("rate"+times_string_map.get(t));

                }
                else {
                    s+=t + times_string_map.get(t)+"\n";
                    System.out.println(t + times_string_map.get(t));
                }
            }


            writer.write(s+"\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

//    private void estimateMultiThetaPerWindow(int geneStart, int geneEnd, int window_length) {
//
//        int geneLength = geneEnd - geneStart;
//        int no_of_windows = (int)Math.floor(geneLength/window_length);
//
//        int window_start = geneStart;
//        for(int i = 0; i < no_of_windows; i++) {
//
//            List<Double> subSampleTime = new ArrayList<Double>();
//            subSampleTime.addAll(distanceMatrix.getSampleTimes());
//
//            System.out.println("window "+ i+ " at pos: "+window_start+" to "+(window_length+window_start));
//            int[][] subMatrix = distanceMatrix.main.subMatrix(distanceMatrix.sequenceMatrix, subSampleTime, window_start, window_start+window_length);
//            //I will also need to get a subMatrix of the sampleTImes
//            //System.out.println("seqs: "+subMatrix.length);
//
//            double [][] pairwiseDistances = distanceMatrix.getHKYDistanceMatrix(subMatrix);
//
//            //System.out.println(pairwiseDistances.length+": "+subSampleTime.size());
//            double[] beta = estimateMultiTheta(pairwiseDistances, subSampleTime);
//            System.out.println();
//            window_start+= window_length;
//
//
//        }
//    }
//
//    private void seqs_per_window(int geneStart, int geneEnd, int window_length) {
//
//        int geneLength = geneEnd - geneStart;
//        int no_of_windows = (int)Math.floor(geneLength/window_length);
//
//        int window_start = geneStart;
//        for(int i = 0; i < no_of_windows; i++) {
//
//            List<Double> subSampleTime = new ArrayList<Double>();
//            subSampleTime.addAll(distanceMatrix.getSampleTimes());
//
//            System.out.println("window "+ i+ " at pos: "+window_start+" to "+(window_length+window_start));
//            int[][] subMatrix = distanceMatrix.main.subMatrix(distanceMatrix.sequenceMatrix, subSampleTime, window_start, window_start+window_length);
//            //I will also need to get a subMatrix of the sampleTImes
//            //System.out.println("seqs: "+subMatrix.length);
//
//            System.out.println();
//            window_start+= window_length;
//
//
//        }
//
//    }
//
//    private void seqs_per_gene(int[][] genes, int window_length) {
//
//        int no_genes = genes.length;
//        for(int g=0; g<no_genes; g++) {
//
//            System.out.println(this.genes[g]);
//            int[] gene_boundary = genes[g];
//            int gene_start = gene_boundary[0]-1;
//            int gene_end = gene_boundary[1]-1;
//            int gene_length = gene_end-gene_start;
//
//            seqs_per_window(gene_start, gene_end, window_length);
//        }
//    }
//
//    public void estimateThetaGenome(int[][] genes, int window_length) {
//
//        int no_genes = genes.length;
//
//        for(int g=0; g<no_genes; g++) {
//
//            System.out.println(this.genes[g]);
//            int[] gene_boundary = genes[g];
//            int gene_start = gene_boundary[0]-1;
//            int gene_end = gene_boundary[1]-1;
//            int gene_length = gene_end-gene_start;
//            int no_sites = gene_length/window_length;
//
//            estimateThetaPerWindow(gene_start, gene_end, window_length);
//
//        }
//
//    }
//
//    public void estimateMultiThetaGenome(int[][] genes, int window_length) {
//
//        int no_genes = genes.length;
//
//        for(int g=0; g<no_genes; g++) {
//            System.out.println(this.genes[g]);
//            int[] gene_boundary = genes[g];
//            int gene_start = gene_boundary[0]-1;
//            int gene_end = gene_boundary[1]-1;
//            int gene_length = gene_end-gene_start;
//            int no_sites = gene_length/window_length;
//
//            estimateMultiThetaPerWindow(gene_start, gene_end, window_length);
//        }
//
//
//
//    }


    //test from SUPGMA paper
    private double[] multipleRegression() {

        double [] y = new double [] {0.0271,0.0582,0.0512,0.0317,0.0547,0.005,0.0153,0.0693,0.0327,0.0875,0.0584,0.0089,0.0383,0.0736,0.0352};

        double [][] x = new double[y.length][];

//        x[0] = new double[]{1,0,0,0,0};
//        x[1] = new double[]{0,1,0,0,0};
//        x[2] = new double[]{0,0,1,0,0};
//        x[3] = new double[]{0,0,0,1,0};
//        x[4] = new double[]{0,0,0,1,0};
//        x[5] = new double[]{0,0,0,1,0};
//        x[6] = new double[]{0,0,0,1,0};
//        x[7] = new double[]{0,0,0,1,1};
//        x[8] = new double[]{0,0,0,1,1};
//        x[9] = new double[]{0,0,0,1,1};
//        x[10] = new double[]{0,0,0,1,1};
//        x[11] = new double[]{0,0,0,0,1};
//        x[12] = new double[]{0,0,0,0,1};
//        x[13] = new double[]{0,0,0,0,1};
//        x[14] = new double[]{0,0,0,0,1};

        x[0] = new double[]{1,0,0};
        x[1] = new double[]{1,0,0};
        x[2] = new double[]{1,0,0};
        x[3] = new double[]{1,1,0};
        x[4] = new double[]{1,1,0};
        x[5] = new double[]{1,1,0};
        x[6] = new double[]{1,1,0};
        x[7] = new double[]{1,1,1};
        x[8] = new double[]{1,1,1};
        x[9] = new double[]{1,1,1};
        x[10] = new double[]{1,1,1};
        x[11] = new double[]{1,0,1};
        x[12] = new double[]{1,0,1};
        x[13] = new double[]{1,0,1};
        x[14] = new double[]{1,0,1};

        OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();

//        regression.newSampleData(y, x);
        regression.setNoIntercept(true);
        regression.newSampleData(y, x);
        double[] beta = regression.estimateRegressionParameters();

        for(double d: beta) {
            System.out.println(d);
        }

        return beta;
    }

//    private void fillNodeMaps(Map<Integer, String> nodeNames, Map<Integer, Double> nodeHeights){
//        //, List<Integer> sampledLineages, List<Double> sampledTimes) {
//
//        for(int i=0; i<sampleNames.size(); i++) {
//            nodeNames.put(i+1, "'"+sampleNames.get(i)+"'");
//            nodeHeights.put(i+1, samplingTimes.get(i));
//
//        }
//    }
//
//    public void getUPGMA(double[][] distmatrix) {
//
//        //while(distmatrix.length>1) {
//
//        String newick = "";
//        Map<Integer, String> nodeNames = new HashMap<Integer, String>();
//        Map<Integer, Double> nodeHeights = new HashMap<Integer, Double>();
//        fillNodeMaps(nodeNames, nodeHeights);
//
//
//        int n_lineages = nodeNames.size();
//        int k = 0;
//        String treeString = "";
//
//
//        //while(distmatrix.length>1) {
//
//        double min = 100;
//        int taxon_i = -1;
//        int taxon_j = -1;
//
//
//        String node1 = "";
//        String node2 = "";
//
//        for (int i = 0; i < distmatrix.length; i++) {
//
//
//            for (int j = 0; j < distmatrix.length - 1; j++) {
//
//                System.out.println(distmatrix[i][j]);
//
//                if (distmatrix[i][j] < min) {
//                    min = distmatrix[i][j];
//                    taxon_i = i;
//                    taxon_j = j;//
//
//                    node1 = sampleNames.get(i);
//                    node2 = sampleNames.get(j+1);
//
//                }
//
//            }
//
//        }
//        String dist = String.valueOf(min / 2);
//        treeString = "("+node1+":"+dist+","+node2+":"+dist+")";
//        nodeNames.put(n_lineages+k+1, treeString);
//        nodeHeights.put(n_lineages+k+1,Double.valueOf(dist));
//
//        double[][] newmat = new double[distmatrix.length-1][distmatrix.length-2];
//        int new_i = 0;
//        for (int i = 0; i < distmatrix.length; i++) {
//
//
//
//            int new_j = 0;
//            for (int j = 1; j < distmatrix.length - 1; j++) {
//
//
//                //if(!(i==taxon_i && j==taxon_j)) {
//
//                System.out.println(i+","+j+","+distmatrix[taxon_i][j]+","+distmatrix[taxon_j][j]);
//
//
//
//
//                newmat[i-1][new_j] = (distmatrix[taxon_i][j]+distmatrix[taxon_i][j])/2;
//                //System.out.println(new_i+","+new_j);
//                new_j++;
//                //}
//
//
//
//            }
//
//        }
//
//        for (int i = 0; i < newmat.length; i++) {
//
//            String s="";
//            for (int j = 0; j < newmat.length-1; j++) {
//
//                s+=newmat[i][j]+",";
//
//            }
//            System.out.println(s);
//        }
//
//        //}
//
//
//
//
// }



    public static void main(String [] args) throws IOException, ImportException {

        SUPGMA s = new SUPGMA("/Users/jayna/Documents/Projects/HIV_ANPI/original_fasta/processed_fasta_by_group/15513_1#61-67_aln_C_processed.fasta", 4);

        s.estimateMultiThetaPerWindow(295,1575, 300);
        s.estimateMultiThetaPerWindow(1776,4520, 300);
        s.estimateMultiThetaPerWindow(5801,7855, 300);
        s.estimateMultiThetaPerWindow(8295,8585, 300);

//        HKYDistanceMatrix distanceMatrix = (HKYDistanceMatrix)s.getDistances('H');
//
//        List<Double> samplingTimes = s.samplingTimes;
//
//        double[] beta = s.estimateMultiTheta(distanceMatrix, samplingTimes);
//        double[][] updatedDistances = s.updateDistanceMatrix(beta, distanceMatrix);
//
//
//        double[][] originalDistances = distanceMatrix.getDistances();
//
//        double[][] testDistances = new double[][]{{0.0,0.2,0.3,0.6,0.6,0.8}, {0.2,0.0,0.3,0.6,0.6,0.8}, {0.3,0.3,0.0,0.5,0.5,0.7}, {0.6,0.6,0.7,0.0,0.3,0.8}, {0.6,0.6,0.7,0.3,0.0,0.8}, {0.8,0.8,0.7,0.8,0.8,0.0}};
//
//        //new double [6]{new double[]{0.0,0.2,0.4,0.6,0.6,0.8}, new double[]{0.2,0.0,0.4,0.6,0.6,0.8}, new double[]{}, new double[]{}, new double[]{}, new double[]{}};
////        for(int i=0; i<originalDistances.length; i++) {
////
////            String ss = "";
////            for (int j = 0; j < originalDistances.length; j++) {
////
////                ss += originalDistances[i][j] + ",";
////
////            }
////            System.out.println(ss);
////        }
//        BasicDistanceMatrix basicDistanceMatrix = new BasicDistanceMatrix(distanceMatrix.getTaxa(), updatedDistances);
//
//
//        ClusteringTreeBuilder treeBuilder = ClusteringTreeBuilder.getBuilder(TreeBuilderFactory.Method.UPGMA, basicDistanceMatrix);
////
////
//        SimpleRootedTree tree = (SimpleRootedTree)treeBuilder.build();
//
//
//
//        Set<Node> externalNodes = tree.getExternalNodes();
//        Set<Node> internalNodes = tree.getInternalNodes();
//
//        internalNodes.remove(tree.getRootNode());
//        double maxDistance = 0.0;
//        for(Node n: externalNodes) {
//
//            int index = distanceMatrix.getTaxa().indexOf(tree.getTaxon(n));
//            Double sampleTime = samplingTimes.get(index);
//
//            Double samplingInterval = Collections.max(samplingTimes) - sampleTime;
//
//            double distance = samplingInterval * beta[beta.length - 1];
//
//
//            if (maxDistance < distance) {
//                maxDistance = distance;
//            }
//
//
//            double currentDepth = tree.getHeight(tree.getRootNode())-tree.getHeight(n);
//            double currentParentDepth = tree.getHeight(tree.getRootNode())-tree.getHeight(tree.getParent(n));
//            double currentLength = currentDepth-currentParentDepth;
//
//            double newDepth = currentDepth-distance;
//
//
//            double newLength = tree.getLength(n) - distance;
//
//
//
//
//            tree.setLength(n, newLength);
//            if(tree.getLength(n) < 0) {
//
//
//                double newHeight = tree.getHeight(n);
//
//                s.fixHeight(tree, tree.getParent(n), tree.getHeight(n));
//
//            }
//
//
//        }





//        BufferedWriter writer = new BufferedWriter(new FileWriter("/Users/jayna/Documents/Projects/AMC_HCV_DATA/fasta_nt/random_subsampled_datasets/sub_004_subsampled_test1_SUPGMA.tre"));
//        NexusExporter exporter = new NexusExporter(writer);
//        exporter.exportTree(tree);
//        writer.close();

    }

    private void fixHeight(SimpleRootedTree tree, Node node, double height) {


        tree.setHeight(node, height);

        if (!tree.isRoot(node)) {
            if (tree.getHeight((tree.getParent(node))) < 0) {

                System.out.println();
                fixHeight(tree, tree.getParent(node), height);

            }

        }
    }
}
