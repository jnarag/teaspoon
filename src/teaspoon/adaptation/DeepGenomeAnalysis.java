package teaspoon.adaptation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import teaspoon.app.utils.MainAlignmentParser;
import teaspoon.app.utils.TeaspoonMethods;
import teaspoon.app.utils.TeaspoonValues;

/**
 * Created by jayna on 11/03/16.
 */
public class DeepGenomeAnalysis implements Analysis {

    String ancestralFile;
    String mainFile;
    TeaspoonValues[][] value_matrix;
    int no_timepoints;
    int no_datasets;
    int bootstraps;

    //deep genome MainAnalysis
    int window_length = 3;
    boolean fixedNR;
    double[] nr;
    double[][] bins;
    double[] L = {0.0, 0.15, 0.75};
    double[] H = {0.15, 0.75, 1.0};
    boolean[] Nvec = {false,true,false};
    int[][]  genes;
    double[] prior = {1.0,1.0,1.0,1.0};
    Map<String, String[]> timepoints_multi;
    String[] timepoints;
    int[] timepoints_per_dataset;
    double[] firstTimepoint;
    String[] datasets = new String[no_datasets];
    int[] map;
    Map<String, Integer> which;

    /**
     * 	
     * @param ancestralFile
     * @param mainFile
     */
    public DeepGenomeAnalysis(String ancestralFile, String mainFile) {

        this.ancestralFile = ancestralFile;
        this.mainFile = mainFile;
        timepoints_multi = new HashMap<String, String[]>();
        nr = new double[no_datasets];
        bins = new double[2][L.length];
        
        for(int i=0;i<L.length;i++){
            bins[0][i]=L[i];
            bins[1][i]=H[i];
        }
    }

    @Override
    public void bmAnalysis() {

        MainAlignmentParser ancestral;
        List<int[][]> main_alignments = new ArrayList<int[][]>();    //list of main_alignments that represent different timepoints...

        no_timepoints = timepoints.length;
        try {
            BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
            ancestral = new MainAlignmentParser(reader1.readLine().trim(),true);
            ancestral.readFASTA();

            BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));
            String output = mainFile;

            BufferedWriter highfreq = null;
            BufferedWriter lowfreq = null;
            BufferedWriter midfreq = null;

            if (fixedNR) {
                //output+="_"+"useFixedNeutralRatio";
                output = output.replace(".txt", "_fixedNR");
                System.out.println(output);
                highfreq = new BufferedWriter(new FileWriter(output + "_highfreq_table.csv"));
                lowfreq = new BufferedWriter(new FileWriter(output + "_lowfreq_table.csv"));
                midfreq = new BufferedWriter(new FileWriter(output + "_midfreq_table.csv"));
            } else {
                output = output.replace(".txt", "_indivNR");
                System.out.println(output);
            }

            while(reader2.ready()) {
                String filename = reader2.readLine().trim();
                MainAlignmentParser m = new MainAlignmentParser(filename, true);
                System.out.println(filename);
                main_alignments.add(m.readFASTA());
            }

            reader1.close();
            reader2.close();
            
            BufferedWriter summaryResults = new BufferedWriter(new FileWriter(output + "_summaryTable_windowlength_"+window_length+".csv"));
            StringBuffer mid = new StringBuffer();
            StringBuffer low = new StringBuffer();
            StringBuffer high = new StringBuffer();
            StringBuilder summary = new StringBuilder();
            StringBuilder neutralRatio = new StringBuilder();

            mid.append("gene,time,window,total_sites_mid,no_silent_sites_mid,no_replacement_sites_mid,r_mid/s_mid,no_of_adaptations\n");
            low.append("gene,time,window,total_sites_low,no_silent_sites_low,no_replacement_sites_low,r_low/s_low,no_of_noneutral_sites\n");
            high.append("gene,time,window,total_sites_high,no_silent_sites_high,no_replacement_sites_high,r_high/s_high,no_of_adaptations\n");
            summary.append("Gene,Timepoint,r_l,r_m,r_h,s_l,s_m,s_h,a_l,a_h,Neutral Ratio,no_of_windows\n");
            neutralRatio.append("Gene,r_m,s_m,Average Neutral Ratio\n");

            int no_genes = genes.length;

            for(int g=0; g<no_genes; g++) {

                int[] gene_boundary = genes[g];
                int gene_start = gene_boundary[0]-1;
                int gene_end = gene_boundary[1];
                int gene_length = gene_end-gene_start;
                int no_sites = (int)Math.round(Math.ceil(gene_length/(double)window_length));

                System.out.println(gene_start+","+gene_end+": "+no_sites);

                DescriptiveStatistics r_m = new DescriptiveStatistics();
                DescriptiveStatistics s_m = new DescriptiveStatistics();

                for (int t = 0; t < no_timepoints; t++) {

                    DescriptiveStatistics r_l = new DescriptiveStatistics();
                    DescriptiveStatistics s_l = new DescriptiveStatistics();
                    DescriptiveStatistics r_h = new DescriptiveStatistics();
                    DescriptiveStatistics s_h = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_l = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_m = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_h = new DescriptiveStatistics();

                    //System.out.println("*****************"+timepoints[t]+"*****************");
                    System.out.println("t:\t" + timepoints[t]);
                    int c = 0; // no of windows with no seqs
                    double d = 0.0;
                    int codon_start = gene_start;
                    int codon_end = codon_start + window_length;
                    
                    for (int i = 0; i < no_sites; i++) {

                        // check if sequences in frame
//                        codon_start += window_length;
//                        codon_end += window_length;

                        if(codon_end > gene_end) {
                            codon_end = gene_end;
                            d = (codon_end-codon_start)/(double)window_length;
                        }
                        
                        //System.out.println(codon_start + "," + codon_end);
                        int[][] site_ans = TeaspoonMethods.subMatrix(ancestral.getSequenceMatrix(), codon_start, codon_end, false);

                        if (site_ans.length > 0) {

                            int[] site_ans_con = ancestral.consensusArray(site_ans);

                            // sub_alignment corresponding to the window - only includes taxa where reads have 95 % coverage
                            int[][] site_main = TeaspoonMethods.subMatrix(main_alignments.get(t), codon_start, codon_end, true);

                            //making sure there are more than 100 sequences in the alignment
                            if (site_main.length >= 100) {

                                //System.out.println(i+","+timepoints[t]+",main_length = 0");
                                BhattMethod bm = new BhattMethod(site_main, site_ans_con);

                                if (fixedNR == true) {
                                    bm.Method(bins, prior, true, Nvec, nr[g]);
                                } else {
                                    bm.Method(bins, prior, true, Nvec);
                                }

                                for (int r = 0; r < bm.getReplacementSubstitutionsCountArray().length; r++) {
                                    if (Double.isNaN(bm.getReplacementSubstitutionsCountArray()[r])) {
                                        bm.getReplacementSubstitutionsCountArray()[r] = 0.0;
                                    }
                                }

                                for (int s = 0; s < bm.getSilentSubstitutionsCountArray().length; s++) {
                                    if (Double.isNaN(bm.getSilentSubstitutionsCountArray()[s])) {
                                        bm.getSilentSubstitutionsCountArray()[s] = 0.0;
                                    }
                                }

                                for (int s = 0; s < bm.getNonNeutralSubstitutions().length; s++) {
                                    if (Double.isNaN(bm.getNonNeutralSubstitutions()[s])) {
                                        bm.getNonNeutralSubstitutions()[s] = 0.0;
                                    }
                                }


                                TeaspoonMethods.record(low, datasets[g], new double[]{g, i, Double.parseDouble(timepoints[t]), 0}, bm);
                                TeaspoonMethods.record(mid, datasets[g], new double[]{g, i, Double.parseDouble(timepoints[t]), 1}, bm);
                                TeaspoonMethods.record(high, datasets[g], new double[]{g, i, Double.parseDouble(timepoints[t]), 2}, bm);

                                if (!Double.isNaN(bm.getReplacementSubstitutionsCountArray()[1])) {
                                    r_m.addValue(bm.getReplacementSubstitutionsCountArray()[1]);
                                }
 
                                if (!Double.isNaN(bm.getSilentSubstitutionsCountArray()[1])) {
                                    s_m.addValue(bm.getSilentSubstitutionsCountArray()[1]);
                                }

                                r_l.addValue(bm.getReplacementSubstitutionsCountArray()[0]);
                                s_l.addValue(bm.getSilentSubstitutionsCountArray()[0]);
//                              r_m.addValue(bm.ReplacementCountArray[1]);
//                              s_m.addValue(bm.SilentCountArray[1]);
                                r_h.addValue(bm.getReplacementSubstitutionsCountArray()[2]);
                                s_h.addValue(bm.getSilentSubstitutionsCountArray()[2]);
                                adapt_l.addValue(bm.getNonNeutralSubstitutions()[0]);
                                adapt_m.addValue(bm.getNonNeutralSubstitutions()[1]);
                                adapt_h.addValue(bm.getNonNeutralSubstitutions()[2]);

                            } else {
                                System.out.println("no seqs at window " + i);
                                if(d>0) {
                                    d = 1;
                                }
                                c++;
                            }
                        } else {
                            System.out.println("no seqs at first timepoint at window " + i);
                            if(d>0) {
                                d = 1;
                            }
                            c++;
                        }
                        codon_start += window_length;
                        codon_end += window_length;
                    }
                    high.append("\n");
                    low.append("\n");
                    mid.append("\n");

                    System.out.println("no of windows: "+(double)((no_sites-1)+d-c));

                    //System.out.println(gene_names[g]+","+timepoints[t]+","+r_h.getSum()+","+adapt_h.getSum()+","+(no_sites-c));
                    summary.append(datasets[g]).append(",").append(timepoints[t]).append(",").append(r_l.getSum()).append(",").append(r_m.getSum()).append(",").append(r_h.getSum()).append(",").append(s_l.getSum()).append(",").append(s_m.getSum()).append(",").append(s_h.getSum()).append(",").append(adapt_l.getSum()).append(",").append(adapt_h.getSum()).append(",").append(r_m.getSum()/s_m.getSum()).append(",").append((no_sites -1-c) + d).append("\n");
                }
                neutralRatio.append(datasets[g]).append(",").append(r_m.getSum()).append(",").append(s_m.getSum()).append(",").append(r_m.getSum() / s_m.getSum()).append("\n");
                System.out.println(">" + datasets[g] + ": r_m = " + r_m.getSum() + ", s_m = " + s_m.getSum() + " average_nr = " + r_m.getSum() / s_m.getSum());
            }

            if(fixedNR) {
                lowfreq.write(low.toString());
                midfreq.write(mid.toString());
                highfreq.write(high.toString());

                highfreq.close();
                lowfreq.close();
                midfreq.close();
            }

            summaryResults.write(summary.toString());
            summaryResults.write("\n\n");
            summaryResults.write(neutralRatio.toString());
            summaryResults.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    @Override
    /*
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#bmAnalysisBootstrap(int)
     */
    public void bmAnalysisBootstrap(int bootstraps) {

        this.bootstraps = bootstraps;
        value_matrix = new TeaspoonValues[timepoints.length][datasets.length];

        for(int i=0; i<value_matrix.length; i++){
            for(int j = 0; j < value_matrix[0].length; j++) {
                value_matrix[i][j] = new TeaspoonValues(bootstraps);
            }
        }

        MainAlignmentParser ancestral = null;
        List<int[][]> main_alignments = new ArrayList<int[][]>();

        try {
            BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
            BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));
            ancestral = new MainAlignmentParser(reader1.readLine().trim(), true);
            ancestral.readFASTA();

            while (reader2.ready()) {
                String filename = reader2.readLine().trim();
                MainAlignmentParser m = new MainAlignmentParser(filename, true);
                System.out.println(filename);
                main_alignments.add(m.readFASTA());
            }            
            reader1.close();
            reader2.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        Map<String, DescriptiveStatistics> dataStatistics = new HashMap<>();

        for(String dataset: datasets) {
            dataStatistics.put(dataset, new DescriptiveStatistics());
        }

        for( int bs=0; bs < bootstraps; bs++) {
            System.out.println("I am on Run  "+(bs+1)+"  of "+bootstraps);

            //list of main_alignments that represent different timepoints...

            no_timepoints = timepoints.length;

            int no_genes = genes.length;

            for (int geneIndex = 0; geneIndex < no_genes; geneIndex++) {

                int[] gene_boundary = genes[geneIndex];
                int gene_start = gene_boundary[0]-1;
                int gene_end = gene_boundary[1];
                int gene_length = gene_end-gene_start;
                //int no_sites = gene_length / window_length;
                int no_sites = (int)Math.round(Math.ceil(gene_length/(double)window_length));

                System.out.println(gene_start + "," + (gene_end) + ": " + no_sites);

                int[] sampler = new int[gene_length / 3];
                //choosing your numCodons
                RandomGenerator generator = new RandomGenerator();
                for (int x = 0; x < sampler.length; x++) {
                    sampler[x] = generator.nextInt(sampler.length - 1);//generator.nextInt(sampler.length-1);
                }

                assert ancestral != null;
                int[][] site_ans = TeaspoonMethods.subMatrix(ancestral.getSequenceMatrix(), gene_start, gene_end, false);

                double[] x_time_points = new double[no_timepoints+1];
                double[] y_adapt_points = new double[no_timepoints+1];
                SimpleRegression regression = new SimpleRegression(true);
                regression.addData(firstTimepoint[geneIndex], 0);

                for (int timepointIndex = 0; timepointIndex < no_timepoints; timepointIndex++) {
                    DescriptiveStatistics r_m = new DescriptiveStatistics();
                    DescriptiveStatistics s_m = new DescriptiveStatistics();
                    DescriptiveStatistics r_l = new DescriptiveStatistics();
                    DescriptiveStatistics s_l = new DescriptiveStatistics();
                    DescriptiveStatistics r_h = new DescriptiveStatistics();
                    DescriptiveStatistics s_h = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_l = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_m = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_h = new DescriptiveStatistics();

                    // TODO refactor c and d to be sensibly named
                    // NOTE used elsewhere for another variable, be careful of scope!!
                    int c = 0; // number of windows with 'no' sequences (e.g. depth < critical, usually '100')

                    int[][] site_main = TeaspoonMethods.subMatrix(main_alignments.get(timepointIndex), gene_start, gene_end, false);

                    //if (site_ans.length > 0) {

                    int[] site_ans_con = ancestral.consensusArray(site_ans);
                    BhattMethod b = new BhattMethod(site_main, site_ans_con);
                    Store s = b.CreateBlocks(3, site_ans_con.length, sampler); //******
                    //BhattMethod bm = new BhattMethod(s.RandomisedIntegerMatrix, s.RandomisedIntegerAncestral);

                    double d = 0.0;
                    int codon_start = 0;
                    int codon_end = codon_start + window_length;   // check if sequences in frame

                    for (int i = 0; i < no_sites; i++) {

                    	if(codon_end > gene_end) {
                            codon_end = gene_end;
                            d = (codon_end-codon_start)/(double)window_length;
                        }

                        int[] temp_ans = Arrays.copyOfRange(s.RandomisedIntegerAncestral, codon_start, codon_end);

                        if (temp_ans.length > 0) {

                            int[][] temp_main = TeaspoonMethods.subMatrix(s.RandomisedIntegerMatrix, codon_start, codon_end, true);

                            if (temp_main.length >= 100) {

                                BhattMethod bm = new BhattMethod(temp_main, temp_ans);

                                bm.Method(bins, prior, true, Nvec, nr[geneIndex]);

                                for (int x = 0; x < bm.getReplacementSubstitutionsCountArray().length; x++) {
                                    if (Double.isNaN(bm.getReplacementSubstitutionsCountArray()[x])) {
                                        bm.getReplacementSubstitutionsCountArray()[x] = 0.0;
                                    }
                                }

                                for (int x = 0; x < bm.getSilentSubstitutionsCountArray().length; x++) {
                                    if (Double.isNaN(bm.getSilentSubstitutionsCountArray()[x])) {
                                        bm.getSilentSubstitutionsCountArray()[x] = 0.0;
                                    }
                                }

                                for (int x = 0; x < bm.getNonNeutralSubstitutions().length; x++) {
                                    if (Double.isNaN(bm.getNonNeutralSubstitutions()[x])) {
                                        bm.getNonNeutralSubstitutions()[x] = 0.0;
                                    }
                                }

                                if (Double.isNaN(bm.Adaptation)) {
                                    bm.Adaptation = 0.0;
                                }

                                r_l.addValue(bm.getReplacementSubstitutionsCountArray()[0]);
                                s_l.addValue(bm.getSilentSubstitutionsCountArray()[0]);
                                r_m.addValue(bm.getReplacementSubstitutionsCountArray()[1]);
                                s_m.addValue(bm.getSilentSubstitutionsCountArray()[1]);
                                r_h.addValue(bm.getReplacementSubstitutionsCountArray()[2]);
                                s_h.addValue(bm.getSilentSubstitutionsCountArray()[2]);
                                adapt_l.addValue(bm.getNonNeutralSubstitutions()[0]);
                                adapt_m.addValue(bm.getNonNeutralSubstitutions()[1]);
                                adapt_h.addValue(bm.Adaptation);

                            } else {
                                System.out.println("no seqs at window " + i);
                                if(d>0) {
                                    d = 1;
                                }
                                c++;
                            }
                        } else {
                            System.out.println("no seqs at window " + i);
                            if(d>0) {
                                d = 1;
                            }
                            c++;
                        }
                        codon_start += window_length;
                        codon_end += window_length;
                    }
                    
                    System.out.println("no of windows: "+(double)((no_sites-1)+d-c));

                    value_matrix[timepointIndex][geneIndex].getReplacementPolymorphsMidClass()[bs] = r_m.getSum() ;/// (((no_sites-1)+d-c) * (window_length / 3));
                    value_matrix[timepointIndex][geneIndex].getSilentPolymorphsMidClass()[bs] = s_m.getSum() ;/// (((no_sites-1)+d-c) * (window_length / 3));
                    value_matrix[timepointIndex][geneIndex].getReplacementSubsHigh()[bs] = r_h.getSum() ;/// (((no_sites-1)+d-c) * (window_length / 3));
                    value_matrix[timepointIndex][geneIndex].getSilentSubsHigh()[bs] = s_h.getSum() ;/// (((no_sites-1)+d-c) * (window_length / 3));
                    value_matrix[timepointIndex][geneIndex].getReplacementPolymorphsLowClass()[bs] = r_l.getSum() ;/// (((no_sites-1)+d-c) * (window_length / 3));
                    value_matrix[timepointIndex][geneIndex].getSilentPolymorphsLowClass()[bs] = s_l.getSum() ;/// (((no_sites-1)+d-c) * (window_length / 3));
                    value_matrix[timepointIndex][geneIndex].getAdaptations()[bs] = adapt_h.getSum() ;/// ((no_sites - c) * (window_length / 3));
                    value_matrix[timepointIndex][geneIndex].setRow(timepoints[timepointIndex]);
                    value_matrix[timepointIndex][geneIndex].setColumn(datasets[geneIndex]);
                    value_matrix[timepointIndex][geneIndex].getNumWindows()[bs] = (no_sites+1 - c);

                    x_time_points[timepointIndex+1]=Double.parseDouble(timepoints[timepointIndex]);
                    y_adapt_points[timepointIndex+1]=adapt_h.getSum();// / ((no_sites - c) * (window_length / 3));

                    if(x_time_points[timepointIndex+1]>0) {
                        regression.addData(Double.parseDouble(timepoints[timepointIndex]), y_adapt_points[timepointIndex+1]);
                    }
                }

                double slope = regression.getSlope();
                dataStatistics.get(datasets[geneIndex]).addValue(slope);
            }
            System.out.println("I am on Run  " + (bs + 1) + "  of " + bootstraps);
        }

        System.out.println();
        System.out.println("************Adaptive substitutions per year************\n");

        for(String aDataset: datasets) {
            System.out.println(aDataset+", "+dataStatistics.get(aDataset).getMean()+", "+dataStatistics.get(aDataset).getPercentile(25)+", "+dataStatistics.get(aDataset).getPercentile(75));
        }

        System.out.println();
        String output = mainFile;

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(output.replace(".txt", "_bootstraps_n_"+bootstraps+"_adaptations.csv")));
            writer.write("dataset, Time,Mean,Median,Lower Quartile,Upper Quartile,Standard Deviation,No of Codons\n");

            for (int datasetIndex = 0; datasetIndex < no_datasets; datasetIndex++) {

            	String s = String.valueOf(datasets[datasetIndex]);
                writer.write(s+","+firstTimepoint[datasetIndex] + "," + 0 + ","+ 0 + "," + 0 + "," + 0 + "," + 0 + "," + 0+"\n");

                for (int t = 0; t < no_timepoints; t++) {

                    DescriptiveStatistics bstraps = new DescriptiveStatistics(value_matrix[t][datasetIndex].getAdaptations());

                    double mean = bstraps.getMean();
                    double lq = bstraps.getPercentile(25);
                    double uq = bstraps.getPercentile(75);
                    double med = bstraps.getPercentile(50);
                    double std = bstraps.getStandardDeviation();

                    if(value_matrix[t][datasetIndex].getNumCodons()>0) {
                        System.out.println(s+","+timepoints[t] + "," + mean + ","+ med + "," + lq + "," + uq + "," + std + "," + value_matrix[t][datasetIndex].getNumCodons());
                        writer.write(s + "," + timepoints[t] + "," + mean + "," + med + "," + lq + "," + uq + "," + std + "," + value_matrix[t][datasetIndex].getNumCodons() + "\n");
                    }
                }
                writer.write("\n");
                System.out.println();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void w3bAnalysis() {

    }

    @Override
    public void w3bAnalysisBootstrap(int bootstraps) {

    }

    @Override
    public void getDiversityStats() {

        this.bootstraps = bootstraps;
        value_matrix = new TeaspoonValues[timepoints.length][datasets.length];
        for (int i = 0; i < value_matrix.length; i++) {
            for (int j = 0; j < value_matrix[0].length; j++) {
                value_matrix[i][j] = new TeaspoonValues(bootstraps);
            }
        }

        MainAlignmentParser ancestral = null;
        List<int[][]> main_alignments = new ArrayList<int[][]>();

        try {
            BufferedReader ancestralReader = new BufferedReader(new FileReader(ancestralFile));
            BufferedReader mainReader = new BufferedReader(new FileReader(mainFile));
            ancestral = new MainAlignmentParser(ancestralReader.readLine().trim(), true);
            ancestral.readFASTA();

            while (mainReader.ready()) {
                String filename = mainReader.readLine().trim();
                MainAlignmentParser m = new MainAlignmentParser(filename, true);
                System.out.println(filename);
                main_alignments.add(m.readFASTA());
            }
            ancestralReader.close();
            mainReader.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        for (int bootstrapIndex = 0; bootstrapIndex < bootstraps; bootstrapIndex++) {
            System.out.println("I am on Run  " + (bootstrapIndex + 1) + "  of " + bootstraps);

            //list of main_alignments that represent different timepoints...

            no_timepoints = timepoints.length;

            int no_genes = genes.length;

            for (int geneIndex = 0; geneIndex < no_genes; geneIndex++) {
                int[] gene_boundary = genes[geneIndex];
                int gene_start = gene_boundary[0] - 1;
                int gene_end = gene_boundary[1] - 1;
                int gene_length = gene_end - gene_start;
                int no_sites = gene_length / window_length;
                System.out.println(gene_start + "," + gene_end + ": " + no_sites);

                for (int timepointIndex = 0; timepointIndex < no_timepoints; timepointIndex++) {
 
                	DescriptiveStatistics r_m = new DescriptiveStatistics();
                    DescriptiveStatistics s_m = new DescriptiveStatistics();
                    DescriptiveStatistics r_l = new DescriptiveStatistics();
                    DescriptiveStatistics s_l = new DescriptiveStatistics();
                    DescriptiveStatistics r_h = new DescriptiveStatistics();
                    DescriptiveStatistics s_h = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_l = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_m = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_h = new DescriptiveStatistics();
                    int c = 0;
                    
                    for (int siteIndex = 0; siteIndex < no_sites; siteIndex++) {

                        int codon_start = gene_start;
                        int codon_end = codon_start + window_length;   // check if sequences in frame
                        codon_start += siteIndex * window_length;
                        codon_end += siteIndex * window_length;

                        //System.out.println(codon_start + "," + codon_end);
                        assert ancestral != null;
                        int[][] site_ans = TeaspoonMethods.subMatrix(ancestral.getSequenceMatrix(), codon_start, codon_end, false);

                        if (site_ans.length > 0) {

                            int[] site_ans_con = ancestral.consensusArray(site_ans);

                            //System.out.println("*****************"+timepoints[t]+"*****************");
                            //System.out.println("t:\t" + timepoints[t]);
                            //main_alignments.get(t).readFASTA();

                            int[][] site_main = TeaspoonMethods.subMatrix(main_alignments.get(timepointIndex), codon_start, codon_end, true);

                            int[] sampler = new int[site_ans_con.length / 3];
                            //choosing your numCodons
                            RandomGenerator generator = new RandomGenerator();

                            for (int x = 0; x < sampler.length; x++) {
                                sampler[x] = generator.nextInt(sampler.length - 1);//generator.nextInt(sampler.length-1);
                            }

                            if (site_main.length > 100) {

                                BhattMethod b = new BhattMethod(site_main, site_ans_con);
                                Store s = b.CreateBlocks(3, site_main[0].length, sampler); //******
                                BhattMethod bm = new BhattMethod(s.RandomisedIntegerMatrix, s.RandomisedIntegerAncestral);
                                bm.Method(bins, prior, true, Nvec, nr[geneIndex]);

                                for (int x = 0; x < bm.getReplacementSubstitutionsCountArray().length; x++) {
                                    if (Double.isNaN(bm.getReplacementSubstitutionsCountArray()[x])) {
                                        bm.getReplacementSubstitutionsCountArray()[x] = 0.0;
                                    }
                                }

                                for (int x = 0; x < bm.getSilentSubstitutionsCountArray().length; x++) {
                                    if (Double.isNaN(bm.getSilentSubstitutionsCountArray()[x])) {
                                        bm.getSilentSubstitutionsCountArray()[x] = 0.0;
                                    }
                                }

                                for (int x = 0; x < bm.getNonNeutralSubstitutions().length; x++) {
                                    if (Double.isNaN(bm.getNonNeutralSubstitutions()[x])) {
                                        bm.getNonNeutralSubstitutions()[x] = 0.0;
                                    }
                                }

                                if (Double.isNaN(bm.Adaptation)) {
                                    bm.Adaptation = 0.0;
                                }

                                r_l.addValue(bm.getReplacementSubstitutionsCountArray()[0]);
                                s_l.addValue(bm.getSilentSubstitutionsCountArray()[0]);
                                r_m.addValue(bm.getReplacementSubstitutionsCountArray()[1]);
                                s_m.addValue(bm.getSilentSubstitutionsCountArray()[1]);
                                r_h.addValue(bm.getReplacementSubstitutionsCountArray()[2]);
                                s_h.addValue(bm.getSilentSubstitutionsCountArray()[2]);
                                adapt_l.addValue(bm.getNonNeutralSubstitutions()[0]);
                                adapt_m.addValue(bm.getNonNeutralSubstitutions()[1]);
                                adapt_h.addValue(bm.Adaptation);

                            } else {
                                System.out.println("no seqs at window " + siteIndex);
                                c++;
                            }
                        } else {
                            System.out.println("no seqs at first timepoint at window " + siteIndex);
                            c++;
                        }
                    }

                    System.out.println("no of windows: " + (no_sites - c));

                    value_matrix[timepointIndex][geneIndex].getReplacementPolymorphsMidClass()[bootstrapIndex] = r_m.getSum();
                    value_matrix[timepointIndex][geneIndex].getSilentPolymorphsMidClass()[bootstrapIndex] = s_m.getSum();
                    value_matrix[timepointIndex][geneIndex].getReplacementSubsHigh()[bootstrapIndex] = r_h.getSum();
                    value_matrix[timepointIndex][geneIndex].getSilentSubsHigh()[bootstrapIndex] = s_h.getSum();
                    value_matrix[timepointIndex][geneIndex].getReplacementPolymorphsLowClass()[bootstrapIndex] = r_l.getSum();
                    value_matrix[timepointIndex][geneIndex].getSilentPolymorphsLowClass()[bootstrapIndex] = s_l.getSum();
                    value_matrix[timepointIndex][geneIndex].getAdaptations()[bootstrapIndex] = adapt_h.getSum();
                    value_matrix[timepointIndex][geneIndex].setRow(timepoints[timepointIndex]);
                    value_matrix[timepointIndex][geneIndex].setColumn(datasets[geneIndex]);
                    value_matrix[timepointIndex][geneIndex].getNumWindows()[bootstrapIndex] = (no_sites - c);
                }
            }
        }
    }

    @Override
    public TeaspoonValues[][] getValue_matrix() {
        return this.value_matrix;
    }

    @Override
    public String[] getDatasets() {
        return this.datasets;
    }

    @Override
    public int[] getTimepoints_per_dataset() {
        return this.timepoints_per_dataset;
    }

    public int getBootstraps() {
        return this.bootstraps;
    }

    public void setCodonMap(int[] map){
        this.map = map;
    }

    public void setWhichMap(Map<String, Integer> which) {
        this.which = which;
    }

    public void getBootstrapsBybmAnalysisByWindow(int bootstraps) {

        this.bootstraps = bootstraps;
        value_matrix = new TeaspoonValues[timepoints.length][datasets.length];
        
        for (int i = 0; i < value_matrix.length; i++) {
            for (int j = 0; j < value_matrix[0].length; j++) {
                value_matrix[i][j] = new TeaspoonValues(bootstraps);
            }
        }

        MainAlignmentParser ancestral = null;
        List<int[][]> main_alignments = new ArrayList<int[][]>();

        try {
            BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
            BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));
            ancestral = new MainAlignmentParser(reader1.readLine().trim(), true);
            ancestral.readFASTA();

            while (reader2.ready()) {
                String filename = reader2.readLine().trim();
                MainAlignmentParser m = new MainAlignmentParser(filename, true);
                System.out.println(filename);
                main_alignments.add(m.readFASTA());
            }
            
            reader1.close();
            reader2.close();
            
        } catch (FileNotFoundException e) {
            e.printStackTrace();  
        } catch (IOException e) {
            e.printStackTrace();  
        }

        for (int bs = 0; bs < bootstraps; bs++) {
        	
            System.out.println("I am on Run  " + (bs + 1) + "  of " + bootstraps);
            //list of main_alignments that represent different timepoints...
            no_timepoints = timepoints.length;
            int no_genes = genes.length;

            for (int g = 0; g < no_genes; g++) {

                int[] gene_boundary = genes[g];
                int gene_start = gene_boundary[0] - 1;
                int gene_end = gene_boundary[1] - 1;
                int gene_length = gene_end - gene_start;
                int no_sites = gene_length / window_length;

                System.out.println(gene_start + "," + gene_end + ": " + no_sites);

                for (int t = 0; t < no_timepoints; t++) {
                    DescriptiveStatistics r_m = new DescriptiveStatistics();
                    DescriptiveStatistics s_m = new DescriptiveStatistics();
                    DescriptiveStatistics r_l = new DescriptiveStatistics();
                    DescriptiveStatistics s_l = new DescriptiveStatistics();
                    DescriptiveStatistics r_h = new DescriptiveStatistics();
                    DescriptiveStatistics s_h = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_l = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_m = new DescriptiveStatistics();
                    DescriptiveStatistics adapt_h = new DescriptiveStatistics();

                    int c = 0;
  
                    for (int i = 0; i < no_sites; i++) {

                        int codon_start = gene_start;
                        int codon_end = codon_start + window_length;   // check if sequences in frame
                        codon_start += i * window_length;
                        codon_end += i * window_length;

                        //System.out.println(codon_start + "," + codon_end);
                        assert ancestral != null;
                        int[][] site_ans = TeaspoonMethods.subMatrix(ancestral.getSequenceMatrix(), codon_start, codon_end, false);

                        if (site_ans.length > 0) {

                            int[] site_ans_con = ancestral.consensusArray(site_ans);
                            //System.out.println("*****************"+timepoints[t]+"*****************");
                            //System.out.println("t:\t" + timepoints[t]);
                            //main_alignments.get(t).readFASTA();

                            int[][] site_main = TeaspoonMethods.subMatrix(main_alignments.get(t), codon_start, codon_end, true);
                            int[] sampler = new int[site_ans_con.length / 3];
                            
                            //choosing your numCodons
                            RandomGenerator generator = new RandomGenerator();

                            for (int x = 0; x < sampler.length; x++) {
                                sampler[x] = generator.nextInt(sampler.length - 1);//generator.nextInt(sampler.length-1);
                            }

                            if (site_main.length > 10) {
                                //System.out.println(i+","+timepoints[t]+",main_length = 0");

                                BhattMethod b = new BhattMethod(site_main, site_ans_con);
                                Store s = b.CreateBlocks(3, site_main[0].length, sampler); //******
                                BhattMethod bm = new BhattMethod(s.RandomisedIntegerMatrix, s.RandomisedIntegerAncestral);

                                //bm.Method(bins, prior, false, Nvec, neutral_ratio[g]);
                                bm.Method(bins, prior, true, Nvec, nr[g]);
//                                if (useFixedNeutralRatio == true) {
//                                	bm.Method(bins, prior, true, Nvec, nr[g]);
//                                } else {
//                                	bm.Method(bins, prior, true, Nvec);
//                                }

                                for (int x = 0; x < bm.getReplacementSubstitutionsCountArray().length; x++) {
                                    if (Double.isNaN(bm.getReplacementSubstitutionsCountArray()[x])) {
                                        bm.getReplacementSubstitutionsCountArray()[x] = 0.0;
                                    }
                                }

                                for (int x = 0; x < bm.getSilentSubstitutionsCountArray().length; x++) {
                                    if (Double.isNaN(bm.getSilentSubstitutionsCountArray()[x])) {
                                        bm.getSilentSubstitutionsCountArray()[x] = 0.0;
                                    }
                                }

                                for (int x = 0; x < bm.getNonNeutralSubstitutions().length; x++) {
                                    if (Double.isNaN(bm.getNonNeutralSubstitutions()[x])) {
                                        bm.getNonNeutralSubstitutions()[x] = 0.0;
                                    }
                                }

                                if (Double.isNaN(bm.Adaptation)) {
                                    bm.Adaptation = 0.0;
                                }

                                //System.out.println(t+","+timepoints[t])
//                                    record(low, new double[]{g, i, Double.parseDouble(timepoints[t]), 0}, bm);
//                                    record(mid, new double[]{g, i, Double.parseDouble(timepoints[t]), 1}, bm);
//                                    record(high, new double[]{g, i, Double.parseDouble(timepoints[t]), 2}, bm);

                                r_l.addValue(bm.getReplacementSubstitutionsCountArray()[0]);
                                s_l.addValue(bm.getSilentSubstitutionsCountArray()[0]);
                                r_m.addValue(bm.getReplacementSubstitutionsCountArray()[1]);
                                s_m.addValue(bm.getSilentSubstitutionsCountArray()[1]);
                                r_h.addValue(bm.getReplacementSubstitutionsCountArray()[2]);
                                s_h.addValue(bm.getSilentSubstitutionsCountArray()[2]);
                                adapt_l.addValue(bm.getNonNeutralSubstitutions()[0]);
                                adapt_m.addValue(bm.getNonNeutralSubstitutions()[1]);
                                adapt_h.addValue(bm.Adaptation);

//                                if (!Double.isNaN(bm.ReplacementCountArray[1])) {
//                                	r_m.addValue(bm.ReplacementCountArray[1]);
//                                }
//                                if (!Double.isNaN(bm.SilentCountArray[1])) {
//                                	s_m.addValue(bm.SilentCountArray[1]);
//                                }

                            } else {
                                System.out.println("no seqs at window " + i);
                                c++;
                            }
                        } else {
                            System.out.println("no seqs at first timepoint at window " + i);
                            c++;
                        }
                    }
                    
                    System.out.println("no of windows: " + (no_sites - c));
//                        high.append("\n");
//                        low.append("\n");
//                        mid.append("\n");
                    value_matrix[t][g].getReplacementPolymorphsMidClass()[bs] = r_m.getSum();
                    value_matrix[t][g].getSilentPolymorphsMidClass()[bs] = s_m.getSum();
                    value_matrix[t][g].getReplacementSubsHigh()[bs] = r_h.getSum();
                    value_matrix[t][g].getSilentSubsHigh()[bs] = s_h.getSum();
                    value_matrix[t][g].getReplacementPolymorphsLowClass()[bs] = r_l.getSum();
                    value_matrix[t][g].getSilentPolymorphsLowClass()[bs] = s_l.getSum();
                    value_matrix[t][g].getAdaptations()[bs] = adapt_h.getSum();
                    value_matrix[t][g].setRow(timepoints[t]);
                    value_matrix[t][g].setColumn(datasets[g]);
                    value_matrix[t][g].getNumWindows()[bs] = (no_sites - c);

                }
                //System.out.println(">" + gene_names[g] + ": r_m = " + r_m.getSum() + ", s_m = " + s_m.getSum() + " average_nr = " + r_m.getSum() / s_m.getSum());
            }
        }
    }


//    public static void main(String [] args) {
//
////        MainAnalysis bs = new MainAnalysis("//Users/jayna/Documents/Projects/AMC_HCV_DATA/ancestral_HCVpacbio_filelist.txt",
////                "/Users/jayna/Documents/Projects/AMC_HCV_DATA/main_HCVpacbio_filelist.txt");
//
////        MainAnalysis bs = new MainAnalysis("/Users/jayna/Documents/Projects/MTE_HIV/ancestral_file.txt",
////                  "/Users/jayna/Documents/Projects/MTE_HIV/main_file.txt");
//
//
//        //deep genome MainAnalysis
//
//        List<DeepGenomeAnalysis> deepGenomeList = new ArrayList<>();
//
//        List<int[][]> gene_coordinates = new ArrayList<>();
//
//        List<double[]> nr_list = new ArrayList<>();
//
//        List<String[]> timepoints_list = new ArrayList<>();
//
//
//
//        gene_coordinates.add(new int[][]{{294,1574},{1775,4519},{5800,7875},{8315,8605}});
//        gene_coordinates.add(new int[][]{{295,1578},{1779,4523},{5811,7811},{8251,8541}});
//        gene_coordinates.add(new int[][]{{294,1574},{1775,4519},{5800,7842},{8282,8572}});
//        gene_coordinates.add(new int[][]{{293,1573},{1774,4518},{5799,7808},{8247,8537}});
//        gene_coordinates.add(new int[][]{{295,1575},{1776,4520},{5801,7855},{8295,8585}});
//        gene_coordinates.add(new int[][]{{294,1574},{1778,4519},{5800,7806},{8246,8536}});
//        gene_coordinates.add(new int[][]{{295,1578},{1797,4541},{5835,7919},{8338,8640}});
//        gene_coordinates.add(new int[][]{{295,1575},{1797,4541},{5822,7846},{8286,8576}});
//
////        gene_coordinates.add(new int[][]{{294,1574},{1775,4519},{5800,7875},{8315,8938}});
////        gene_coordinates.add(new int[][]{{295,1578},{1779,4523},{5811,7811},{8251,8952}});
////        gene_coordinates.add(new int[][]{{294,1574},{1775,4519},{5800,7842},{8282,8983}});
////        gene_coordinates.add(new int[][]{{293,1573},{1774,4518},{5799,7808},{8247,8870}});
////        gene_coordinates.add(new int[][]{{295,1575},{1776,4520},{5801,7855},{8295,8918}});
////        gene_coordinates.add(new int[][]{{294,1574},{1778,4519},{5800,7806},{8246,8869}});
////        gene_coordinates.add(new int[][]{{295,1578},{1797,4541},{5835,7919},{8338,8973}});
////        gene_coordinates.add(new int[][]{{295,1575},{1797,4541},{5822,7846},{8286,8909}});
//
//
//        timepoints_list.add(new String[]{"2008.6192"});
//        timepoints_list.add(new String[]{"2005.9534","2008.5178","2011.0301"});
//        timepoints_list.add(new String[]{"2006.5644","2008.7479","2011.4932"});
//        timepoints_list.add(new String[]{"2009.7699"});
//        timepoints_list.add(new String[]{"2004.9726","2007.1945","2009.8027","2011.1644"});
//        timepoints_list.add(new String[]{"2009.5342"});
//        timepoints_list.add(new String[]{"2010.0301"});
//        timepoints_list.add(new String[]{"2010.1452"});
//
//
//        nr_list.add(new double[] {0.59754419,0.422687674,4.69888515,0.596970359});
//        nr_list.add(new double[] {0.695251345,0.487694516,2.169390943,2.126155712});
//        nr_list.add(new double[] {0.383574361,0.625901721,2.590430972,0.68590778});
//        nr_list.add(new double[] {0.467184508,0.544443491,4.477416339,0.505629913});
//        nr_list.add(new double[] {0.972356426,0.782089794,2.482752792,3.700956057});
//        nr_list.add(new double[] {0.766413237,3.178117435,2.75673462,135.4243521});
//        nr_list.add(new double[] {1.304859458,0.740904326,3.070492898,2.129520815});
//        nr_list.add(new double[] {0.193130206,0.305294462,4.402938816,0.009999977});
//
////        nr_list.add(new double[] {0.59754419,0.422687674,4.69888515,0.549805231});
////        nr_list.add(new double[] {0.695251345,0.487694516,2.169390943,2.007446765});
////        nr_list.add(new double[] {0.383574361,0.625901721,2.590430972,1.128349745});
////        nr_list.add(new double[] {0.467184508,0.544443491,4.477416339,1.235627431});
////        nr_list.add(new double[] {0.972356426,0.782089794,2.482752792,2.284056361});
////        nr_list.add(new double[] {0.766413237,3.178117435,2.75673462,7.961451291});
////        nr_list.add(new double[] {1.304859458,0.740904326,3.070492898,1.790799242});
////        nr_list.add(new double[] {0.193130206,0.305294462,4.402938816,0.015501522});
//
//        deepGenomeList.add(new DeepGenomeAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS081.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS081.txt"));
//        deepGenomeList.add(new DeepGenomeAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS114.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS114.txt"));
//        deepGenomeList.add(new DeepGenomeAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS133.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS133.txt"));
//        deepGenomeList.add(new DeepGenomeAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS364.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS364.txt"));
//        deepGenomeList.add(new DeepGenomeAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS380.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS380.txt"));
//        deepGenomeList.add(new DeepGenomeAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS468.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS468.txt"));
//        deepGenomeList.add(new DeepGenomeAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS559.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS559.txt"));
//        deepGenomeList.add(new DeepGenomeAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS385.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS385.txt"));
//
//
//
//        for(int i=0; i< deepGenomeList.size(); i++) {
//
//            //int i = 1;
//            DeepGenomeAnalysis bs = deepGenomeList.get(i);
//
//            bs.genes = gene_coordinates.get(i);
//            bs.timepoints = timepoints_list.get(i);
//            bs.fixedNR = true;
//            bs.window_length = 3;
//
//
//
//            bs.neutral_ratio = nr_list.get(i);
//            bs.fixedNR = true;
//            //bs.runDeepGenomeAnalysis();
////            bs.runBootstrapDeepGenomeAnalysis(100);
////            String outfile1 = bs.mainFile;
////
////
////            outfile1 = outfile1.replace("main_file_P","P");
////            outfile1 = outfile1.replace(".txt","_adapt_bs.csv");
////            bs.writeoutBootstrapResults(outfile1,"a");
////
////
////            String outfile2 = outfile1.replace("_adapt", "_rh");
////            bs.writeoutBootstrapResults(outfile2,"r_high");
////
////            String outfile3 = outfile2.replace("_rh", "_sh");
////            bs.writeoutBootstrapResults(outfile3,"s_high");
//
//            //System.out.println(outfile1);
//            //System.out.println(outfile2);
//        }
//
//
////      MainAnalysis bs = new MainAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS081.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS081.txt");
////        MainAnalysis bs = new MainAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS114.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS114.txt");
////        MainAnalysis bs = new MainAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS133.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS133.txt");
////        MainAnalysis bs = new MainAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS364.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS364.txt");
////        MainAnalysis bs = new MainAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS380.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS380.txt");
////        MainAnalysis bs = new MainAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS468.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS468.txt");
////        MainAnalysis bs = new MainAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS559.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS559.txt");
////        MainAnalysis bs = new MainAnalysis("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/ancestral_file_PS385.txt", "/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/main_file_PS385.txt");
//
////        bs.runDeepGenomeAnalysis();
////        bs.runBootstrapDeepGenomeAnalysis(100);
//
////        bs.writeoutBootstrapResults("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/PS133_adapt_bs.csv","a");
////        bs.writeoutBootstrapResults("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/PS133_rh_bs.csv","r_high");
//
//        //bs.writeoutBootstrapResults("/Users/jayna/Documents/Projects/HIV_deep_genome/seqs/adapt_bs.csv","a");
//
//
//    }
}
