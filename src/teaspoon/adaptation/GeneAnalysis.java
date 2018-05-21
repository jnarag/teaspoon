package teaspoon.adaptation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
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
public class GeneAnalysis implements Analysis {

    String ancestralFile = "";
    String mainFile = "";
    TeaspoonValues[][] value_matrix;
    int no_timepoints;
    int no_datasets;
    int bootstraps;

    double[] nr;
    double[] sm;
    double[] rm;
    double[] L = {0.0, 0.15, 0.75};
    double[] H = {0.15, 0.75, 1.0};
    double[] low = {0.0, 0.15};
    double[] mid = {0.15, 0.75};
    double[] high = {0.75, 1.0};
    boolean[] Nvec = {false,true,false};
    double[] prior = {1.0,1.0,1.0,1.0};
    double[][] bins;
    //this could be no of site for next-gen data
    String[] timepoints;
    int[] timepoints_per_dataset;
    String[] datasets = new String[no_datasets];
    double[] firstTimepoint;
    int n_seq_limit = 10;

    Map<String, String[]> timepoints_multi;
    boolean fixedNR;
    int[] gene_boundary = new int[] {445-1,924};
    boolean byGene = false;
    int [] map;
    Map<String, Integer> which = new HashMap<String, Integer>();

    /**
     * 
     * @param ancestralFile
     * @param mainFile
     */
    public GeneAnalysis(String ancestralFile, String mainFile) {

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

    /**
     * Analyse using Bhatt method, no bootstrap
     * @param ancestorAverage
     */
    public void bmAnalysis(boolean ancestorAverage) {
        try {
            BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
            BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));

            String output = mainFile;

            if (fixedNR) {
                //output+="_"+"useFixedNeutralRatio";
                output = output.replace(".txt", "_fixedNR");
            } else {
                output = output.replace(".txt", "_indivNR");
                System.out.println(output);
            }

            BufferedWriter highfreq = new BufferedWriter(new FileWriter(output + "_highfreq_table.csv"));
            BufferedWriter lowfreq = new BufferedWriter(new FileWriter(output + "_lowfreq_table.csv"));
            BufferedWriter midfreq = new BufferedWriter(new FileWriter(output + "_midfreq_table.csv"));

            StringBuffer sb1 = new StringBuffer();
            StringBuffer mid = new StringBuffer();
            StringBuffer low = new StringBuffer();
            StringBuffer high = new StringBuffer();

            sb1.append("patient,timepoint,range,total_sites,no_silent_sites,no_replacement_sites,Replacement/Silent Ratio,No_of_NonNeutral_changes\n");
            mid.append("gene,time,total_sites_mid,no_silent_sites_mid,no_replacement_sites_mid,r_mid/s_mid,no_of_adaptations\n");
            low.append("gene,time,total_sites_low,no_silent_sites_low,no_replacement_sites_low,r_mid/s_mid,no_of_noneutral_sites\n");
            high.append("gene,time,total_sites_high,no_silent_sites_high,no_replacement_sites_high,r_mid/s_mid,no_of_adaptations\n");

            for (int t = 0; t < no_datasets; t++) {
            	
                no_timepoints = timepoints_per_dataset[t];
                if (reader1.ready()) {

                    String ancesfilename = reader1.readLine().trim();
                    MainAlignmentParser ra = new MainAlignmentParser(ancesfilename);

                    int[][] ancestralMatrix = ra.readNEXUS();

                    DescriptiveStatistics r_m = new DescriptiveStatistics();
                    DescriptiveStatistics s_m = new DescriptiveStatistics();

                    DescriptiveStatistics ans_stats = new DescriptiveStatistics();
                    for (int d = 0; d < no_timepoints; d++) {

                        if (reader2.ready()) {

                        	String mainFile = reader2.readLine().trim();
                            MainAlignmentParser r = new MainAlignmentParser(mainFile);

                            int[][] main = r.readNEXUS();
                            for (int a = 0; a < ancestralMatrix.length; a++) {

                                int[] ans = ancestralMatrix[a];
                                int[] ans_tmp = ans;

                                if (which.containsKey(datasets[t])) {
                                    System.out.println(which.get(datasets[t]));
                                    int[][] tmp_1 = TeaspoonMethods.subsetter(main, map, which.get(datasets[t]));
                                    int[] tmp_2 = TeaspoonMethods.subsetter(ans_tmp, map, which.get(datasets[t]));
                                    main = tmp_1;
                                    ans_tmp = tmp_2;
                                }

                                BhattMethod bm = new BhattMethod(main, ans_tmp);

                                if (fixedNR == true) {

                                    bm.Method(bins, prior, true, Nvec, nr[t]);
                                    double r_h = bm.getReplacementSubstitutionsCountArray()[2];
                                    double s_h = bm.getSilentSubstitutionsCountArray()[2];
//                                adaptations = r_h - (neutral_ratio[t]*s_h)*(1+1.0/(s_h));
                                    //(1.0+(1.0/s_mid[t]));
                                    //*
                                } else {
                                    bm.Method(bins, prior, true, Nvec);
                                }
                                //System.out.println(d + ": " + mainFile + " : " + main[0].length + " neutralratio:" + bm.neutralratio);

                                TeaspoonMethods.record(low, datasets[t], new double[]{t, Double.parseDouble(timepoints[d]), 0}, bm);
                                TeaspoonMethods.record(mid, datasets[t], new double[]{t, Double.parseDouble(timepoints[d]), 1}, bm);
                                TeaspoonMethods.record(high, datasets[t], new double[]{t, Double.parseDouble(timepoints[d]), 2}, bm);

                                if (!Double.isNaN(bm.getReplacementSubstitutionsCountArray()[1])) {
                                    r_m.addValue(bm.getReplacementSubstitutionsCountArray()[1]);
                                }

                                if (!Double.isNaN(bm.getSilentSubstitutionsCountArray()[1])) {
                                    s_m.addValue(bm.getSilentSubstitutionsCountArray()[1]);
                                }
                               
                                if (!Double.isNaN(bm.getReplacementSubstitutionsCountArray()[0]) && !Double.isNaN(bm.getSilentSubstitutionsCountArray()[0])) {
                                    ans_stats.addValue(bm.getReplacementSubstitutionsCountArray()[0]+bm.getSilentSubstitutionsCountArray()[0]);
                                }
                                
//                           sb1.append("\n");
                                //System.out.println(">" + datasets[t] + ": r_m = " + r_m.getSum() + ", s_m = " + s_m.getSum() + " average_nr = " + r_m.getSum() / s_m.getSum());
                                high.append("\n");
                                low.append("\n");
                                mid.append("\n");
                            }
                        }
                        System.out.println(datasets[t]+","+ans_stats.getMean()+","+Math.pow(ans_stats.getVariance(), 0.5));
                    }
                }
            }

            lowfreq.write(low.toString());
            midfreq.write(mid.toString());
            highfreq.write(high.toString());

            reader1.close();
            reader2.close();
            highfreq.close();
            lowfreq.close();
            midfreq.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();  
        } catch (IOException e) {
            e.printStackTrace();  
        }
    }

    /**
     * Analyse using Bhatt method, no bootstrap
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#bmAnalysis()
     */
    public void bmAnalysis(){

        no_datasets = datasets.length;

        try {

            BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
            BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));

            String output = mainFile;

            if (fixedNR) {
                //output+="_"+"useFixedNeutralRatio";
                output = output.replace(".txt", "_fixedNR");
            } else {
                output = output.replace(".txt", "_indivNR");
                System.out.println(output);
            }
            BufferedWriter highfreq = new BufferedWriter(new FileWriter(output + "_highfreq_table.csv"));
            BufferedWriter lowfreq = new BufferedWriter(new FileWriter(output + "_lowfreq_table.csv"));
            BufferedWriter midfreq = new BufferedWriter(new FileWriter(output + "_midfreq_table.csv"));

            StringBuffer sb1 = new StringBuffer();
            StringBuffer mid = new StringBuffer();
            StringBuffer low = new StringBuffer();
            StringBuffer high = new StringBuffer();

            sb1.append("patient,timepoint,range,total_sites,no_silent_sites,no_replacement_sites,Replacement/Silent Ratio,No_of_NonNeutral_changes\n");
            mid.append("gene,time,total_sites_mid,no_silent_sites_mid,no_replacement_sites_mid,r_mid/s_mid,no_of_adaptations\n");
            low.append("gene,time,total_sites_low,no_silent_sites_low,no_replacement_sites_low,r_low/s_low,no_of_noneutral_sites\n");
            high.append("gene,time,total_sites_high,no_silent_sites_high,no_replacement_sites_high,r_high/s_high,no_of_adaptations\n");

            for (int dd = 0; dd < no_datasets; dd++) {

                no_timepoints = timepoints_per_dataset[dd];
                timepoints = timepoints_multi.get(datasets[dd]);

                System.out.println(no_timepoints);
                if (reader1.ready()) {

                    String ancesfilename = reader1.readLine().trim();
                    MainAlignmentParser ra = new MainAlignmentParser(ancesfilename);

                    int[][] ancestralMatrix = ra.readFASTA();

                    DescriptiveStatistics r_m = new DescriptiveStatistics();
                    DescriptiveStatistics s_m = new DescriptiveStatistics();

                    for (int tt = 0; tt < no_timepoints; tt++) {

                        if (reader2.ready()) {
                            String mainFile = reader2.readLine().trim();

                            MainAlignmentParser r = new MainAlignmentParser(mainFile);

                            int[][] main = r.readFASTA();

                            while((main == null || main.length < 10) && reader2.ready()) {

                                mainFile = reader2.readLine().trim();
                                r = new MainAlignmentParser(mainFile);
                                main = r.readFASTA();
                                tt+=1;

                                if(tt == no_timepoints) {
                                    break;
                                }
                            }

                            if(byGene) {
                                ancestralMatrix = ra.subMatrix(gene_boundary[0], gene_boundary[1]);
                                main = r.subMatrix(gene_boundary[0], gene_boundary[1]);
                            }


                            int[] ans_tmp = ra.consensusArray(ancestralMatrix);

                            if (which.containsKey(datasets[dd])) {
                                int[][] tmp_1 = TeaspoonMethods.subsetter(main, map, which.get(datasets[dd]));
                                int[] tmp_2 = TeaspoonMethods.subsetter(ans_tmp, map, which.get(datasets[dd]));
                                main = tmp_1;
                                ans_tmp = tmp_2;
                            }

                            BhattMethod bm = new BhattMethod(main, ans_tmp);

                            if (fixedNR) {
                                bm.Method(bins, prior, true, Nvec, nr[dd]);
                            } else {
                                bm.Method(bins, prior, true, Nvec);
                            }

                            System.out.println(tt + ": " + mainFile + " : " + main[0].length + " : "+main.length+" neutralratio:" + bm.neutralratio);
                            TeaspoonMethods.record(low, datasets[dd],new double[]{dd, Double.parseDouble(timepoints[tt]), 0}, bm);
                            TeaspoonMethods.record(mid, datasets[dd],new double[]{dd, Double.parseDouble(timepoints[tt]), 1}, bm);
                            TeaspoonMethods.record(high, datasets[dd],new double[]{dd, Double.parseDouble(timepoints[tt]), 2}, bm);

                            if (!Double.isNaN(bm.getReplacementSubstitutionsCountArray()[1])) {
                                r_m.addValue(bm.getReplacementSubstitutionsCountArray()[1]);
                            }
                            
                            if (!Double.isNaN(bm.getSilentSubstitutionsCountArray()[1])) {
                                s_m.addValue(bm.getSilentSubstitutionsCountArray()[1]);
                            }
                        }
                    }

                    System.out.println(">" + datasets[dd] + ": r_m = " + r_m.getSum() + ", s_m = " + s_m.getSum() + " average_nr = " + r_m.getSum() / s_m.getSum());
                    high.append("\n");
                    low.append("\n");
                    mid.append("\n");
                }
            }

            lowfreq.write(low.toString());
            midfreq.write(mid.toString());
            highfreq.write(high.toString());

            reader1.close();
            reader2.close();
            highfreq.close();
            lowfreq.close();
            midfreq.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    
    @Override
    /**
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#bmAnalysisBootstrap(int)
     */
    public void bmAnalysisBootstrap(int bootstraps) {

        this.bootstraps = bootstraps;
        value_matrix = new TeaspoonValues[this.no_timepoints][this.no_datasets];

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
            ancestral = new MainAlignmentParser(reader1.readLine().trim(), false);
            ancestral.readFASTA();
            while (reader2.ready()) {
                String filename = reader2.readLine().trim();
                MainAlignmentParser m = new MainAlignmentParser(filename, false);
                System.out.println(filename);
                main_alignments.add(m.readFASTA());
            }
            reader1.close();
            reader2.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        assert ancestral != null;
        int[][] ancestralMatrix = ancestral.getSequenceMatrix();
        int[] ans = ancestral.consensusArray(ancestralMatrix);

//        DescriptiveStatistics regression_PRD = new DescriptiveStatistics();
//        DescriptiveStatistics regression_YRD1 = new DescriptiveStatistics();
//        DescriptiveStatistics regression_YRD22 = new DescriptiveStatistics();
//        DescriptiveStatistics regression_YRD21 = new DescriptiveStatistics();

        Map<String, DescriptiveStatistics> dataStatistics = new HashMap<>();

        for(String d: datasets) {
            dataStatistics.put(d, new DescriptiveStatistics());
        }

        for (int bs = 0; bs < bootstraps; bs++) {
            //choosing your numCodons
            RandomGenerator generator = new RandomGenerator();
            int[] sampler = new int[ans.length / 3];

            for (int x = 0; x < sampler.length; x++) {
                sampler[x] = generator.nextInt(sampler.length - 1);//generator.nextInt(sampler.length-1);
            }

            int c = 0;

            for (int t = 0; t < no_datasets; t++) {
                int [] ans_tmp = ans;
//                System.out.println(neutral_ratio[t]);

                no_timepoints = timepoints_per_dataset[t];
                timepoints = timepoints_multi.get(datasets[t]);

                double[] x_time_points = new double[no_timepoints+1];
                double[] y_adapt_points = new double[no_timepoints+1];

                if (which.containsKey(datasets[t])) {

                    int[] tmp_2 = TeaspoonMethods.subsetter(ans_tmp, map, which.get(datasets[t]));

                    ans_tmp = tmp_2;
                    sampler = new int[ans_tmp.length / 3];

                    for (int x = 0; x < sampler.length; x++) {
                        sampler[x] = generator.nextInt(sampler.length - 1);//generator.nextInt(sampler.length-1);
                    }
                }

                SimpleRegression regression = new SimpleRegression(true);
                regression.addData(firstTimepoint[t], 0); //2013.33 should be first timepoint...
                //System.out.println(2013.33+","+0);
                for (int d = 0; d < no_timepoints; d++) {

                    int[][] main = main_alignments.get(c);
                    //if(reader2.ready()) {
                    //String mainFile = reader2.readLine().trim();
                    //MainAlignmentParser r = new MainAlignmentParser(mainFile);
                    //int [] ans_tmp = ancestral.consensusArray(ancestralMatrix);
                    //System.out.println(c);
                    x_time_points[0] = firstTimepoint[t];

                    //System.out.println(d+","+main);
                    while((main == null || main.length < 10)) {
                        d+=1;
                        c++;
                        if(d == no_timepoints && c== main_alignments.size()) {
                            break;
                        }
                        //System.out.println(main);
                        main = main_alignments.get(c);
                        //System.out.println(">"+main);
                        //System.out.println(c);
                    }
                    //System.out.println(c);
                    //System.out.println(d+","+main);

                    if (which.containsKey(datasets[t])) {

                        //int[] list = codonlist.get(datasets[t]);
                        int[][] tmp_1 = TeaspoonMethods.subsetter(main, map, which.get(datasets[t]));
                        //int[] tmp_2 = methods.Subsetter(ans_tmp, map, which.get(datasets[t]));

                        main = tmp_1;
                        //ans_tmp = tmp_2;
                    }

                    //int[] tmp_2 = ancestral.consensusArray(ancestralMatrix);
                    //System.out.println(main);
                    BhattMethod b = new BhattMethod(main, ans_tmp); //******
                    Store s = b.CreateBlocks(3, main[0].length, sampler); //******
                    BhattMethod bm = new BhattMethod(s.RandomisedIntegerMatrix, s.RandomisedIntegerAncestral);

                    bm.Method(bins, prior, true, Nvec, nr[t]);

                    //System.out.println(timepoints[d]+","+datasets[t]);
                    value_matrix[d][t].setRow(timepoints[d]);
                    value_matrix[d][t].setColumn(datasets[t]);
                    value_matrix[d][t].setNumCodons(sampler.length);
                    value_matrix[d][t].getReplacementPolymorphsMidClass()[bs] = bm.getReplacementSubstitutionsCountArray()[1];
                    value_matrix[d][t].getSilentPolymorphsMidClass()[bs] = bm.getSilentSubstitutionsCountArray()[1];
                    value_matrix[d][t].getReplacementSubsHigh()[bs] = bm.getReplacementSubstitutionsCountArray()[2];
                    value_matrix[d][t].getSilentSubsHigh()[bs] = bm.getSilentSubstitutionsCountArray()[2];
                    value_matrix[d][t].getSilentPolymorphsLowClass()[bs] = bm.getSilentSubstitutionsCountArray()[0];
                    value_matrix[d][t].getReplacementPolymorphsLowClass()[bs] = bm.getReplacementSubstitutionsCountArray()[0];

                    if (bm.getReplacementSubstitutionsCountArray()[2] > 0) {
                        value_matrix[d][t].getSilentToReplacementPolymorphsHighClassRatio()[bs] = bm.getSilentSubstitutionsCountArray()[2] / bm.getReplacementSubstitutionsCountArray()[2];
                    }

                    if (bm.getSilentSubstitutionsCountArray()[2] > 0) {
                        value_matrix[d][t].getReplacementToSilentPolymorphsHighClassRatio()[bs] = bm.getReplacementSubstitutionsCountArray()[2] / bm.getSilentSubstitutionsCountArray()[2];
                    }

                    if (bm.getSilentSubstitutionsCountArray()[0] > 0) {
                        value_matrix[d][t].getReplacementToSilentPolymorphsLowClassRatio()[bs] = bm.getReplacementSubstitutionsCountArray()[0] / bm.getSilentSubstitutionsCountArray()[0];
                    }
                    
                    if (bm.getSilentSubstitutionsCountArray()[1] > 0) {
                        value_matrix[d][t].getNeutralRatio()[bs] = bm.neutralratio;
                    }

                    if (Double.isNaN(bm.Adaptation)) {
                        bm.Adaptation = 0.0;
                    }
                    
                    value_matrix[d][t].getAdaptations()[bs] = bm.Adaptation;
                    DiversityStats diversityStats = new DiversityStats(bm.integer_matrix);
                    double[] wattersonEstimates = diversityStats.wattersonEstimates();
                    value_matrix[d][t].getTajimasD()[bs] = diversityStats.TajimasD();

                    if (Double.isNaN(value_matrix[d][t].getTajimasD()[bs])) {
                        value_matrix[d][t].getTajimasD()[bs] = 0.0;
                    }

                    value_matrix[d][t].getWattersonS()[bs] = wattersonEstimates[0];
                    value_matrix[d][t].getWattersonPi()[bs] = wattersonEstimates[1];
                    value_matrix[d][t].getTheta()[bs] = diversityStats.theta();

                    c++;

                    x_time_points[d+1]=Double.parseDouble(timepoints[d]);
                    y_adapt_points[d+1]=bm.Adaptation;

                    if(x_time_points[d+1]>0) {
                        regression.addData(Double.parseDouble(timepoints[d]), y_adapt_points[d+1]);
                    }
                }
                double slope = regression.getSlope();
                dataStatistics.get(datasets[t]).addValue(slope);
            }
            System.out.println("I am on Run  " + (bs + 1) + "  of " + bootstraps);
        }

        System.out.println();
        System.out.println("************Adaptive substitutions per year************\n");

        for(String d: datasets) {
            System.out.println(d+", "+dataStatistics.get(d).getMean()+", "+dataStatistics.get(d).getPercentile(25)+", "+dataStatistics.get(d).getPercentile(75));
        }
//        System.out.println(+regression_PRD.getMean()+", "+regression_PRD.getPercentile(25)+", "+regression_PRD.getPercentile(75));
//        System.out.println("YRD1 "+regression_YRD1.getMean()+", "+regression_YRD1.getPercentile(25)+", "+regression_YRD1.getPercentile(75));
//        System.out.println("YRD2.2 "+regression_YRD22.getMean()+", "+regression_YRD22.getPercentile(25)+", "+regression_YRD22.getPercentile(75));
//        System.out.println("YRD2.1 "+regression_YRD21.getMean()+", "+regression_YRD21.getPercentile(25)+", "+regression_YRD21.getPercentile(75));

        System.out.println();
        String output = mainFile;

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(output.replace(".txt", "_bootstraps_n_"+bootstraps+"_adaptations.csv")));
            writer.write("dataset, Time,Mean,Median,Lower Quartile,Upper Quartile,Standard Deviation,No of Codons\n");

            for (int d = 0; d < no_datasets; d++) {

            	String s = String.valueOf(datasets[d]);
                writer.write(s+","+firstTimepoint[d] + "," + 0 + ","+ 0 + "," + 0 + "," + 0 + "," + 0 + "," + 0+"\n");

                for (int t = 0; t < no_timepoints; t++) {
                	
                    DescriptiveStatistics bstraps = new DescriptiveStatistics(value_matrix[t][d].getAdaptations());
                    double mean = bstraps.getMean();
                    double lq = bstraps.getPercentile(25);
                    double uq = bstraps.getPercentile(75);
                    double med = bstraps.getPercentile(50);
                    double std = bstraps.getStandardDeviation();

                    if(value_matrix[t][d].getNumCodons()>0) {
                        System.out.println(s+","+timepoints[t] + "," + mean + ","+ med + "," + lq + "," + uq + "," + std + "," + value_matrix[t][d].getNumCodons());
                        writer.write(s + "," + timepoints[t] + "," + mean + "," + med + "," + lq + "," + uq + "," + std + "," + value_matrix[t][d].getNumCodons() + "\n");
                    }
//                System.out.println();
//                System.out.println(d + "," + t + "," + value_matrix[t][d].r_high[d]);
//                for (int b = 0; b < 10; b++) {
//                    s += "," + value_matrix[t][d].adaptations[b];
//                }
                }
                writer.write("\n");
                System.out.println();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#w3bAnalysis()
     */
    public void w3bAnalysis() {

        no_datasets = datasets.length;

        try {
            BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
            BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));
            String output = mainFile.replace(".txt", "_output.txt");

            if (fixedNR) {
                //output+="_"+"useFixedNeutralRatio";
                output = output.replace(".txt", "_fixedNR");
            } else {
                output = output.replace(".txt", "_indivNR");
                System.out.println(output);
            }

            BufferedWriter highfreq = new BufferedWriter(new FileWriter(output + "_highfreq_table.csv"));
            BufferedWriter lowfreq = new BufferedWriter(new FileWriter(output + "_lowfreq_table.csv"));
            BufferedWriter midfreq = new BufferedWriter(new FileWriter(output + "_midfreq_table.csv"));

            StringBuffer sb1 = new StringBuffer();
            StringBuffer mid = new StringBuffer();
            StringBuffer low = new StringBuffer();
            StringBuffer high = new StringBuffer();

            sb1.append("patient,timepoint,range,total_sites,no_silent_sites,no_replacement_sites,Replacement/Silent Ratio,No_of_NonNeutral_changes\n");
            mid.append("gene,time,total_sites_mid,no_silent_sites_mid,no_replacement_sites_mid,r_mid/s_mid,no_of_adaptations\n");
            low.append("gene,time,total_sites_low,no_silent_sites_low,no_replacement_sites_low,r_low/s_low,no_of_noneutral_sites\n");
            high.append("gene,time,total_sites_high,no_silent_sites_high,no_replacement_sites_high,r_high/s_high,no_of_adaptations\n");

//            BufferedWriter writer4 = new BufferedWriter(new FileWriter(output));
//            StringBuilder sb1 = new StringBuilder();
//            StringBuilder sb2 = new StringBuilder();
//            StringBuilder sb3 = new StringBuilder();
//            StringBuilder sb4 = new StringBuilder();
//            StringBuilder sb5 = new StringBuilder();
//            sb5.append("dataset,range,total_sites,no_silent_sites,no_replacement_sites\n");

            for(int dd=0; dd<no_datasets; dd++) {

                no_timepoints = timepoints_per_dataset[dd];
                timepoints = timepoints_multi.get(datasets[dd]);

                if(reader1.ready()) {

                    String ancesfilename = reader1.readLine().trim();
                    MainAlignmentParser ra = new MainAlignmentParser(ancesfilename);

                    int[][] ancestralMatrix = ra.readFASTA();
                    int[] ans = ra.consensusArray(ancestralMatrix);

//                    sb1.append(datasets[t]);
//                    sb2.append(datasets[t]);
//                    sb3.append(datasets[t]);
//                    sb4.append(datasets[t]);
//                    System.out.println(datasets[t]);
                    //no_of_codons[t] = Math.round(ans.length/3);

                    DescriptiveStatistics r_m = new DescriptiveStatistics();
                    DescriptiveStatistics s_m = new DescriptiveStatistics();
                    for(int tt=0; tt < no_timepoints; tt++) {

                        if(reader2.ready()) {

                        	String mainFile = reader2.readLine().trim();
                            MainAlignmentParser r = new MainAlignmentParser(mainFile);
                            int[][] main = r.readFASTA();
                            while((main == null || main.length < 10) && reader2.ready()) {

                                mainFile = reader2.readLine().trim();
                                r = new MainAlignmentParser(mainFile);
                                main = r.readFASTA();
                                tt+=1;
                                if(tt == no_timepoints) {
                                    break;
                                }
                            }

                            int[] ans_tmp = ans;

                            if(which.containsKey(datasets[dd])) {

                                //int[] list = codonlist.get(datasets[t]);
                                int [][] tmp_1 = TeaspoonMethods.subsetter(main, map, which.get(datasets[dd]));
                                int [] tmp_2 = TeaspoonMethods.subsetter(ans_tmp, map, which.get(datasets[dd]));
                                main = tmp_1;
                                ans_tmp = tmp_2;
                            }

                            Williamson3bin w = new Williamson3bin(main,ans_tmp);

                            if (fixedNR) {
                                w.williamson3bin_method(nr[dd], this.low, this.mid, this.high);
                            } else {
                                w.williamson3bin_method(this.low, this.mid, this.high);
                            }

//                            sb1.append(",").append(w.mid_R);
//                            sb2.append(",").append(w.mid_S);
//                            sb3.append(",").append(w.Nr);
//                            sb4.append(",").append(tt).append(",").append(w.Adapt).append("\n");//
//                            sb5.append(datasets[dd]).append(",low,").append(w.low_R + w.low_S).append(",").append(w.low_S).append(",").append(w.low_R).append("\n");
//                            sb5.append(datasets[dd]).append(",mid,").append(w.mid_R + w.mid_S).append(",").append(w.mid_S).append(",").append(w.mid_R).append("\n");
//                            sb5.append(datasets[dd]).append(",high,").append(w.high_R + w.high_S).append(",").append(w.high_S).append(",").append(w.high_R).append("\n");
//                            sb5.append("\n");

                            System.out.println(tt + ": " + mainFile + " : " + main[0].length + " : "+main.length+" neutralratio:" + w.Nr);

                            TeaspoonMethods.record(low, datasets[dd],new double[]{dd, Double.parseDouble(timepoints[tt]), 0}, w);
                            TeaspoonMethods.record(mid, datasets[dd],new double[]{dd, Double.parseDouble(timepoints[tt]), 1}, w);
                            TeaspoonMethods.record(high, datasets[dd],new double[]{dd, Double.parseDouble(timepoints[tt]), 2}, w);

                            if (!Double.isNaN(w.getMidR())) {
                                r_m.addValue(w.getMidR());
                            }

                            if (!Double.isNaN(w.getMidS())) {
                                s_m.addValue(w.getMidS());
                            }
                        }
                    }
                    
                    System.out.println(">" + datasets[dd] + ": r_m = " + r_m.getSum() + ", s_m = " + s_m.getSum() + " average_nr = " + r_m.getSum() / s_m.getSum());
                    high.append("\n");
                    low.append("\n");
                    mid.append("\n");
                }
//                sb1.append("\n");
//                sb2.append("\n");
//                sb3.append("\n");
//                sb4.append("\n");
            }

            lowfreq.write(low.toString());
            midfreq.write(mid.toString());
            highfreq.write(high.toString());

            reader1.close();
            reader2.close();
            highfreq.close();
            midfreq.close();
            lowfreq.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    @Override
    /*
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#w3bAnalysisBootstrap(int)
     */
    public void w3bAnalysisBootstrap(int bootstraps) {

        this.bootstraps = bootstraps;
        value_matrix = new TeaspoonValues[this.no_timepoints][this.no_datasets];

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
            ancestral = new MainAlignmentParser(reader1.readLine().trim(), false);
            ancestral.readFASTA();

            while (reader2.ready()) {
                String filename = reader2.readLine().trim();
                MainAlignmentParser m = new MainAlignmentParser(filename, false);
                System.out.println(filename);
                main_alignments.add(m.readFASTA());
            }
            
            reader1.close();
            reader2.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        assert ancestral != null;
        int[][] ancestralMatrix = ancestral.getSequenceMatrix();
        int[] ans = ancestral.consensusArray(ancestralMatrix);
        Map<String, DescriptiveStatistics> dataStatistics = new HashMap<>();

        for(String d: datasets) {
            dataStatistics.put(d, new DescriptiveStatistics());
        }

        for (int bs = 0; bs < bootstraps; bs++) {

            //int c = 0;
//            if (which.containsKey(datasets[0])) {
//                int[] tmp_2 = methods.Subsetter(ans, map, which.get(datasets[0]));
//                ans = tmp_2;
//
//            }
            //choosing your numCodons
            //int [] ans = ancestral.consensusArray(ancestralMatrix);

            RandomGenerator generator = new RandomGenerator();
            int[] sampler = new int[ans.length / 3];

            for (int x = 0; x < sampler.length; x++) {
                sampler[x] = generator.nextInt(sampler.length - 1);//generator.nextInt(sampler.length-1);
            }

            int c = 0;

            for (int t = 0; t < no_datasets; t++) {

                int [] ans_tmp = ans;
//                System.out.println(neutral_ratio[t]);

                no_timepoints = timepoints_per_dataset[t];
                timepoints = timepoints_multi.get(datasets[t]);

                double[] x_time_points = new double[no_timepoints+1];
                double[] y_adapt_points = new double[no_timepoints+1];
                //               int[] ans_tmp = ans;

                if (which.containsKey(datasets[t])) {

                    int[] tmp_2 = TeaspoonMethods.subsetter(ans_tmp, map, which.get(datasets[t]));

                    ans_tmp = tmp_2;
                    sampler = new int[ans_tmp.length / 3];

                    for (int x = 0; x < sampler.length; x++) {
                        sampler[x] = generator.nextInt(sampler.length - 1);//generator.nextInt(sampler.length-1);
                    }
                }

                //System.out.println(datasets[t]);
                SimpleRegression regression = new SimpleRegression(true);
                regression.addData(firstTimepoint[t], 0); //2013.33 should be first timepoint...
                //System.out.println(2013.33+","+0);
                for (int d = 0; d < no_timepoints; d++) {

                    int[][] main = main_alignments.get(c);
                    //if(reader2.ready()) {
                    //String mainFile = reader2.readLine().trim();

                    //MainAlignmentParser r = new MainAlignmentParser(mainFile);
                    //int [] ans_tmp = ancestral.consensusArray(ancestralMatrix);
                    //System.out.println(c);
                    x_time_points[0] = 2013.33;

                    //System.out.println(d+","+main);
                    while((main == null || main.length < 10)) {
                        d+=1;
                        c++;
                        if(d == no_timepoints && c== main_alignments.size()) {
                            break;
                        }
                        //System.out.println(main);
                        main = main_alignments.get(c);
                        //System.out.println(">"+main);
                        //System.out.println(c);
                    }
                    //System.out.println(c);
                    //System.out.println(d+","+main);

                    if (which.containsKey(datasets[t])) {

                        //int[] list = codonlist.get(datasets[t]);
                        int[][] tmp_1 = TeaspoonMethods.subsetter(main, map, which.get(datasets[t]));
                        //int[] tmp_2 = methods.Subsetter(ans_tmp, map, which.get(datasets[t]));
                        main = tmp_1;
                        //ans_tmp = tmp_2;
                    }

                    //int[] tmp_2 = ancestral.consensusArray(ancestralMatrix);

                    //System.out.println(main);
                    Williamson3bin w = new Williamson3bin(main, ans_tmp); //******
                    Store s = w.CreateBlocks(3, main[0].length, sampler); //******
                    Williamson3bin w3b = new Williamson3bin(s.RandomisedIntegerMatrix, s.RandomisedIntegerAncestral);

                    w3b.williamson3bin_method(nr[t], this.low, this.mid, this.high);

                    //System.out.println(timepoints[d]+","+datasets[t]);
                    value_matrix[d][t].setRow(timepoints[d]);
                    value_matrix[d][t].setColumn(datasets[t]);
                    value_matrix[d][t].setNumCodons(sampler.length);

                    value_matrix[d][t].getReplacementPolymorphsMidClass()[bs] = w3b.getMidR();
                    value_matrix[d][t].getSilentPolymorphsMidClass()[bs] = w3b.getMidS();
                    value_matrix[d][t].getReplacementSubsHigh()[bs] = w3b.getHighR();
                    value_matrix[d][t].getSilentSubsHigh()[bs] = w3b.getHighS();
                    value_matrix[d][t].getSilentPolymorphsLowClass()[bs] = w3b.getLowS();
                    value_matrix[d][t].getReplacementPolymorphsLowClass()[bs] = w3b.getLowR();

                    if (w3b.getHighR() > 0) {
                        value_matrix[d][t].getSilentToReplacementPolymorphsHighClassRatio()[bs] = w3b.getHighS() / w3b.getHighR();
                    }

                    if (w3b.getHighS() > 0) {
                        value_matrix[d][t].getReplacementToSilentPolymorphsHighClassRatio()[bs] = w3b.getHighR() / w3b.getHighS();
                    }

                    if (w3b.getLowS() > 0) {
                        value_matrix[d][t].getReplacementToSilentPolymorphsLowClassRatio()[bs] = w3b.getLowR()  / w3b.getLowS() ;

                    }
                    
                    if (w3b.getMidS()  > 0) {
                        value_matrix[d][t].getNeutralRatio()[bs] = w3b.Nr;
                    }

                    if (Double.isNaN(w3b.getAdaptiveMutations())) {
                        w3b.setAdaptiveMutations(0.0);
                    }
                    
                    value_matrix[d][t].getAdaptations()[bs] = w3b.getAdaptiveMutations();
                    DiversityStats diversityStats = new DiversityStats(w3b.integer_matrix);
                    double[] wattersonEstimates = diversityStats.wattersonEstimates();
                    value_matrix[d][t].getTajimasD()[bs] = diversityStats.TajimasD();
                    
                    if (Double.isNaN(value_matrix[d][t].getTajimasD()[bs])) {
                        value_matrix[d][t].getTajimasD()[bs] = 0.0;
                    }
                    
                    value_matrix[d][t].getWattersonS()[bs] = wattersonEstimates[0];
                    value_matrix[d][t].getWattersonPi()[bs] = wattersonEstimates[1];
                    value_matrix[d][t].getTheta()[bs] = diversityStats.theta();

                    c++;

                    x_time_points[d+1]=Double.parseDouble(timepoints[d]);
                    y_adapt_points[d+1]=w3b.getAdaptiveMutations();

                    if(x_time_points[d+1]>0) {
                        regression.addData(Double.parseDouble(timepoints[d]), y_adapt_points[d+1]);
                    }
                }

                double slope = regression.getSlope();
                dataStatistics.get(datasets[t]).addValue(slope);
            }
            System.out.println("I am on Run  " + (bs + 1) + "  of " + bootstraps);
        }

        System.out.println();
        System.out.println("************Adaptive substitutions per year************\n");

        for(String d: datasets) {
            System.out.println(d+", "+dataStatistics.get(d).getMean()+", "+dataStatistics.get(d).getPercentile(25)+", "+dataStatistics.get(d).getPercentile(75));
        }
        
//        System.out.println(+regression_PRD.getMean()+", "+regression_PRD.getPercentile(25)+", "+regression_PRD.getPercentile(75));
//        System.out.println("YRD1 "+regression_YRD1.getMean()+", "+regression_YRD1.getPercentile(25)+", "+regression_YRD1.getPercentile(75));
//        System.out.println("YRD2.2 "+regression_YRD22.getMean()+", "+regression_YRD22.getPercentile(25)+", "+regression_YRD22.getPercentile(75));
//        System.out.println("YRD2.1 "+regression_YRD21.getMean()+", "+regression_YRD21.getPercentile(25)+", "+regression_YRD21.getPercentile(75));

        System.out.println();

        String output = mainFile;
        
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(output.replace(".txt", "_bootstraps_n_"+bootstraps+"_adaptations.csv")));

            writer.write("dataset, Time,Mean,Median,Lower Quartile,Upper Quartile,Standard Deviation,No of Codons\n");

            for (int d = 0; d < no_datasets; d++) {
                String s = String.valueOf(datasets[d]);

                writer.write(s+","+firstTimepoint[d] + "," + 0 + ","+ 0 + "," + 0 + "," + 0 + "," + 0 + "," + 0+"\n");

                for (int t = 0; t < no_timepoints; t++) {

                    DescriptiveStatistics bstraps = new DescriptiveStatistics(value_matrix[t][d].getAdaptations());

                    double mean = bstraps.getMean();
                    double lq = bstraps.getPercentile(25);
                    double uq = bstraps.getPercentile(75);
                    double med = bstraps.getPercentile(50);
                    double std = bstraps.getStandardDeviation();

                    if(value_matrix[t][d].getNumCodons()>0) {
                        System.out.println(s+","+timepoints[t] + "," + mean + ","+ med + "," + lq + "," + uq + "," + std + "," + value_matrix[t][d].getNumCodons());
                        writer.write(s + "," + timepoints[t] + "," + mean + "," + med + "," + lq + "," + uq + "," + std + "," + value_matrix[t][d].getNumCodons() + "\n");
                    }
//                System.out.println();
//                System.out.println(d + "," + t + "," + value_matrix[t][d].r_high[d]);
//                for (int b = 0; b < 10; b++) {
//                    s += "," + value_matrix[t][d].adaptations[b];
//                }

                }
                writer.write("\n");
                System.out.println();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /*
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#getDiversityStats()
     */
    public void getDiversityStats() {

        List<int[][]> main_alignments = new ArrayList<int[][]>();    //list of main_alignments that represent different timepoints...
        MainAlignmentParser ancestral;

        //no_timepoints = timepoints.length+1;

        try {
            BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
            ancestral = new MainAlignmentParser(reader1.readLine().trim(),true);
            main_alignments.add(ancestral.readFASTA());
            BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));
            String output = mainFile;

            while(reader2.ready()) {
                String filename = reader2.readLine().trim();
                MainAlignmentParser m = new MainAlignmentParser(filename, true);
                //System.out.println(filename);
                main_alignments.add(m.readFASTA());
            }

            BufferedWriter summaryResults = new BufferedWriter(new FileWriter(output.replace(".txt", "_diversityStatsTable.csv")));
            StringBuilder summary = new StringBuilder();

            summary.append("Dataset,Timepoint,Average Pairwise Diversity (pi),No of segregating sites (S),Tajima's D,no_of_seqs\n");
            //int no_datasets = datasets.length;

            int count = 0;
            
            for(int d=0; d<no_datasets; d++) {

                no_timepoints = timepoints_per_dataset[d];

                timepoints = timepoints_multi.get(datasets[d]);

                System.out.println(datasets[d]);

                for (int t = 0; t < no_timepoints; t++) {

                    double pi = Double.POSITIVE_INFINITY;
                    double S = Double.POSITIVE_INFINITY;
                    double tajimasD = Double.POSITIVE_INFINITY;
                    int[][] site_main = main_alignments.get(count);

                    if(byGene) {
                        site_main = TeaspoonMethods.subMatrix(site_main, gene_boundary[0], gene_boundary[1], false);
                    }
                    
                    count++;

                    if (site_main.length > 100) {
                        DiversityStats stats = new DiversityStats(site_main);
                        pi = stats.numberOfPairwiseDifferences();
                        S = stats.numberOfSegregatingSites();
                        tajimasD = stats.TajimasD();
                    }
                    summary.append(datasets[d]).append(",").append(timepoints[t]).append(",").append(pi).append(",").append(S).append(",").append(tajimasD).append(",").append(site_main.length).append("\n");
                }
            }
            System.out.println(summary);
            summaryResults.write(summary.toString());
            summaryResults.write("\n");
            
            reader1.close();
            reader2.close();
            summaryResults.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    /*
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#getValue_matrix()
     */
    public TeaspoonValues[][] getValue_matrix() {
        return this.value_matrix;
    }

    @Override
    /*
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#getDatasets()
     */
    public String[] getDatasets() {
        return this.datasets;
    }

    @Override
    /*
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#getTimepoints_per_dataset()
     */
    public int[] getTimepoints_per_dataset() {
        return this.timepoints_per_dataset;
    }

    /*
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#getBootstraps()
     */
    public int getBootstraps() {
        return this.bootstraps;
    }

    /*
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#setCodonMap(int[])
     */
    public void setCodonMap(int[] map){
        this.map = map;
    }

    @Override
    /*
     * (non-Javadoc)
     * @see teaspoon.adaptation.Analysis#setWhichMap(java.util.Map)
     */
    public void setWhichMap(Map<String, Integer> which) {
        this.which = which;
    }


//    public void getBootstrapsByW3Bin() {
//        value_matrix = new TeaspoonValues[this.no_timepoints][this.no_datasets];
//        for(int i=0; i<value_matrix.length; i++){
//
//            for(int j = 0; j < value_matrix[0].length; j++) {
//
//                value_matrix[i][j] = new TeaspoonValues(bootstraps);
//
//            }
//        }
//
//        MainAlignmentParser ancestral = null;
//        List<int[][]> main_alignments = new ArrayList<int[][]>();
//
//        try {
//            BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
//            ancestral = new MainAlignmentParser(reader1.readLine().trim(), true);
//            ancestral.readFASTA();
//
//            BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));
//
////                }
//            while (reader2.ready()) {
//
//                String filename = reader2.readLine().trim();
//                MainAlignmentParser m = new MainAlignmentParser(filename, true);
//                System.out.println(filename);
//                main_alignments.add(m.readFASTA());
//            }
//        } catch (FileNotFoundException e) {
//            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//        } catch (IOException e) {
//            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//        }
//
//        int n = 0;
//
//        for( int bs=0; bs < bootstraps; bs++) {
////            try {
////                BufferedReader reader1 = new BufferedReader(new FileReader("ancestral_filelist_nr.txt"));
////                BufferedReader reader2 = new BufferedReader(new FileReader("main_filelist_nr.txt"));
//
//            for(int t=0; t<no_datasets; t++) {
//
//                System.out.println(t);
//
////                    if(reader1.ready()) {
//
////                        String ancesfilename = reader1.readLine().trim();
//
//
//                //MainAlignmentParser ra = new MainAlignmentParser(ancesfilename);
//                int[][] ancestralMatrix = ancestral.sequenceMatrix;
//                int[] ans = ancestral.consensusArray(ancestralMatrix);
//                int[] sampler = new int[ans.length/3];
//
//                //choosing your numCodons
//                RandomGenerator generator = new RandomGenerator();
//                for(int x=0; x< sampler.length; x++) {
//                    sampler[x] = generator.nextInt(sampler.length-1);
//                }
//
//                //no_of_codons[t] = Math.round(ans.length/3);
//
//                for(int d=0; d < no_timepoints; d++) {
//
//                    //if(reader2.ready()) {
//                    //String mainFile = reader2.readLine().trim();
//                    n++;
//                    //MainAlignmentParser r = new MainAlignmentParser(mainFile);
//                    int[][] main = main_alignments.get(d);
//
//                    int[] ans_tmp = ans;
//                    if(which.containsKey(datasets[t])) {
//
//                        //int[] list = codonlist.get(datasets[t]);
//                        int [][] tmp_1 = methods.Subsetter(main, map, which.get(datasets[t]));
//                        int [] tmp_2 = methods.Subsetter(ans_tmp, map, which.get(datasets[t]));
//
//                        main = tmp_1;
//                        ans_tmp = tmp_2;
//
//                        sampler = new int[ans_tmp.length/3];
//                        for(int x=0; x< sampler.length; x++) {
//                            sampler[x] = generator.nextInt(sampler.length-1);
//                        }
//                    }
//
//                    Williamson3bin ww = new Williamson3bin(main,ans_tmp);
//                    //creating matrices with randomly chosen numCodons
//                    Store s = ww.CreateBlocks(3,main[0].length,sampler);
//                    Williamson3bin w = new Williamson3bin(s.RandomisedIntegerMatrix,s.RandomisedIntegerAncestral);
//                    w.williamson3bin_method(low, mid, high);
//
//                    value_matrix[d][t].row = timepoints[d];
//                    //value_matrix[d][t].column = datasets[t];
//                    value_matrix[d][t].dataset = datasets[t];
//                    value_matrix[d][t].codons = sampler.length;
//
//                    value_matrix[d][t].r_mid[bs] = w.mid_R;
//                    value_matrix[d][t].s_mid[bs] = w.mid_S;
//                    value_matrix[d][t].adaptations[bs] = w.Adapt;
//                    value_matrix[d][t].neutral_ratio[bs] = w.Nr;
//                }
//            }
//            System.out.println("I am on Run  "+(bs+1)+"  of "+bootstraps);
//        }
//
//
//
////            }
////            catch (FileNotFoundException e) {
////                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
////            } catch (IOException e) {
////                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
////            }
//
//
//
//    }
//
//    // add an option whether to fix neutral_ratio or estimate
//    public void getBootstrapsByBm(int bootstraps) {
//
//        this.bootstraps = bootstraps;
//        value_matrix = new TeaspoonValues[this.no_timepoints][this.no_datasets];
//
//
//        for (int i = 0; i < value_matrix.length; i++) {
//
//            for (int j = 0; j < value_matrix[0].length; j++) {
//
//                value_matrix[i][j] = new TeaspoonValues(bootstraps);
//
//            }
//        }
//
//
//        MainAlignmentParser ancestral = null;
//        List<int[][]> main_alignments = new ArrayList<int[][]>();
//
//        try {
//            BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
//            ancestral = new MainAlignmentParser(reader1.readLine().trim(), true);
//            ancestral.readFASTA();
//
//            BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));
//
////                }
//            while (reader2.ready()) {
//
//                String filename = reader2.readLine().trim();
//                MainAlignmentParser m = new MainAlignmentParser(filename, true);
//                System.out.println(filename);
//                main_alignments.add(m.readFASTA());
//            }
//        } catch (IOException e) {
//            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//        }
//
//
//        for (int bs = 0; bs < bootstraps; bs++) {
////            try {
////                BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
////                BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));
//
//            for (int t = 0; t < no_datasets; t++) {
//
//                //System.out.println(t);
//
//                //if(reader1.ready()) {
//
//                //String ancesfilename = reader1.readLine().trim();
//
//
//                //MainAlignmentParser ra = new MainAlignmentParser(ancesfilename);
//                assert ancestral != null;
//                int[][] ancestralMatrix = ancestral.sequenceMatrix;
//                int[] ans = ancestral.consensusArray(ancestralMatrix);
//                int[] sampler = new int[ans.length / 3];
//
//                //choosing your numCodons
//                RandomGenerator generator = new RandomGenerator();
//                for (int x = 0; x < sampler.length; x++) {
//                    sampler[x] = generator.nextInt(sampler.length - 1);//generator.nextInt(sampler.length-1);
//                }
//
//                //no_of_codons[t] = Math.round(ans.length/3);
//
//                no_timepoints = timepoints_per_dataset[t];
//                timepoints = timepoints_multi.get(datasets[t]);
//                //System.out.println(timepoints.length);
//
//
//                for (int d = 0; d < no_timepoints; d++) {
//
//                    //if(reader2.ready()) {
//                    //String mainFile = reader2.readLine().trim();
//
//                    //MainAlignmentParser r = new MainAlignmentParser(mainFile);
//                    int[][] main = main_alignments.get(d);
//
//                    int[] ans_tmp = ans;
//                    if (which.containsKey(datasets[t])) {
//
//                        //int[] list = codonlist.get(datasets[t]);
//                        int[][] tmp_1 = methods.Subsetter(main, map, which.get(datasets[t]));
//                        int[] tmp_2 = methods.Subsetter(ans_tmp, map, which.get(datasets[t]));
//
//                        main = tmp_1;
//                        ans_tmp = tmp_2;
//
////                                    sampler = new int[ans_tmp.length/3];
////                                    for(int x=0; x< sampler.length; x++) {
////                                        sampler[x] = generator.nextInt(sampler.length-1);
////                                    }
//
//                    }
//
//                    BhattMethod b = new BhattMethod(main, ans_tmp); //******
//                    Store s = b.CreateBlocks(3, main[0].length, sampler); //******
//                    BhattMethod bm = new BhattMethod(s.RandomisedIntegerMatrix, s.RandomisedIntegerAncestral);
//
//                    bm.Method(bins, prior, true, Nvec, neutral_ratio[t]);
//
//
//                    //System.out.println(timepoints[d]+","+datasets[t]);
//                    value_matrix[d][t].row = timepoints[d];
//                    value_matrix[d][t].column = datasets[t];
//                    value_matrix[d][t].codons = sampler.length;
//
//                    value_matrix[d][t].r_mid[bs] = bm.ReplacementCountArray[1];
//                    value_matrix[d][t].s_mid[bs] = bm.SilentCountArray[1];
//                    value_matrix[d][t].r_high[bs] = bm.ReplacementCountArray[2];
//                    value_matrix[d][t].s_high[bs] = bm.SilentCountArray[2];
//                    value_matrix[d][t].s_low[bs] = bm.SilentCountArray[0];
//                    value_matrix[d][t].r_low[bs] = bm.ReplacementCountArray[0];
//
//
//                    if (bm.ReplacementCountArray[2] > 0) {
//                        value_matrix[d][t].sr_ratio_high[bs] = bm.SilentCountArray[2] / bm.ReplacementCountArray[2];
//                    }
//
//                    if (bm.SilentCountArray[2] > 0) {
//                        value_matrix[d][t].rs_ratio_high[bs] = bm.ReplacementCountArray[2] / bm.SilentCountArray[2];
//                    }
//
//                    if (bm.SilentCountArray[0] > 0) {
//                        value_matrix[d][t].rs_low[bs] = bm.ReplacementCountArray[0] / bm.SilentCountArray[0];
//
//                    }
//                    if (bm.SilentCountArray[1] > 0) {
//                        value_matrix[d][t].neutral_ratio[bs] = bm.neutralratio;
//
//                    }
//
//                    if (Double.isNaN(bm.Adaptation)) {
//                        bm.Adaptation = 0.0;
//                    }
//                    value_matrix[d][t].adaptations[bs] = bm.Adaptation;
//
//                    DiversityStats diversityStats = new DiversityStats(bm.integer_matrix);
//                    double[] wattersonEstimates = diversityStats.wattersonEstimates();
//
//                    value_matrix[d][t].tajimas_D[bs] = diversityStats.TajimasD();
//                    if (Double.isNaN(value_matrix[d][t].tajimas_D[bs])) {
//                        value_matrix[d][t].tajimas_D[bs] = 0.0;
//                    }
//                    value_matrix[d][t].watterson_S[bs] = wattersonEstimates[0];
//                    value_matrix[d][t].watterson_pi[bs] = wattersonEstimates[1];
//                    value_matrix[d][t].theta[bs] = diversityStats.theta();
//
//
//                    //}
//                    //}
//                    //}
//                    //}
//
//                }
//            }
//            System.out.println("I am on Run  " + (bs + 1) + "  of " + bootstraps);
//
//
//        }
//    }
//
//    //}
//
//    //}
//
//    public void getBootstrapByW3b(int bootstraps) {
//
//        this.bootstraps = bootstraps;
//        value_matrix = new TeaspoonValues[this.no_timepoints][this.no_datasets];
//        for(int i=0; i<value_matrix.length; i++){
//
//            for(int j = 0; j < value_matrix[0].length; j++) {
//
//                value_matrix[i][j] = new TeaspoonValues(bootstraps);
//
//            }
//        }
//
//
//        for( int bs=0; bs < bootstraps; bs++) {
//            try {
//                BufferedReader reader1 = new BufferedReader(new FileReader(ancestralFile));
//                BufferedReader reader2 = new BufferedReader(new FileReader(mainFile));
//
//
//                for(int t=0; t<no_datasets; t++) {
//
//
//                    if(reader1.ready()) {
//
//                        String ancesfilename = reader1.readLine().trim();
//                        MainAlignmentParser ra = new MainAlignmentParser(ancesfilename);
//                        int[][] tmp = ra.readNEXUS();
//                        int[] ans = ra.consensusArray(tmp);
//                        int[] sampler = new int[ans.length/3];
//
//                        //choosing your numCodons
//                        RandomGenerator generator = new RandomGenerator();
//                        for(int x=0; x< sampler.length; x++) {
//                            sampler[x] = generator.nextInt(sampler.length-1);//generator.nextInt(sampler.length-1);
//                        }
//
//                        no_timepoints = timepoints_per_dataset[t];
//                        timepoints = timepoints_multi.get(datasets[t]);
//
//                        for(int d=0; d < no_timepoints; d++) {
//
//                            if(reader2.ready()) {
//                                String mainFile = reader2.readLine().trim();
//                                MainAlignmentParser r = new MainAlignmentParser(mainFile);
//                                int[][] main = r.readNEXUS();
//
//                                Williamson3bin w = new Williamson3bin(main,ans);
//                                Store s = w.CreateBlocks(3, main[0].length, sampler);
//                                Williamson3bin w3b = new Williamson3bin(s.RandomisedIntegerMatrix, s.RandomisedIntegerAncestral);
//
//                                w3b.williamson3bin_method(neutral_ratio[t], low, mid, high);
//                                value_matrix[d][t].row = timepoints[d];
//                                value_matrix[d][t].column = datasets[t];
//                                value_matrix[d][t].adaptations[bs] = w3b.Adapt;
//                                value_matrix[d][t].row = timepoints[d];
//                                value_matrix[d][t].column = datasets[t];
//                                value_matrix[d][t].codons = sampler.length;
//
//                                value_matrix[d][t].r_mid[bs] = w3b.mid_R;
//                                value_matrix[d][t].s_mid[bs] = w3b.mid_S;
//                                value_matrix[d][t].r_high[bs] = w3b.high_R;
//                                value_matrix[d][t].s_high[bs] = w3b.high_R;
//                                value_matrix[d][t].s_low[bs] = w3b.low_S;
//                                value_matrix[d][t].r_low[bs] = w3b.low_R;
//                                if(w3b.high_R>0) {
//                                    value_matrix[d][t].sr_ratio_high[bs] = w3b.high_S/w3b.high_R;
//                                }
//
//                                if(w3b.high_S>0) {
//                                    value_matrix[d][t].rs_ratio_high[bs] = w3b.high_R/w3b.high_S;
//                                }
//
//                                if(w3b.low_S>0) {
//                                    value_matrix[d][t].rs_low[bs] = w3b.low_R/w3b.low_S;
//
//                                }
//
//                                if(Double.isNaN(w3b.Adapt)){
//                                    w3b.Adapt = 0.0;
//                                }
//                                value_matrix[d][t].adaptations[bs] = w3b.Adapt;
//                                value_matrix[d][t].neutral_ratio[bs] = w3b.Nr;
//
//                                DiversityStats diversityStats = new DiversityStats(w3b.integer_matrix);
//                                double[] wattersonEstimates = diversityStats.wattersonEstimates();
//
//                                value_matrix[d][t].tajimas_D[bs] = diversityStats.TajimasD();
//                                if(Double.isNaN(value_matrix[d][t].tajimas_D[bs])) {
//                                    value_matrix[d][t].tajimas_D[bs]  = 0.0;
//                                }
//                                value_matrix[d][t].watterson_S[bs] = wattersonEstimates[0];
//                                value_matrix[d][t].watterson_pi[bs] = wattersonEstimates[1];
//                                value_matrix[d][t].theta[bs] = diversityStats.theta();
//
//                            }
//                        }
//                    }
//                }
//
//            }
//            catch (FileNotFoundException e) {
//                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            } catch (IOException e) {
//                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
//
//        }
//
//    }
}
