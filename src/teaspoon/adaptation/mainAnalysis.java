package teaspoon.adaptation;


import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.io.BufferedWriter;
import java.util.HashMap;

/**
 * Created with IntelliJ IDEA.
 * User: jayna
 * Date: 13/08/2013
 * Time: 13:58
 * To change this template use File | Settings | File Templates.
 */
public class mainAnalysis {

    String ancestralfilename;
    String mainfilename;
    String [] datasets;
    int[][] gene_coordinates;
    boolean fixedNR;
    double[] nr;
    String [] timepoints;


    public mainAnalysis() {

        this.fixedNR = false;
        this.nr = null;
//        this.ancestralfilename = ancestralFile;
//        this.mainfilename = mainFile;
//
//        timepoints_multi = new HashMap<String, String[]> ();
//
//        nr = new double[no_datasets];
//
//        bins = new double[2][L.length];
//        for(int i=0;i<L.length;i++){
//            bins[0][i]=L[i];
//            bins[1][i]=H[i];
//        }


    }

    public void runBM_oneTimepoint() {

        analyseGene analysis = new analyseGene(ancestralfilename, mainfilename);

        analysis.fixedNR = false;
        analysis.nr = new double[] {4.813952552};

        analysis.datasets = new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27"};
        analysis.no_datasets = analysis.datasets.length;
        analysis.timepoints_multi = new HashMap<String, String[]>();
        analysis.timepoints = new String[]{"1"};

        //write method to set timepoints_multi, timepoints, timepoints_per_dataset, no_timepoints

        for(int i=0; i<analysis.datasets.length; i++) {
            analysis.timepoints_multi.put(analysis.datasets[i], analysis.timepoints);
        }

        analysis.timepoints_per_dataset = new int[analysis.no_datasets];
        int max_timepoints = 0;

        for(int i=0; i < analysis.no_datasets; i++) {

            analysis.timepoints_per_dataset[i] = analysis.timepoints_multi.get(analysis.datasets[i]).length;
            if(analysis.timepoints_per_dataset[i]>max_timepoints) {
                max_timepoints = analysis.timepoints_multi.get(analysis.datasets[i]).length;
            }
        }


        analysis.no_timepoints = max_timepoints;
        analysis.bmAnalysis();
    }

    public void runBM_multipleTimepoints() {

        analyseGene analysis = new analyseGene(ancestralfilename, mainfilename);

        analysis.fixedNR = this.fixedNR;
        analysis.nr = this.nr;
//        analysis.nr = new double[] {0.018};
//        analysis.nr = new double[]{0.22909292816799146, 0.2624104494960995,0.14893244343173687,0.017836580830472952};
                //nr = new double[]{0.2378342136223165, 0.25807525692320593, 0.14893244343173687}; exclude YRD2.1
                //nr = new double[]{0.22909292816799146, 0.2624104494960995,0.06757118528744181};
//
//        //fluB
//        analysis.datasets = new String[] {"HA_yam_pink", "HA_yam_purple", "HA_vic_red", "HA_vic_blue"};
//
//        analysis.timepoints_multi.put("HA_yam_pink", new String [] {"1990","1992","1994","1996","1998","2000","2002","2004","2006","2008","2010","2012","2014"});
//        analysis.timepoints_multi.put("HA_yam_purple", new String [] {"1990","1992","1994","1996","1998","2000","2002","2004","2006","2008","2010","2012","2014"});
//        analysis.timepoints_multi.put("HA_vic_red", new String [] {"1990","1992","1994","1996","1998","2000","2002","2004","2006","2008","2010","2012","2014"});
//        analysis.timepoints_multi.put("HA_vic_blue", new String [] {"1990","1992","1994","1996","1998","2000","2002","2004","2006","2008","2010","2012","2014"});
//
//        analysis.nr = new double[] {0.4086278689480355, 0.4006102116532152, 0.3332121714520498, 0.27419848421258597}; //full HA m-ratio
//        //analysis.nr = new double[] {0.1596024989480543, 0.19231470087595393, 0.09949753779601428, 0.1625399773384175};  //int HA m-ratio
//        analysis.which.put("HA_yam_pink",1);
//        analysis.which.put("HA_yam_purple",1);
//        analysis.which.put("HA_vic_red",1);
//        analysis.which.put("HA_vic_blue",1);
//
        StructuralMap structuralMap = new StructuralMap();
//
//        analysis.map = structuralMap.FluA_H7;
//
//        analysis.which.put("PRD1",0);
//        analysis.which.put("YRD1",0);
//        analysis.which.put("YRD2.2",0);
//        analysis.which.put("YRD2.1",0);



        //HCV_pacbio

        //analysis.datasets = new String[] {"YRD2.1"};
//        analysis.datasets = new String[] {"PRD1", "YRD1", "YRD2.2", "YRD2.1"};
//        analysis.timepoints_multi.put("PRD1", new String[] {"1", "2", "3", "4"});
//        analysis.timepoints_multi.put("YRD1", new String[] {"1", "2"});
//        analysis.timepoints_multi.put("YRD2.2", new String[] {"1", "2", "3", "4"});
//        analysis.timepoints_multi.put("YRD2.1", new String[] {"1", "2", "3", "4"});


        analysis.datasets = new String[] {"PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"};
        analysis.timepoints_multi.put("PB2", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
        analysis.timepoints_multi.put("PB1", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
        analysis.timepoints_multi.put("PA", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
        analysis.timepoints_multi.put("HA", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
        analysis.timepoints_multi.put("NP", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
        analysis.timepoints_multi.put("NA", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
        analysis.timepoints_multi.put("MP", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
        analysis.timepoints_multi.put("NS", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});

//        analysis.timepoints_multi.put("p_4", new String[] {"2002.4274", "2003.0192", "2005.2055", "2006.0932", "2006.4521", "2007.1233", "2008.1366", "2009.1973", "2010.2849", "2012.4891", "2013.4521"});
//        //timepoints_multi.put("p_37", new String[] {"2002.377", "2002.6055", "2002.874", "2003.2", "2004.2077", "2005.4027", "2006.3178", "2007.3507"});
//        analysis.timepoints_multi.put("p_37", new String[] {"2002.0877", "2002.3370", "2002.6055", "2002.874", "2003.2", "2004.2077", "2005.4027", "2006.3178", "2007.3507"});
//        analysis.timepoints_multi.put("p_53", new String[] {"2001.5041", "2001.8685", "2002.1178", "2002.3479", "2002.874", "2003.2877", "2005.4027", "2008.347", "2009.1890"});
//        analysis.timepoints_multi.put("p_61", new String[] {"2002.0192", "2002.3863", "2002.6329", "2003.0986", "2003.7644", "2004.6667", "2005.011", "2006.0164", "2007.0274", "2008.2923", "2008.9481", "2009.074", "2009.2466", "2009.5425", "2010.3753", "2012.1858", "2012.8115", "2013.6438"});


        analysis.no_datasets = analysis.datasets.length;
        analysis.timepoints_per_dataset = new int[analysis.no_datasets];


        int max_timepoints = 0;
        for(int i=0; i < analysis.no_datasets; i++) {

            analysis.timepoints_per_dataset[i] = analysis.timepoints_multi.get(analysis.datasets[i]).length;
            if(analysis.timepoints_per_dataset[i]>max_timepoints) {
                max_timepoints = analysis.timepoints_multi.get(analysis.datasets[i]).length;
            }
            //System.out.println(bs.timepoints_multi.get(bs.datasets[i]).length);
        }


        analysis.no_timepoints = max_timepoints;
        analysis.bmAnalysis();

//       analysis.bmAnalysisBootstrap(100);

//        Value[][] value_matrix = analysis.getValue_matrix();
//        //for(int b=0; b < 100; b++) {
//
//        for (int d = 0; d < analysis.no_datasets; d++) {
//            String s = String.valueOf(analysis.datasets[d]);
//
//
//
//            for (int t = 0; t < analysis.no_timepoints; t++) {
//
//
//                DescriptiveStatistics bstraps = new DescriptiveStatistics(value_matrix[t][d].adaptations);
//
//                double lq = bstraps.getPercentile(25);
//                double uq = bstraps.getPercentile(75);
//                double med = bstraps.getPercentile(50);
//                double std = bstraps.getStandardDeviation();
//
//
//                System.out.println(s+","+t + "," + med + "," + lq + "," + uq + "," + std + "," + value_matrix[t][d].codons);
////                System.out.println();
////                System.out.println(d + "," + t + "," + value_matrix[t][d].rh[d]);
//
////                for (int b = 0; b < 10; b++) {
////
////                    s += "," + value_matrix[t][d].adaptations[b];
////
////                }
//
//            }
//            System.out.println();
//        }




    }

    public void runDeepGenomeAnalysis() {

        analyseDeepGenome analysis = new analyseDeepGenome(ancestralfilename, mainfilename);

        analysis.datasets = new String[] {"gag", "pol", "env", "nef"};

        analysis.bmAnalysis();


    }

    public void runBootstrapDeepGenomeAnalysis(int bootstraps) {

        analyseDeepGenome analysis = new analyseDeepGenome(ancestralfilename, mainfilename);

        analysis.datasets = new String[] {"gag", "pol", "env", "nef"};

        analysis.no_datasets = analysis.datasets.length;

        analysis.timepoints_per_dataset = new int[analysis.no_datasets];
        for(int i=0; i < analysis.no_datasets; i++) {

            analysis.timepoints_per_dataset[i] = analysis.timepoints.length;

        }

        analysis.fixedNR = true;

    }

    public void readParams(String filename) {

        this.ancestralfilename = filename;
        this.mainfilename = filename;
    }

    public static void main(String [] args) {


        mainAnalysis mainAnalysis = new mainAnalysis();

//        mainAnalysis.ancestralfilename = "/Users/jayna/Dropbox/H7N9_HA_Jing/ancestralfile.txt";
//        mainAnalysis.mainfilename = "/Users/jayna/Dropbox/H7N9_HA_Jing/mainfile4.txt";




//        adaptarate_H3N8
        mainAnalysis.ancestralfilename = args[0];
        mainAnalysis.mainfilename = args[1];



        if(args.length>2) {
            mainAnalysis.fixedNR = Boolean.parseBoolean(args[2]);
            mainAnalysis.nr = new double[] { Double.parseDouble(args[3]), Double.parseDouble(args[4]), Double.parseDouble(args[5]), Double.parseDouble(args[6]),
                    Double.parseDouble(args[7]), Double.parseDouble(args[8]), Double.parseDouble(args[9]), Double.parseDouble(args[10])};
        }




        mainAnalysis.runBM_multipleTimepoints();

    }



}