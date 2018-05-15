package teaspoon.adaptation;


import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: jayna
 * Date: 13/08/2013
 * Time: 13:58
 * To change this template use File | Settings | File Templates.
 * @author @jnarag
 * @version 0.0.1
 */
public class mainAnalysis {

    String ancestralfilename;
    String mainfilename;
    String [] datasets;
    int[][] gene_coordinates;
    boolean useFixedNeutralRatio;
    double[] neutralRatio;
    String [] timepoints;
    Map<String, String[]> timepoints_multi;
    int[][] genes;
    double[] firstTimepoint;


    /**
     * Default no-arg constructor.
     */
    public mainAnalysis() {

        this.useFixedNeutralRatio = false;
        this.neutralRatio = null;
        this.timepoints_multi = new HashMap<>();
        /*
         * This was commented out by Jayna previous to v0.1
        this.ancestralfilename = ancestralFile;
        this.mainfilename = mainFile;

<<<<<<< HEAD
        timepoints_multi = new HashMap<String, String[]> ();
=======
//        this.ancestralfilename = ancestralFile;
//        this.mainfilename = mainFile;
//
//        timepoints_multi = new HashMap<String, String[]> ();
//
//        neutral_ratio = new double[no_datasets];
//
//        bins = new double[2][L.length];
//        for(int i=0;i<L.length;i++){
//            bins[0][i]=L[i];
//            bins[1][i]=H[i];
//        }
>>>>>>> 3556e679f42cf8e89b9e1c373c236414aad689b7

        neutralRatio = new double[no_datasets];

        bins = new double[2][L.length];
        for(int i=0;i<L.length;i++){
            bins[0][i]=L[i];
            bins[1][i]=H[i];
        }
         */
    }

    /**
     * The simplest type of analysis.
     * @see teaspoon.adaptation.analyseGene
     * @see teaspoon.adaptation.mainAnalysis#runBootstrapOneTimepoint
     */
    public void runOneTimepoint() {

        analyseGene analysis = new analyseGene(ancestralfilename, mainfilename);

        analysis.fixedNR = this.useFixedNeutralRatio;
        analysis.nr = this.neutralRatio;

        analysis.datasets = this.datasets;
        analysis.no_datasets = analysis.datasets.length;
        analysis.timepoints_multi = new HashMap<String, String[]>();
        analysis.timepoints = this.timepoints;

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

    /**
     * A single-timepoint analysis with bootstrap replicates. 
     * <br/>Assumes parameters are set elsewhere in the object, e.g. 
     * <pre>
     * ancestralfilename, mainfilename, useFixedNeutralRatio, neutralRatio, datasets, timepoints
     * </pre>
     * @param bootstraps - number of bootstraps (? &ge;1 assumed?)
     * @return analysis an analyseGene object containing the results
     * @see teaspoon.adaptation.analyseGene
     * @see teaspoon.adaptation.mainAnalysis#runOneTimepoint
     */
    public analyseGene runBootstrapOneTimepoint(int bootstraps) {

        analyseGene analysis = new analyseGene(ancestralfilename, mainfilename);

        analysis.fixedNR = this.useFixedNeutralRatio;
        analysis.nr = this.neutralRatio;

        analysis.datasets = this.datasets;
        analysis.no_datasets = analysis.datasets.length;
        analysis.timepoints_multi = new HashMap<String, String[]>();
        analysis.timepoints = this.timepoints;

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
        analysis.bmAnalysisBootstrap(bootstraps);

        return analysis;

    }

    public void runMultipleTimepoints() {

        analyseGene analysis = new analyseGene(ancestralfilename, mainfilename);

        analysis.fixedNR = this.useFixedNeutralRatio;
        analysis.nr = this.neutralRatio;
        analysis.timepoints_multi = this.timepoints_multi;
        analysis.datasets = this.datasets;

//        analysis.neutral_ratio = new double[] {0.018};
//        analysis.neutral_ratio = new double[]{0.22909292816799146, 0.2624104494960995,0.14893244343173687,0.017836580830472952};
                //neutral_ratio = new double[]{0.2378342136223165, 0.25807525692320593, 0.14893244343173687}; exclude YRD2.1
                //neutral_ratio = new double[]{0.22909292816799146, 0.2624104494960995,0.06757118528744181};
//
//        //fluB
//        analysis.datasets = new String[] {"HA_yam_pink", "HA_yam_purple", "HA_vic_red", "HA_vic_blue"};
//
//        analysis.timepoints_multi.put("HA_yam_pink", new String [] {"1990","1992","1994","1996","1998","2000","2002","2004","2006","2008","2010","2012","2014"});
//        analysis.timepoints_multi.put("HA_yam_purple", new String [] {"1990","1992","1994","1996","1998","2000","2002","2004","2006","2008","2010","2012","2014"});
//        analysis.timepoints_multi.put("HA_vic_red", new String [] {"1990","1992","1994","1996","1998","2000","2002","2004","2006","2008","2010","2012","2014"});
//        analysis.timepoints_multi.put("HA_vic_blue", new String [] {"1990","1992","1994","1996","1998","2000","2002","2004","2006","2008","2010","2012","2014"});
//
//        analysis.neutral_ratio = new double[] {0.4086278689480355, 0.4006102116532152, 0.3332121714520498, 0.27419848421258597}; //full HA m-ratio
//        //analysis.neutral_ratio = new double[] {0.1596024989480543, 0.19231470087595393, 0.09949753779601428, 0.1625399773384175};  //int HA m-ratio
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


//        analysis.datasets = new String[] {"PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"};
//        analysis.timepoints_multi.put("PB2", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
//        analysis.timepoints_multi.put("PB1", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
//        analysis.timepoints_multi.put("PA", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
//        analysis.timepoints_multi.put("HA", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
//        analysis.timepoints_multi.put("NP", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
//        analysis.timepoints_multi.put("NA", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
//        analysis.timepoints_multi.put("MP", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
//        analysis.timepoints_multi.put("NS", new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});

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
////                System.out.println(d + "," + t + "," + value_matrix[t][d].r_high[d]);
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

    public analyseGene runBootstrapMultipleTimepoints(int bootstraps) {

        analyseGene analysis = new analyseGene(ancestralfilename, mainfilename);

        analysis.fixedNR = this.useFixedNeutralRatio;
        analysis.nr = this.neutralRatio;
        analysis.timepoints_multi = this.timepoints_multi;
        analysis.datasets = this.datasets;

        analysis.no_datasets = analysis.datasets.length;
        analysis.timepoints_per_dataset = new int[analysis.no_datasets];
        analysis.firstTimepoint = this.firstTimepoint;


        int max_timepoints = 0;
        for(int i=0; i < analysis.no_datasets; i++) {

            analysis.timepoints_per_dataset[i] = analysis.timepoints_multi.get(analysis.datasets[i]).length;
            if(analysis.timepoints_per_dataset[i]>max_timepoints) {
                max_timepoints = analysis.timepoints_multi.get(analysis.datasets[i]).length;
            }
        }


        analysis.no_timepoints = max_timepoints;
        analysis.bmAnalysisBootstrap(bootstraps);

        return analysis;

    }

    public void runDeepGenomeAnalysis() {

        analyseDeepGenome analysis = new analyseDeepGenome(ancestralfilename, mainfilename);

        analysis.datasets = this.datasets;
        analysis.genes = this.genes;
        analysis.timepoints = this.timepoints;
        analysis.firstTimepoint = this.firstTimepoint;
        analysis.fixedNR = this.useFixedNeutralRatio;
        analysis.nr = this.neutralRatio;
        analysis.bmAnalysis();
    }

    public analyseDeepGenome runBootstrapDeepGenomeAnalysis(int bootstraps) {

        analyseDeepGenome analysis = new analyseDeepGenome(ancestralfilename, mainfilename);

        analysis.datasets = this.datasets;
        analysis.no_datasets = analysis.datasets.length;
        analysis.timepoints_per_dataset = new int[analysis.no_datasets];
        analysis.genes = this.genes;
        analysis.timepoints = this.timepoints;
        analysis.firstTimepoint = this.firstTimepoint;

        for(int i=0; i < analysis.no_datasets; i++) {
            analysis.timepoints_per_dataset[i] = analysis.timepoints.length;
        }

        analysis.fixedNR = this.useFixedNeutralRatio;
        analysis.nr = this.neutralRatio;
        analysis.bmAnalysisBootstrap(bootstraps);

        return analysis;
    }

    public void readParams(String filename) {

        this.ancestralfilename = filename;
        this.mainfilename = filename;
    }

    public static void main(String [] args) {

        mainAnalysis mainAnalysis = new mainAnalysis();
        /*
         * initialise 'which analysis' flag as multi ('option 1 multiple timepoints')
         * values:
         * 	'multi-hcv'
         * 	'multi-flu'
         * 	'deep'
         * 	'one'
         */
        String whichAnalysis = "multi-hcv";		
        
        if(args.length==3){
            // read main and ancestral files from the command line args
        	mainAnalysis.ancestralfilename 	= args[0];
            mainAnalysis.mainfilename 		= args[1];
            whichAnalysis					= args[2].toLowerCase();
        }else{
        	// Jayna's hardcoded
        	mainAnalysis.ancestralfilename	= "./ancestralfile.txt";
            mainAnalysis.mainfilename		= "./mainfile.txt";
        }
        
        // For convenience, run a switch for the four analyses...
        switch(whichAnalysis){
        	// takes one of: {"multi-hcv","multi-flu","deep","one"}; default multi-hcv; fall-through break (no action);
        	case("multi-hcv"):{
                /** #1 run multiple timepoints, HCV **/

        		//HCV data
                mainAnalysis.datasets = new String[]{"p_53"};
                //dates are in units of year and relative to the first sample timepoint, which is zero or 0 years.
                mainAnalysis.timepoints_multi.put("p_53", new String[] {"0.3644", "0.6137", "0.8438", "1.3699", "1.7836", "3.8986", "6.8429", "7.6849"});
                mainAnalysis.firstTimepoint = new double[]{0.0};
                mainAnalysis.useFixedNeutralRatio = true;
                mainAnalysis.neutralRatio = new double[] {0.7186788};
                mainAnalysis.runMultipleTimepoints();
                analyseGene analysis = mainAnalysis.runBootstrapMultipleTimepoints(1);

                /* don't fall-through to next block(!) */
        		break;
       	}

        	case("multi-flu"):{
                /** #1 run multiple timepoints, flu **/

        		//initializing datasets
        		mainAnalysis.datasets = new String[]{"PRD1", "YRD1", "YRD2.2", "YRD2.1"};
        		//initializing timepoints per datasets
        		mainAnalysis.timepoints_multi.put("PRD1", new String[] {"2014.17", "2015.17", "2016.17", "2017.17"});
        		mainAnalysis.timepoints_multi.put("YRD1", new String[] {"2014.17", "2015.17"});
        		mainAnalysis.timepoints_multi.put("YRD2.2", new String[] {"2014.17", "2015.17", "2016.17", "2017.17"});
        		mainAnalysis.timepoints_multi.put("YRD2.1", new String[] {"2014.17", "2015.17", "2016.17", "2017.17"});
        		mainAnalysis.firstTimepoint = new double[]{2013.33,2013.33,2013.33,2013.33};
        		//estimating (useFixedNeutralRatio = false) or fixing the neutral ratio (useFixedNeutralRatio = true)
        		mainAnalysis.useFixedNeutralRatio = true;
        		//if useFixedNeutralRatio then need to specify neutral ratios for each dataset (neutral_ratio), which equals r_m/s_m
        		mainAnalysis.neutralRatio = new double[] {0.22909292816799146, 0.2624104494960995,0.14893244343173687,0.017836580830472952};
        		analyseGene analysis = mainAnalysis.runBootstrapMultipleTimepoints(10);      

        		/* don't fall-through to next block(!) */
        		break;
        	}

        	case("deep"):{
        		/** #3 run deepgenome analysis **/
        		mainAnalysis.datasets = new String[] {"gag", "pol", "env", "nef"};
        		mainAnalysis.genes = new int[][]{{294,1574},{1775,4519},{5800,7875},{8315,8605}};
        		mainAnalysis.timepoints = new String[]{"2006.5644","2008.7479","2011.4932"};
        		mainAnalysis.firstTimepoint = new double[]{2004.1315, 2004.1315, 2004.1315, 2004.1315};
        		mainAnalysis.useFixedNeutralRatio = true;
        		mainAnalysis.neutralRatio = new double[] {0.59754419,0.422687674,4.69888515,0.596970359};
        		mainAnalysis.runDeepGenomeAnalysis();
        		analyseDeepGenome analysis = mainAnalysis.runBootstrapDeepGenomeAnalysis(1);

        		/* don't fall-through to next block(!) */
        		break;
        	}

        	case("one"):{
                /** #2 run one timepoint analysis **/
                mainAnalysis.datasets = new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27"};
                mainAnalysis.useFixedNeutralRatio = true;
                mainAnalysis.neutralRatio = new double[] {4.813952552};
                mainAnalysis.timepoints = new String[]{"1"};
                mainAnalysis.runOneTimepoint();

                /* don't fall-through to next block(!) */
        		break;
        	}

        	default:{
        		System.err.println("No analysis flag set; use one of [multi-hcv | multi-flu | deep | one]");
        		break;
        	}
        }
     }
}


//analysis.writeoutBootstrapResults("/Users/jayna/Documents/Projects/HIV_ANPI/adaptation/PS133_adapt_bs.csv","a")
