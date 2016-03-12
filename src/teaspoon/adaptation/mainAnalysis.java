package teaspoon.adaptation;


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

//        this.ancestralFile = ancestralFile;
//        this.mainFile = mainFile;
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

        analysis.fixedNR = false;
        analysis.datasets = new String[] {"p_4", "p_37", "p_53", "p_61"};
        analysis.timepoints_multi.put("p_4", new String[] {"2002.4274", "2003.0192", "2005.2055", "2006.0932", "2006.4521", "2007.1233", "2008.1366", "2009.1973", "2010.2849", "2012.4891", "2013.4521"});
        //timepoints_multi.put("p_37", new String[] {"2002.377", "2002.6055", "2002.874", "2003.2", "2004.2077", "2005.4027", "2006.3178", "2007.3507"});
        analysis.timepoints_multi.put("p_37", new String[] {"2002.0877", "2002.3370", "2002.6055", "2002.874", "2003.2", "2004.2077", "2005.4027", "2006.3178", "2007.3507"});
        analysis.timepoints_multi.put("p_53", new String[] {"2001.5041", "2001.8685", "2002.1178", "2002.3479", "2002.874", "2003.2877", "2005.4027", "2008.347", "2009.1890"});
        analysis.timepoints_multi.put("p_61", new String[] {"2002.0192", "2002.3863", "2002.6329", "2003.0986", "2003.7644", "2004.6667", "2005.011", "2006.0164", "2007.0274", "2008.9481", "2009.074", "2009.5425", "2010.3753", "2012.1858","2013.6438"});


        analysis.no_datasets = analysis.datasets.length;
        analysis.timepoints_per_dataset = new int[analysis.no_datasets];


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

        mainAnalysis.ancestralfilename = "/Users/jayna/Documents/Projects/AMC_HCV_DATA/ancestral_HCVpacbio_filelist.txt";
        mainAnalysis.mainfilename = "/Users/jayna/Documents/Projects/AMC_HCV_DATA/main_HCVpacbio_filelist.txt";

        mainAnalysis.runBM_multipleTimepoints();

    }



}