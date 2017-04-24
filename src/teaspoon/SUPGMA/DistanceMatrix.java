package teaspoon.SUPGMA;


import java.util.List;

//import jebl.evolution.taxa.Taxon;
import jebl.evolution.taxa.Taxon;
import teaspoon.adaptation.Read_main;

/**
 * Created by jayna on 12/02/2014.
 */
public class DistanceMatrix {

    int [][] sequenceMatrix;
    List<Double> sampleTimes;
    List<String> sampleNames;
    List<Taxon> taxa;
    Read_main main;

    double constA;
    double constB;
    double constC;

    protected static final double MAX_DISTANCE = 1000.0;

    public DistanceMatrix(int[][] sequenceMatrix) {

        this.sequenceMatrix = sequenceMatrix;


    }


    public DistanceMatrix(String filename) {

        main = new Read_main(filename, true);
        //re-write method so only rely on readFasta to read in sequences (convert nexus/phylip files to fasta)
        sequenceMatrix = main.readFASTA();
        sampleTimes = main.sampleTimes;
        sampleNames = main.sampleNames;
        taxa = main.taxa;



    }

    public List<Taxon> getTaxa() {
        return taxa;
    }

    public double[][] getRawDistances() {

        double[][] rawDistances = new double[sequenceMatrix.length][sequenceMatrix.length];


        for(int i=0; i< sequenceMatrix.length-1; i++) {

            int[] seq1 = sequenceMatrix[i];

            for(int j=(i+1); j<sequenceMatrix.length; j++) {

                int seq_diff = 0;


                int[] seq2 = sequenceMatrix[j];

                for(int s=0; s<seq1.length; s++) {

                    int base_i = seq1[s];
                    int base_j = seq2[s];

                    if(base_i != 5 || base_j != 5) {   //ignore gaps
                        if(base_i != base_j) {

                            //System.out.println(s+": "+base_i+","+base_j);
                            seq_diff++;

                        }
                    }
                }

                double perSiteDiff = (double)seq_diff/seq1.length;

                rawDistances[i][j] = perSiteDiff;



            }
        }
        return rawDistances;
    }

    public double[][] getJCDistances() {

        double[][] JCCorrectedDistances = new double[sequenceMatrix.length][sequenceMatrix.length];


        for(int i=0; i< sequenceMatrix.length-1; i++) {

            int[] seq1 = sequenceMatrix[i];

            for(int j=(i+1); j<sequenceMatrix.length; j++) {

                int seq_diff = 0;


                int[] seq2 = sequenceMatrix[j];

                for(int s=0; s<seq1.length; s++) {

                    int base_i = seq1[s];
                    int base_j = seq2[s];

                    if(base_i != 5 || base_j != 5) {   //ignore gaps
                        if(base_i != base_j) {

                            //System.out.println(s+": "+base_i+","+base_j);
                            seq_diff++;

                        }
                    }
                }

                double perSiteDiff = (double)seq_diff/seq1.length;
                double correctedDiff = 0.0;

                if(perSiteDiff != 0.0){
                    // use constants?
                    correctedDiff =-(3.0/4.0)*Math.log(1.0-(4.0/3.0)*perSiteDiff);
                }
                JCCorrectedDistances[i][j] = correctedDiff;



            }
        }
        return JCCorrectedDistances;

    }

    public int calculateMatrixDimension() {

        int count = 0;

        int dimensions = (int)Math.floor((Math.pow(sequenceMatrix.length, 2)-sequenceMatrix.length)*0.5);

        for(int i = 0; i< sequenceMatrix.length-1; i++) {

            for(int j = (i+1); j < sequenceMatrix.length; j++) {
                if(sampleTimes.get(i)>=sampleTimes.get(j)) {


                    count++;
                }
            }

        }
        //System.out.println("number of times t_m > t_n "+count);
        return dimensions;
    }

    public double[][] getHKYDistanceMatrix() {


        getHKYconstants();

        double [][] HKYcorrectedDistances = new double[sequenceMatrix.length][sequenceMatrix.length];

        for(int i = 0; i< sequenceMatrix.length-1; i++) {

            for(int j = (i+1); j < sequenceMatrix.length; j++) {

                int taxon1 = i;
                int taxon2 = j;

                double correctedDistance = getHKYDistance(taxon1, taxon2);

                HKYcorrectedDistances[i][j] = correctedDistance;

            }
        }

        return HKYcorrectedDistances;
    }

    public double[][] getHKYDistanceMatrix(int[][] sequenceMatrix) {

        double [][] HKYcorrectedDistances = new double[sequenceMatrix.length][sequenceMatrix.length];
        getHKYconstants();

        int count = 1;
        for(int i = 0; i< sequenceMatrix.length-1; i++) {

            for(int j = (i+1); j < sequenceMatrix.length; j++) {

                //System.out.println(i+","+j);

                HKYcorrectedDistances[i][j] = getHKYDistance(i, j);

            }
        }

        return HKYcorrectedDistances;
    }

    public List<Double> getSampleTimes() {

        return sampleTimes;
    }

    public List<String> getSampleNames() {

        return sampleNames;
    }

    public int[][] getSequenceMatrix() {

        return sequenceMatrix;
    }
    private double getHKYDistance(int taxon1, int taxon2) {

        double sumTs = 0.0; // total weight of all columns that have a transition for these taxa
        double sumTv = 0.0; // total weight of all columns that have a transversion for these taxa
        double sumWeight = 0.0; // total weight of all columns (ignoring those with ambiguities, but
        // including identical columns (which have neither a transition nor a transversion) )
        boolean noGapsPairFound = false;

        int[] seq1 = sequenceMatrix[taxon1];
        int[] seq2 = sequenceMatrix[taxon2];
        for(int i=0; i < seq1.length; i++ ) {

            int state1 = seq1[i];
            int state2 = seq2[i];

            // ignore any ambiguous or gaps
            if( state1 == 5 || state2 == 5 ) {
                continue;
            } else {
                noGapsPairFound = true;
            }

            double weight = 1.0;
            // acgt
            if ( state1 != state2 ) {
                if (isTransition(state1, state2) ) {
                    // it's a transition
                    sumTs += weight;
                } else {
                    // it's a transversion
                    sumTv += weight;
                }
            }
            sumWeight += weight; // this also includes the columns with state1 == state2
        }

        if( sumWeight <= 0.0 ) {
            return 0.0;
        }

        while( true ) {

            double P = sumTs / sumWeight;
            double Q = sumTv / sumWeight;

            double a = 1.0 - (P / (2.0 * constA)) - (((constA - constB) * Q) / (2.0 * constA * constC));

            if( a <= 0 ) {
                // minimum number of sites whose removal restores consistency. see comments
                // in TamuraNei.
                final int adjustment = (int)(1 + (sumWeight * -a) / (1/(2.0*constA)  - 1));
                sumTs -= adjustment;
                if( sumTs < 0) {
                    break;
                }
                sumWeight -= adjustment;
                continue;
            }

            double b = 1.0 - (Q / (2.0 * constC));
            if( b < 0 ) {
                break;
            }

            final double distance = -(2.0 * constA * Math.log(a)) + (2.0 * (constA - constB - constC) * Math.log(b));

            return Math.min(distance, MAX_DISTANCE);
        }

        return MAX_DISTANCE;
    }

    private void getHKYconstants() {


        double [] freqs = getFrequencies(sequenceMatrix); //get frequencies

        double freqA = freqs[0];
        double freqC = freqs[1];
        double freqG = freqs[2];
        double freqT = freqs[3];

        double freqR = freqA+freqG;
        double freqY = freqC+freqT;

        constA = ((freqA * freqG) / freqR) + ((freqC * freqT) / freqY);
        constB = (freqA * freqG) + (freqC * freqT);
        constC = (freqR * freqY);


    }

    private double[] getFrequencies(int[][] sequenceMatrix) {

        int A = 0;
        int C = 0;
        int G = 0;
        int T = 0;

        for(int i=0; i< sequenceMatrix.length; i++) {

            for(int j=0; j<sequenceMatrix[i].length; j++) {

                int site = sequenceMatrix[i][j];

                switch(site) {

                    case 1:
                        A++;
                        break;
                    case 2:
                        C++;
                        break;
                    case 3:
                        G++;
                        break;
                    case 4:
                        T++;
                        break;


                }

            }

        }

        double total = A+C+G+T;
        double freqA = A/total;
        double freqC = C/total;
        double freqG = G/total;
        double freqT = T/total;


        return new double[] {freqA, freqC, freqG, freqT};

    }

    private boolean isTransition(int state1, int state2) {

        if(state1 == 1 && state2 == 3) {

            return true;
        }
        else if(state1 == 2 && state2 == 4) {

            return true;
        }
        else if(state1 == 3 && state2 == 1) {
            return true;
        }
        else if(state1 == 4 && state2 == 2) {
            return true;
        }
        else{
            return false;
        }
    }


 //   public static void main(String [] args) {

//        String filename = "cate_MSHackbsp10g10d_part_day_0.fasta";
//
//        teaspoon.SUPGMA.DistanceMatrix dist = new teaspoon.SUPGMA.DistanceMatrix(filename);
//
//        double[][] rawDist = dist.getHKYDistanceMatrix();
//        StringBuilder builder = new StringBuilder();
//
//        for(int i=0; i<rawDist.length; i++) {
//
//
//            for(int j=0; j < rawDist.length; j++) {
//
//                builder.append(rawDist[i][j]).append(",");
//
//            }
//            builder.append("\n");
//        }
//
//        System.out.println(builder.toString());
//    }



}
