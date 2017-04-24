package teaspoon.adaptation; /**
 * Created with IntelliJ IDEA.
 * User: jayna
 * Date: 10/07/2013
 * Time: 15:49
 * To change this template use File | Settings | File Templates.
 */

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class CodeTest {


    public void runMethodByTime() {
        Read_main ancestral = new Read_main("/Users/jayna/Documents/Projects/adapt-a-rate/from_sam/Bhatt Holmes Pybus - The genomic rate of molecular adaptation of the human influenza A virus - Alignments/H1N1/processed_1977_1979.M2.H1N1.nex");


        try {
            BufferedReader reader = new BufferedReader(new FileReader("filelist.txt"));

            int count = 1;
            while(reader.ready()) {

                Read_main main = new Read_main(reader.readLine().trim());

                int [][] ancestralMatrix = ancestral.readNEXUS();
                int [][] mainMatrix = main.readNEXUS();

                int [] ans = ancestral.consensusArray(ancestralMatrix);
                int [] main_consensus = main.consensusArray(mainMatrix);

//                System.out.println(mainMatrix[0].length);
//                System.out.println(ancestralMatrix[0].length) ;
//
//                teaspoon.adaptation.McDonaldKreitman mk = new teaspoon.adaptation.McDonaldKreitman(mainMatrix, ans);
//                double fixed = mk.howManyFixed();
//                double poly = mk.howManySingPoly();
//
//                double [][] contingencyTable = mk.createContingencyNew();
//
//                System.out.println(fixed);
//                System.out.println(poly);
//
//                teaspoon.adaptation.Williamson w = new teaspoon.adaptation.Williamson(mainMatrix, ans);
//
//                fixed = w.howManyFixed();
//                poly = w.howManySingPoly();
//                double multiPoly = w.howManyMultiPoly();
//
//                System.out.println(fixed);
//                System.out.println(poly);
//                System.out.println(multiPoly);

                BhattMethod bm = new BhattMethod(mainMatrix, ans);

                double[] L = {0.0, 0.15, 0.75};//{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
                double[] H = {0.15, 0.75, 1.0};//{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
                boolean[] Nvec = {false,true,false};//{false,false,false,true,true,true,true,false,false,false};

                double[][] bins = new double[2][L.length];
                for(int i=0;i<L.length;i++){
                    bins[0][i]=L[i];
                    bins[1][i]=H[i];
                }
                double[] prior = {1.0,1.0,1.0,1.0};

                boolean[] neutralVec = {false,false,false,false,true,true,false,false,false,false};
                Williamson3bin w3b = new Williamson3bin(mainMatrix,ans);

                double [] low = {0.0, 0.15};
                double [] mid = {0.15, 0.75};
                double [] high = {0.15, 1.0};

                w3b.williamson3bin_method(low, mid, high);
                w3b.williamson3bin_method(1.0, low, mid, high);
//
//
//                System.out.println(w3b.Adapt);
                System.out.println("mid R: "+w3b.mid_R);
                System.out.println("mid S: "+w3b.mid_S);


                bm.Method(bins,prior,true, Nvec);

                System.out.println("Replacement/Silent Ratio");
                bm.print(bm.ReplacementSilentRatio);
                System.out.println("Replacement Count Array");
                bm.print(bm.ReplacementCountArray);
                System.out.println("Silent Count Array");
                bm.print(bm.SilentCountArray);
                System.out.println("Total Count Array");
                bm.print(bm.TotalCountArray);

//
                System.out.println("Non-neutral substitutions");
                bm.print(bm.NonNeutralSubstitutions);
                System.out.println("neutral ratio");
                System.out.println(bm.neutralratio);
//        //System.out.println(bm.neutralbin);
                System.out.println("adaptive sites");
                System.out.println(count +" : "+bm.Adaptation);
                count++;
                System.out.println();
                System.out.println();
            }


        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }
    public static void main (String [] args) {

        CodeTest test = new CodeTest();
        test.runMethodByTime();

    }

    public void execMultiFiles(String fileList) {

        try {
            BufferedReader reader = new BufferedReader(new FileReader(fileList));

            while(reader.ready()) {

                String file = reader.readLine().trim();
                Read_main main = new Read_main(file);




            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
}
