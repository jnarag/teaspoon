package jebl.evolution.alignments;

import java.util.Random;


/**
 * Date: 17/01/2006
 * Time: 08:18:32
 *
 * @author Joseph Heled
 * @version $Id: JackknifedAlignment.java 940 2008-08-26 00:32:47Z stevensh $
 *
 */
public class JackknifedAlignment extends ResampledAlignment {
    public JackknifedAlignment(Alignment srcAlignment, double percent){
        this(srcAlignment, percent, new Random());
    }

    public JackknifedAlignment(Alignment srcAlignment, double percent, long seed){
        this(srcAlignment, percent, new Random(seed));
    }

     public JackknifedAlignment(Alignment srcAlignment, double percent, Random r) {
        final int nSites = srcAlignment.getSiteCount();
        final int nNewSites = (int)Math.ceil(nSites * percent);
        int[] sites = new int[nSites];

        for(int n = 0; n < nSites; ++n) {
            sites[n] = n;
        }

         shuffle(sites, r);

        int[] newSites = new int[nNewSites];
        System.arraycopy(sites, 0, newSites, 0, nNewSites);
        init(srcAlignment, newSites);
    }

    /**
     * Shuffles an array.
     */
    private void shuffle(int[] array, Random random) {
        int l = array.length;
        for (int i = 0; i < l; i++) {
            int index = random.nextInt(l-i) + i;
            int temp = array[index];
            array[index] = array[i];
            array[i] = temp;
        }
    }
}
