package jebl.evolution.alignments;

import java.util.Random;


/**
 * Date: 15/01/2006
 * Time: 10:13:50
 *
 * @author Joseph Heled
 * @version $Id: BootstrappedAlignment.java 940 2008-08-26 00:32:47Z stevensh $
 *
 */
public class BootstrappedAlignment extends ResampledAlignment {

    public BootstrappedAlignment(Alignment srcAlignment, long seed) {
        this(srcAlignment, new Random(seed));
    }
    
    public  BootstrappedAlignment(Alignment srcAlignment, Random r) {
        final int nSites = srcAlignment.getSiteCount();
        int[] sites = new int[nSites];

        for(int n = 0; n < nSites; ++n) {
            sites[n] = r.nextInt(nSites);
        }

        init(srcAlignment, sites);
    }

    public BootstrappedAlignment(Alignment srcAlignment) {
        this(srcAlignment, new Random());
    }
}
