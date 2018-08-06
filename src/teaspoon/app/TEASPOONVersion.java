/**
 * 
 */
package teaspoon.app;


/**
 * 
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 2 Aug 2018
 * @version 0.1
 */
public class TEASPOONVersion  {

    /**
     * Version string: assumed to be in format x.x.x
     */
    private static final String VERSION = "0.1.4";

    private static final String DATE_STRING = "2010-2018";

    private static final boolean IS_PRERELEASE = true;

    private static final String TEASPOON_WEBPAGE = "https://github.com/jnarag/teaspoon/blob/master/README.md";
    
    private static final String TEASPOON_SOURCE = "https://github.com/jnarag/teaspoon";
    
    private static final String CITATION = "Jayna Raghwani, Samir Bhatt, Oliver G. Pybus, 2016. Faster Adaptation in Smaller Populations: Counterintuitive Evolution of HIV during Childhood Infection. PLoS Comp. Biol. Jan 7th 2016";
    
    private static final String CITATION_URL = "http://dx.doi.org/10.1371/journal.pcbi.1004694";

	public TEASPOONVersion() {
		// TODO Auto-generated constructor stub
	}

	public String getVersion() {
		// TODO Auto-generated method stub
		return VERSION;
	}

	/* (non-Javadoc)
	 * @see beast.app.util.Version#getVersionString()
	 */
	public String getVersionString() {
        return "v" + VERSION + (IS_PRERELEASE ? " Prerelease" : "");
	}

	/* (non-Javadoc)
	 * @see beast.app.util.Version#getDateString()
	 */
	public String getDateString() {
		return DATE_STRING;
	}

	/* (non-Javadoc)
	 * @see beast.app.util.Version#getCredits()
	 */
	public String[] getCredits() {
        return new String[]{
        		this.getVersionString(),
        		"",
                "Designed and developed by",
                "Jayna Raghwani, Joe Parker & Oliver G. Pybus",
                "",
                "Department of Zoology, University of Oxford, South Parks Road, Oxford, United Kingdom",
                "jayna.raghwani@zoo.ox.ac.uk",
                "",
                "joe@kitserve.org.uk",
                "",
                "oliver.pybus@zoo.org.uk",
                "",
                "Including code and analyses designed by",
                "Samir Bhatt",
                "",
                "Downloads, Help & Resources:",
                TEASPOON_WEBPAGE,
                "",
                "Source code distributed under the GNU Lesser General Public License:",
                TEASPOON_SOURCE,
                "",
                "Please cite this software if used in any academic work:",
                CITATION,
                CITATION_URL,
                "",
                "TEASPOON GUI developers:",
                "Joe Parker, Jayna Raghwani",
                "",
                "",
                "Incorporating third-party code under various licences:",
                "Colt - Open Source Libraries for High Performance Scientific and Technical Computing (colt.jar), ",
                "JAMA - Java Matrix Algebra (1.0.3, 2012.11.09; Jama-1.0.3.jar), ",
                "JSC - Java Statistical Classes v1.0 (jsc-1.jar), ",
                "JFreeChart (jfreechart-1.0.14.jar), ",
                "XChart (xchart-3.3.1.jar), ",
                "JEBL - Java Evolutionary Biology Library v0.4 (jebl-0.4.jar)"};
	}

	public String getTextCredits() {
	    String sStr = "TEASPOON:Tools for Evolutionary Analysis of Serially-sampled POpulatiONs\n";
        for (String s : getCredits()) {
        	sStr += "\n" + s;
        }
        return sStr;
	}
	public String getHTMLCredits() {
        String sStr = "<H2>TEASPOON:</H2><i><b>T</b>ools for <b>E</b>volutionary <b>A</b>nalysis of <b>S</b>erially-sampled <b>PO</b>pulati <b>ON</b>s</i><hr/>";
        for (String s : getCredits()) {
            if (s.contains("@")) {
                s = "<a href=\"mailto:" + s + "\">" + s + "</a><br>";
            }
            if (s.contains("jar")) {
                s = "<i>" + s + "</i>";
            }
            if (s.contains("http")) {
                sStr += "<a href=\"" + s + "\">" + s + "</a><br>";
            } else {
                sStr += "<p>" + s + "</p>";
            }
        }
        return sStr;
    }

    public String getMajorVersion() {
        return VERSION.substring(0, VERSION.lastIndexOf("."));
    }
}
