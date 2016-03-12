package teaspoon.adaptation;

public class RunFile {
	public static void main(String[] args) { 
		if(args[0].equals("TD")){
			OutputXmlOther x = new  OutputXmlOther();
			x.OutputToXmlTajima();
	//		x.OutputToXmlCounts();
	//		x.OutputToXmlEstimators();
			OutputXmlCoalescent c = new OutputXmlCoalescent("n");
	//		c.OutputToXmlEstimators();
			c.OutputToXmlTajima();		
			x.OutputXmlSeqNames();
		}
		if(args[0].equals("FL")){
			OutputXmlOther x = new  OutputXmlOther();
			x.OutputToXmlCounts();
			x.OutputToXmlEstimators();
			x.OutputToXmlFuLi();
			OutputXmlCoalescent c = new OutputXmlCoalescent("n");
			c.OutputToXmlEstimators();
			c.OutputToXmlFuLi();	
			x.OutputXmlSeqNames();
		}
		if(args[0].equals("FW")){
			OutputXmlOther x = new  OutputXmlOther();
			x.OutputToXmlFayWu(args[1]);
			OutputXmlCoalescent c = new OutputXmlCoalescent("n");
			c.OutputToXmlFayWu();
			x.OutputXmlSeqNames();
		}
		if(args[0].equals("MK")){
			OutputXmlOther x = new  OutputXmlOther();
			x.OutputToXmlMK(args[1]);
			x.OutputXmlSeqNames();
		}
		if(args[0].equals("W")){
			OutputXmlWilliamson W = new OutputXmlWilliamson();
			W.OutputWilliamson(args[1]);
			W.OutputXmlSeqNames();
		}
		if(args[0].equals("EW")){
			OutputXmlWilliamson W = new OutputXmlWilliamson();
			W.OutputEyreWalker(args[1]);
			W.OutputXmlSeqNames();
		}
		if(args[0].equals("NAMES")){
			OutputXmlOther x = new  OutputXmlOther();
			x.OutputXmlSeqNames();
		}
		if(args[0].equals("PB")){
			OutputXml x = new OutputXml();
			x.OutputToXml(args[1],Double.valueOf(args[2]));
			x.OutputXmlSeqNames();
		}
	}
}
//System.out.println(System.getProperty("java.class.path"));

