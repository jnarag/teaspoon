package teaspoon.adaptation;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

import teaspoon.app.utils.AncestralAlignmentParser;
import teaspoon.app.utils.AncestralAlignmentParserWithConsensus;
import teaspoon.app.utils.MainAlignmentParser;

public class OutputXmlCoalescent extends OutputXmlOther{
	private String[] estimator_headers = new String[5];
	private String[] tajima_headers = new String[2];
	private String[] fuli_headers = new String[4];
	private String[] faywu_headers = new String[2];
	private String answer;

	public OutputXmlCoalescent(String needconsensus){
		super();	
		estimator_headers[0] = "Number of Segregating Sites";estimator_headers[1] = "ThetaS";
		estimator_headers[2] = "ThetaK";estimator_headers[3] = "Thetan";estimator_headers[4] = "Number of Singletons";
		tajima_headers[0] = "Low Confidence - 5%";tajima_headers[1] = "High Confidence - 95%";
		fuli_headers[0] = "D* Low Confidence - 5%";fuli_headers[1] = "D* High Confidence - 95%";
		fuli_headers[2] = "F* Low Confidence - 5%";fuli_headers[3] = "F* High Confidence - 95%";
		faywu_headers[0] = "D teaspoon.adaptation.Value"; faywu_headers[1] = "H teaspoon.adaptation.Value";
		answer = needconsensus;
	}


	public void OutputToXmlEstimators(){
		try{
			Element root = new Element("Estimators");	
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<estimator_headers.length;i++){
				sbuf.append(estimator_headers[i]);							// use string buffer
				sbuf.append("     ");
			}
			headerNames.addContent(sbuf.toString());					
			root.addContent(headerNames);	
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				MainAlignmentParser read = new MainAlignmentParser(data.get(m).toString());
				int[][] seq = read.read();
				int sample_size = seq.length;
				DiversityStats TD = new DiversityStats(seq);
				double theta  = TD.theta();						
				Element sequence = new Element("Simulation");	
				sequence.setAttribute("Sequence", String.valueOf(m));
				sequence.setAttribute("Sample-Size", String.valueOf(sample_size));
				sequence.setAttribute("Theta", String.valueOf(theta));
				CoalescentSim CS = new CoalescentSim(sample_size,theta);
				double[] values_array = CS.CoalescentEstimatorsMean();	// watterson estimates
				float[] printed_values = new float[values_array.length];
				for(int i=0;i<values_array.length;i++){
					printed_values[i] = (float) values_array[i];
				}
				Element values = new Element("Values");	
				for(int i=0;i<values_array.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i]));
					values.addContent(node);
				}
				sequence.addContent(values);
				root.addContent(sequence);	
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/CoalEst.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void OutputToXmlTajima(){
		try{
			Element root = new Element("Tajimas-D-Confidence-Interval");	
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<tajima_headers.length;i++){
				sbuf.append(tajima_headers[i]);							// use string buffer
				sbuf.append("     ");
			}
			headerNames.addContent(sbuf.toString());					
			root.addContent(headerNames);
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				MainAlignmentParser read = new MainAlignmentParser(data.get(m).toString());
				int[][] seq = read.read();
				int sample_size = seq.length;
				DiversityStats TD = new DiversityStats(seq);
				double theta  = TD.theta();	
				Element sequence = new Element("Simulation");	
				sequence.setAttribute("Sequence", String.valueOf(m));
				sequence.setAttribute("Sample-Size", String.valueOf(sample_size));
				sequence.setAttribute("Theta", String.valueOf(theta));
				CoalescentSim CS = new CoalescentSim(sample_size,theta);
				double[] values_array = CS.TajimaCS();	// watterson estimates
				float[] printed_values = new float[values_array.length];
				for(int i=0;i<values_array.length;i++){
					printed_values[i] = (float) values_array[i];
				}
				Element values = new Element("Values");	
				for(int i=0;i<values_array.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i]));
					values.addContent(node);
				}
				sequence.addContent(values);
				root.addContent(sequence);
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/CoalTajima.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void OutputToXmlFuLi(){
		try{
			Element root = new Element("Fu-Li-Confidence-Interval");	
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<fuli_headers.length;i++){
				sbuf.append(fuli_headers[i]);							// use string buffer
				sbuf.append("     ");
			}
			headerNames.addContent(sbuf.toString());					
			root.addContent(headerNames);	
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				MainAlignmentParser read = new MainAlignmentParser(data.get(m).toString());
				int[][] seq = read.read();
				int sample_size = seq.length;
				DiversityStats TD = new DiversityStats(seq);
				double theta  = TD.theta();	
				Element sequence = new Element("Simulation");	
				sequence.setAttribute("Sequence", String.valueOf(m));
				sequence.setAttribute("Sample-Size", String.valueOf(sample_size));
				sequence.setAttribute("Theta", String.valueOf(theta));
				CoalescentSim CS = new CoalescentSim(sample_size,theta);
				double[] values_array = CS.FuAndLiCS();	// watterson estimates
				float[] printed_values = new float[values_array.length];
				for(int i=0;i<values_array.length;i++){
					printed_values[i] = (float) values_array[i];
				}
				Element values = new Element("Values");	
				for(int i=0;i<values_array.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i]));
					values.addContent(node);
				}	
				sequence.addContent(values);
				root.addContent(sequence);		
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/CoalFuLi.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void OutputToXmlFayWu(){
		try{
			Element root = new Element("Fay-Wu-Confidence-Interval");	
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<faywu_headers.length;i++){
				sbuf.append(faywu_headers[i]);							// use string buffer
				sbuf.append("     ");
			}
			headerNames.addContent(sbuf.toString());					
			root.addContent(headerNames);	
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				MainAlignmentParser read = new MainAlignmentParser(data.get(m).toString());
				// this code is to check weather the ancestral sequence is an alignment or a single sequence
				int[] ans;
				if(answer.equals("y")){	// if is consensus
					AncestralAlignmentParserWithConsensus readans = new AncestralAlignmentParserWithConsensus(ansdata.get(m).toString());
					ans = readans.read();
				} else {	//else 
					AncestralAlignmentParser readans = new AncestralAlignmentParser(ansdata.get(m).toString());
					ans = readans.read();
				}
				int[][] seq = read.read();
				int sample_size = seq.length;
				FayAndWu  FW  = new FayAndWu(seq,ans);
				DiversityStats TD = new DiversityStats(seq);
				double theta  = TD.theta();	
				double[] actualvalue = FW.FayWu();
				Element sequence = new Element("Simulation");	
				sequence.setAttribute("Sequence", String.valueOf(m));
				sequence.setAttribute("Sample-Size", String.valueOf(sample_size));
				sequence.setAttribute("Theta", String.valueOf(theta));
				CoalescentSim CS = new CoalescentSim(sample_size,theta);
				double[] values_array = CS.FayAndWuCS(actualvalue);	// watterson estimates
				float[] printed_values = new float[values_array.length];
				for(int i=0;i<values_array.length;i++){
					printed_values[i] = (float) values_array[i];
				}
				Element values = new Element("Values");	
				for(int i=0;i<values_array.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i]));
					values.addContent(node);
				}	
				sequence.addContent(values);
				root.addContent(sequence);		
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/CoalFayWu.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}catch (Exception ex) {
			ex.printStackTrace();
		}
	}
}
