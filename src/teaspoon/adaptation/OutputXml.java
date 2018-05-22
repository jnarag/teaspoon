package teaspoon.adaptation;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.output.*;

import teaspoon.app.utils.AncestralAlignmentParser;
import teaspoon.app.utils.AncestralAlignmentParserWithConsensus;
import teaspoon.app.utils.MainAlignmentParser;



public class OutputXml extends OutputXmlOther{

	private String[] headers = new String[9];
//	Constructers	



	public OutputXml(){
		super();	
		headers[0] = "Low-bound";headers[1] = "High-bound"; headers[2] = "Silent";headers[3] = "Replacemnt";
		headers[4] = "Total";headers[5]="R/S-Ratio";headers[6]="Z-Score";
		headers[7] = "Num-NN"; headers[8] = "Prop-NN";
	}



	public void OutputToXml(String needconsensus, Double number_bins) {

		// Start Formatting
		try{
			Element root = new Element("Pybus-Bhatt-results");		// root element
			Element headerNames = new Element("Headers");			// print header names for columns
			// add headers
			StringBuffer sbuf = new StringBuffer();
			for(int i=0;i<headers.length;i++){
				sbuf.append(headers[i]);							// use string buffer
				sbuf.append("  ");
			}
			headerNames.addContent(sbuf.toString());					
			Element Alignments = new Element("Alignments");
			Alignments.addContent(String.valueOf(data.size()));
			root.addContent(Alignments);
			root.addContent(headerNames);							// add header names to root
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				// output sequence name and number
				Element sequence = new Element("Sequence");			// sequence name
				sequence.setAttribute("Location", data.get(m).toString());	// add sequence name
				Element datamat = new Element("Values");		// siteData matrix
				MainAlignmentParser read = new MainAlignmentParser(data.get(m).toString());
				// this code is to check weather the ancestral sequence is an alignment or a single sequence
				int[] ans;
				if(needconsensus.equals("y")){	// if is consensus
					AncestralAlignmentParserWithConsensus readans = new AncestralAlignmentParserWithConsensus(ansdata.get(m).toString());
					ans = readans.read();
				} else {	//else 
					AncestralAlignmentParser readans = new AncestralAlignmentParser(ansdata.get(m).toString());
					ans = readans.read();
				}
				int[][] seq = read.read();
				SiteEstMulti Pb = new SiteEstMulti(seq,ans);
				Pb.SetNumBins(3.0);
				boolean[] vec = {false,false,false,false,true,true,false,false,false,false};
				double[][] matrix = Pb.Value_Matrix(vec, new double [10][], new double [10], false); //Pb.Value_Matrix(vec);	// values
				double[][] intervals = Pb.intervals();
				float[][] printed_intervals = new float[intervals.length][intervals[0].length];
				float[][] printed_matrix = new float[matrix.length][matrix[0].length];
				//  create float values
				for(int i=0;i<matrix.length;i++){
					for(int j=0;j<matrix[0].length;j++){
						Double temp = new Double(matrix[i][j]);
						printed_matrix[i][j] = temp.floatValue();	// printed matrix with float values
					}
				}
				for(int i=0;i<printed_intervals.length;i++){
					for(int j=0;j<printed_intervals[0].length;j++){
						Double temp = new Double(intervals[i][j]);
						printed_intervals[i][j] = temp.floatValue();
					}
				}
				
/*				
				Element node = new Element("V"+String.valueOf(i));
				node.addContent(String.valueOf(printed_matrix[i][j));
				values.addContent(node);
*/			
					
				for(int j=0;j<printed_matrix[0].length;j++){
					Element row = new Element("Row");	
					
					Element node1 = new Element("V"+String.valueOf(0));
					node1.addContent(String.valueOf(printed_intervals[0][j]));
					Element node2 = new Element("V"+String.valueOf(1));
					row.addContent(node1);
					row.addContent(node2);
					node2.addContent(String.valueOf(printed_intervals[1][j]));
					// add matrix content
					for(int i=0;i<printed_matrix.length;i++){
						Element node = new Element("V"+String.valueOf(i+2));
						node.addContent(String.valueOf(printed_matrix[i][j]));
						row.addContent(node);
					}
					row.setAttribute("Index",String.valueOf(j));
					datamat.addContent(row);						// add matrix
				}
				sequence.addContent(datamat);
				root.addContent(sequence);						// add to root
			}			
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/PB.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}



}