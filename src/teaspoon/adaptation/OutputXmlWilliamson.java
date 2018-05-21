package teaspoon.adaptation;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

import teaspoon.app.utils.AncestralAlignmentParser;
import teaspoon.app.utils.AncestralAlignmentParserWithConsensus;
import teaspoon.app.utils.MainAlignmentParser;



public class OutputXmlWilliamson extends OutputXmlOther {
	public String[] headers = new String[7];


	public  OutputXmlWilliamson(){
		super();
		headers[0] = "Number Silent";headers[1] = "Number Replacement";headers[2] = "Total";headers[3] = "Silent/Replacement Ratio";
		headers[4] = "Neutral Ratio"; headers[5] = "Number Non Neutral";headers[6] = "Proportion Non-Neutral";

	}

	public void OutputWilliamson(String needconsensus){
		try{
			Element root = new Element("teaspoon.adaptation.Williamson-teaspoon.adaptation.Test-1991-Results");
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
			root.addContent(headerNames);	
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				// output sequence name and number
				Element sequence = new Element("Sequence");			// sequence name
				sequence.setAttribute("Location", data.get(m).toString());	// add sequence name
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
				Williamson W = new Williamson(seq,ans);
				double[][] value_matrix = W.williamson_method();
				float[][] printed_values = new float[value_matrix.length][value_matrix[0].length];
				for(int i=0;i<value_matrix.length;i++){
					for(int j=0;j<value_matrix[0].length;j++){
						printed_values[i][j] = (float) value_matrix[i][j];
					}
				}
				Element values = new Element("Values");	
				Element bin1 = new Element("B1");
				for(int i=0;i<value_matrix[0].length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[0][i]));
					bin1.addContent(node);
				}
				values.addContent(bin1);
				Element bin2 = new Element("B2");
				for(int i=0;i<value_matrix[0].length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[1][i]));
					bin2.addContent(node);
				}
				values.addContent(bin2);
				Element bin3 = new Element("B3");
				for(int i=0;i<value_matrix[0].length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[2][i]));
					bin3.addContent(node);
				}
				values.addContent(bin3);	
				Element bin4 = new Element("B4");
				for(int i=0;i<value_matrix[0].length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[3][i]));
					bin4.addContent(node);
				}
				values.addContent(bin4);	
				sequence.addContent(values);
				root.addContent(sequence);
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/teaspoon.adaptation.Williamson.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}  catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void OutputEyreWalker(String needconsensus){
		try{
			Element root = new Element("Eyre-Walker-teaspoon.adaptation.Test-1991-Results");
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
			root.addContent(headerNames);	
			for(int m=0;m<data.size();m++){		// Loop through all tiles ****
				// output sequence name and number
				Element sequence = new Element("Sequence");			// sequence name
				sequence.setAttribute("Location", data.get(m).toString());	// add sequence name
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
				Williamson W = new Williamson(seq,ans);
				double[][] value_matrix = W.eyrewalker_method();
				float[][] printed_values = new float[value_matrix.length][value_matrix[0].length];
				for(int i=0;i<value_matrix.length;i++){
					for(int j=0;j<value_matrix[0].length;j++){
						printed_values[i][j] = (float) value_matrix[i][j];
					}
				}
				Element values = new Element("Values");	
				Element bin1 = new Element("B1");
				for(int i=0;i<value_matrix.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i][0]));
					bin1.addContent(node);
				}
				values.addContent(bin1);
				Element bin2 = new Element("B2");
				for(int i=0;i<value_matrix.length;i++){
					Element node = new Element("V"+String.valueOf(i));
					node.addContent(String.valueOf(printed_values[i][1]));
					bin2.addContent(node);
				}
				values.addContent(bin2);
				sequence.addContent(values);
				root.addContent(sequence);			
			}
			Document document = new Document(root);				// create document on root
			XMLOutputter out = new XMLOutputter(Format.getPrettyFormat());		// use pretty formatting
			java.io.FileWriter writer = new java.io.FileWriter("/Users/sam/Desktop/Pybus bhatt web program/xmls/EW.xml");	// write to results.xml
			out.output(document, writer);
			writer.flush();
			writer.close();

		}  catch (Exception ex) {
			ex.printStackTrace();
		}
	}
}

