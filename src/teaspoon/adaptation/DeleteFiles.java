package teaspoon.adaptation;

import java.io.File;

public class DeleteFiles {

public static void main(String[] args) { 
		File input = new File("/Users/sam/Desktop/Pybus bhatt web program/seq/");												// Input 
		File[] data = input.listFiles();
		for(int i=0;i<data.length;i++){
			data[i].delete();
		}
		input.delete();
		
		File input1 = new File("/Users/sam/Desktop/Pybus bhatt web program/xmls/");		
		File[] data1 = input1.listFiles();
		for(int i=0;i<data1.length;i++){
			data1[i].delete();
		}
		input1.delete();
		
		File input2 = new File("/Users/sam/Desktop/Pybus bhatt web program/ans/");		
		File[] data2 = input2.listFiles();
		for(int i=0;i<data2.length;i++){
			data2[i].delete();
		}
		input2.delete();
		
	}
	
}
