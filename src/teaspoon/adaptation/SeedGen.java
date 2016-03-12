package teaspoon.adaptation;

import java.util.Random;

public class SeedGen {


	public static void main(String[] args) { 
		 Random generator = new Random();
		 Long r = Math.abs(generator.nextLong()) ;
		 System.out.println(r);
	}

	
}
