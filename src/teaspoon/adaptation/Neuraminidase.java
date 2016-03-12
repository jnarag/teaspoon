package teaspoon.adaptation;

public class Neuraminidase {
	 // where 0 is internal
	// and 1 is surface
	public int[] badlisth1n1 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,1,0,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,0,0,0,1,0,1,1,1,0,1,1,0,0,0};
    public int[]  badlisth3n2 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,1,0,0,1,1,0,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,0,1,1,1,1,1,0,1,1,1,1,1,0,1,0,0,0,1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0,1,0,1,1,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,1,0,1,0,1,1,1,1,0,1,1,1};
	public int[] badlisth1n1_NoStopCodon = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,1,0,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,1,1,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,0,0,0,1,0,1,1,1,0,1,1,0,0};
	public int[] FluB = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,1,1,1,1,1,0,1,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,0,1,0,1,0,0,0,0,1,1,0,1,1,1,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,1,0,1,0,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,1,1,1,1,1,0,1,0,0,0,1,0,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,0,0};
			public Neuraminidase(){
		
	}
}
