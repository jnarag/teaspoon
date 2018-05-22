package teaspoon.adaptation;

// siteData structure for a coalescet node
public class Node {
	double time;	// stores time of coalescent
	int nmuts;		// stores number of mutations
	int index;		// stores teaspoon.adaptation.Node number
	Node ancestor;	// pointer to ancestor
	Node desc1;		// pointer to descendent 1
	Node desc2;		// pointer to descendent 2
	

	public Node(int newindex){
		time = 0;
		nmuts = 0;
		desc1=null;
		desc2=null;
		ancestor = null;
		index = newindex;
	}
	

}


