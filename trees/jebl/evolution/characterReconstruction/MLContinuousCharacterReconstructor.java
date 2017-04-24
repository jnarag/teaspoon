package jebl.evolution.characterReconstruction;

import java.util.HashMap;
import java.util.Map;

import jebl.math.MatrixCalc;

/**
 * This is the main class that controls the ancestral calculations by the Schluter ML method
 * The theoretical foundations are found in Schulter et al. 1997 Evolution.
 * @author Stephen A. Smith
 */

public class MLContinuousCharacterReconstructor {
    /**
     * a basic constructor which just defines the tree and character
     * @param conchar the continuous character to reconstruct
     * @param tree the tree on which to reconstruct the character
     */
    public MLContinuousCharacterReconstructor(jebl.evolution.characters.ContinuousCharacter conchar, jebl.evolution.trees.SimpleRootedTree tree){
        this.conchar = conchar;
        this.tree = tree;
    }

    /**
     * public
     */
    /**
     * @param tree the tree for the mainAnalysis, unnecessary if using the same tree used from the constructor
     */
    public void setTree(jebl.evolution.trees.SimpleRootedTree tree){
        this.tree = tree;
    }

    /**
     * @param chonchar the character for the analsysis, unnecessary if using the same tree used from the constructor
     */
    public void setChar(jebl.evolution.characters.ContinuousCharacter conchar){
        this.conchar = conchar;
    }

    /**
     * initiates the mainAnalysis
     */
    public void reconstruct(){
        phenotypes = new HashMap<jebl.evolution.graphs.Node, Double>();
        phenotypesSE = new HashMap<jebl.evolution.graphs.Node, Double>();
        internalNodeNums = new HashMap<jebl.evolution.graphs.Node, Integer>();
        beginFill();
        beginStats();
    }

    /**
     * @return the phenotypes estimated from the mainAnalysis
     */
    public Map<jebl.evolution.graphs.Node, Double> getPhenotypes(){return phenotypes;}

    /**
     * @return the phenotype standard errors estimated from the mainAnalysis
     */
    public Map<jebl.evolution.graphs.Node, Double> getPhenotypesSE(){return phenotypesSE;}

    /**
     * private
     */
    /**
     * initiate some variables and arrays that are used throughout the mainAnalysis
     */
    private void beginFill(){
        //init some things to make the code a little easier
        Object [] TinternalNodes = tree.getInternalNodes().toArray();
        numInNodes = TinternalNodes.length;
        internalNodes = new jebl.evolution.graphs.Node[numInNodes];
        Object [] TccTaxa = conchar.getTaxa().toArray();
        for(int i=0;i<TccTaxa.length;i++){
            phenotypes.put(tree.getNode((jebl.evolution.taxa.Taxon)TccTaxa[i]),(Double)conchar.getValue((jebl.evolution.taxa.Taxon)TccTaxa[i]));
            phenotypesSE.put(tree.getNode((jebl.evolution.taxa.Taxon)TccTaxa[i]),(Double)conchar.getSE((jebl.evolution.taxa.Taxon)TccTaxa[i]));
        }
        for(int i=0;i<numInNodes;i++){
            internalNodes[i] = (jebl.evolution.graphs.Node)TinternalNodes[i];
            phenotypes.put(internalNodes[i],0.0);
            phenotypesSE.put(internalNodes[i],0.0);
            internalNodeNums.put(internalNodes[i], i);
        }
    }

    /**
     * the meat of the mainAnalysis which calculates the estimates of both ML estimates for each node
     * as well as the SE estimates for each node
     */
    private void beginStats(){
        //MLE
        double mlsos = getMLSofS(tree);
        //estimate BETA
        double df = numInNodes;//degrees of freedom, with one character this is num of tips - 1
        double [] tempSE = new double [numInNodes];
        for(int i=0; i < numInNodes; i++){
                tempSE[i]= getSEest(i);
                phenotypesSE.put(internalNodes[i],Math.sqrt(2*mlsos/(df*tempSE[i])));
        }
    }

    /**
     * evaluate ML residual sum of squares
     */
    private double getMLSofS(jebl.evolution.trees.Tree tp) {
        computeMLest();
        //global qsum for ease
        qsum = 0;
        addToQ(tree.getRootNode()); //sum of squares
        return qsum;
    }

    /**
     * compute the ML estimates for all nodes
     */
    private void computeMLest(){
        fullM = new double [numInNodes][numInNodes];
        double [] fullRhs = new double [numInNodes];
        fullMcp=fullM;
        fullVcp=fullRhs;
        // calculate q matrix
        for(int i=0; i<numInNodes; i++){
                jebl.evolution.graphs.Node p=internalNodes[i];
                doQCalc(tree.getChildren(p).get(0),i);
                doQCalc(tree.getChildren(p).get(1),i);
        }
        //save a copy for SE calculation
        SEfullM = MatrixCalc.copyMatrix(fullM);
        //Cholesky factorization and solution
        double [][] CHfact=fullM;
        try{
                CHfact=MatrixCalc.choleskyFactor(CHfact);
            }catch(jebl.math.MatrixCalcException.NotSquareMatrixException npe){}
            catch (jebl.math.MatrixCalcException.PositiveDefiniteException pde){};
            double [] rhs=fullRhs;
            double [] mle = null;
            try{
                mle=MatrixCalc.choleskySolve(CHfact,rhs);
            }catch(jebl.math.MatrixCalcException.NotSquareMatrixException npe){};
            // set values at nodes
            for(int i=0; i<numInNodes; i++){
                jebl.evolution.graphs.Node p=internalNodes[i];
                phenotypes.put(p,mle[i]);
        }
    }

    /**
     * compute and return the standard error estimates for a particular node
     */
    private double getSEest(int nodeNum){
        int ir = nodeNum;
        double qpp = SEfullM[ir][ir];
        double [][] oneLessRow = MatrixCalc.deleteMatrixRow(SEfullM,ir);
        double [] grabDcol = MatrixCalc.getColumn(oneLessRow,ir);
        double [][] reducedM = MatrixCalc.deleteMatrixColumn(oneLessRow,ir);
        try{
                reducedM = MatrixCalc.choleskyFactor(reducedM);
        }catch(jebl.math.MatrixCalcException.NotSquareMatrixException npe){}
        catch (jebl.math.MatrixCalcException.PositiveDefiniteException pde){};
            double [] CHsol = null;
            try{
                CHsol = MatrixCalc.choleskySolve(reducedM,grabDcol);
            }catch(jebl.math.MatrixCalcException.NotSquareMatrixException npe){};
        double tempSE = qpp - MatrixCalc.innerProduct(grabDcol,CHsol,0);
        return tempSE;
    }

    /**
     *  add in values corresponding to one node
     */
    private void doQCalc(jebl.evolution.graphs.Node q,int ip){
        double tbl=2/tree.getLength(q);
        fullMcp[ip][ip] += tbl;
        if(tree.isExternal(q)){
            fullVcp[ip] += phenotypes.get(q)*tbl;
        }
        else {
            int iq=internalNodeNums.get(q);
            fullMcp[ip][iq] -= tbl;
            fullMcp[iq][ip] -= tbl;
            fullMcp[iq][iq] += tbl;
        }
    }

    private void addToQ(jebl.evolution.graphs.Node p){
        if(tree.isExternal(p)) return;
        addToQ(tree.getChildren(p).get(0));
        addToQ(tree.getChildren(p).get(1));
        twoQBL(p,tree.getChildren(p).get(0));
        twoQBL(p,tree.getChildren(p).get(1));
    }

    private void twoQBL(jebl.evolution.graphs.Node p, jebl.evolution.graphs.Node q){
        double temp = phenotypes.get(p) - phenotypes.get(q);
        qsum += temp*temp / tree.getLength(q);
    }

    private int numInNodes;
    /**
     * internal node array
     */
    private jebl.evolution.graphs.Node[] internalNodes;
    /**
     * could use the above numbers, but figured this method was more stable
     */
    private Map<jebl.evolution.graphs.Node, Integer> internalNodeNums;
    private double qsum;
    /**
     * all the working phenotypes
     */
    private double [][] fullM;
    /**
     * full matrix copy used in q calculation
     */
    private double [][] fullMcp;
    /**
     * full vector copy
     */
    private double [] fullVcp;
    /**
     * saved full_mat which will be used for calculating the SE
     */
    private double [][] SEfullM;
    /**
     * final phenotypes
     */
    private Map<jebl.evolution.graphs.Node, Double> phenotypes;
    /**
     * final phenotype standard errors
     */
    private Map<jebl.evolution.graphs.Node, Double> phenotypesSE;
    /**
     * tree
     */
    private jebl.evolution.trees.SimpleRootedTree tree;
    /**
     * character
     */
    private jebl.evolution.characters.ContinuousCharacter conchar;
}
