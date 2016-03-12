package jebl.gui.trees.treeviewer.decorators;

import jebl.evolution.trees.Tree;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.graphs.Node;

import java.awt.*;

/**
 * @author Steven Stones-Havas
 * @version $Id$
 *          <p/>
 *          Created on 15/07/2008 13:19:33
 */
public class RgbBranchDecorator implements BranchDecorator{
    public Paint getBranchPaint(Tree tree, Node node) {
        if(!(tree instanceof RootedTree)){
            assert false;
            return Color.black;
        }
        RootedTree rootedTree = (RootedTree)tree;
        Node newNode = rootedTree.getParent(node);
        if(newNode == null)
            return Color.black;
        Object colorString = newNode.getAttribute("nodeColor");
        if(colorString == null)
            return Color.black;
        String rgb = colorString.toString();
        return getColorFromString(rgb);
    }

    public static Color getColorFromString(String rgb) {
        String[] rgbStrings = rgb.split(",");
        if(rgbStrings.length != 3){
            assert false;
            return Color.black;
        }
        try{
            int r = Integer.parseInt(rgbStrings[0]);
            int g = Integer.parseInt(rgbStrings[1]);
            int b = Integer.parseInt(rgbStrings[2]);
            return new Color(r,g,b);
        }
        catch(NumberFormatException ex){
            ex.printStackTrace();
            assert false : ex.getMessage();
        }
        return Color.black;
    }
}
