package jebl.gui.trees.treeviewer.decorators;

import java.awt.Color;
import java.awt.Paint;

import jebl.evolution.trees.Tree;

/**
 * @author Steven Stones-Havas
 * @version $Id$
 *          <p/>
 *          Created on 15/07/2008 13:19:33
 */
public class RgbBranchDecorator implements BranchDecorator {
    public Paint getBranchPaint(Tree tree, jebl.evolution.graphs.Node node) {
        if(!(tree instanceof jebl.evolution.trees.RootedTree)){
            assert false;
            return Color.black;
        }
        jebl.evolution.trees.RootedTree rootedTree = (jebl.evolution.trees.RootedTree)tree;
        jebl.evolution.graphs.Node newNode = rootedTree.getParent(node);
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
