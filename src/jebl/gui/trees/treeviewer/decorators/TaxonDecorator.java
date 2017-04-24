package jebl.gui.trees.treeviewer.decorators;

import java.awt.Font;
import java.awt.Paint;

import jebl.evolution.taxa.Taxon;

/**
 * @author Andrew Rambaut
 * @version $Id: TaxonDecorator.java 181 2006-01-23 17:31:10Z rambaut $
 */
public interface TaxonDecorator {
    Paint getTaxonPaint(Taxon taxon);
    Font getTaxonFont(Taxon taxon, Font font);
}
