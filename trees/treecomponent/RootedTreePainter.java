/*
 * RootedTreePainter.java
 *
 * (c) 2005 JEBL Development Team
 *
 * This package is distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package jebl.gui.trees.treecomponent;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Stroke;
import java.awt.geom.Point2D;

/**
 * @author Alexei Drummond
 *
 * @version $Id: RootedTreePainter.java 309 2006-05-02 08:22:06Z rambaut $
 */
public interface RootedTreePainter  {

	/**
	*	Set line style
	*/
	void setLineStyle(Stroke lineStroke, Paint linePaint);

	/**
	*	Set hilight style
	*/
	void setHilightStyle(Stroke hilightStroke, Paint hilightPaint);

	/**
	 *	Set label style.
	 */
	void setLabelStyle(Font labelFont, Paint labelPaint);

	/**
	 *	Set hilight label style.
	 */
	void setHilightLabelStyle(Font hilightLabelFont, Paint hilightLabelPaint);

	/**
	 * Do the actual painting.
	 */
	void paintTree(Graphics2D g, Dimension size, jebl.evolution.trees.RootedTree tree);

	/**
	*	Find the node under point. Returns -1 if not found.
	*/
	public jebl.evolution.graphs.Node findNodeAtPoint(Point2D point);

}
