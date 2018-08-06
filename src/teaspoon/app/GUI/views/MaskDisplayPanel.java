/**
 * 
 */
package teaspoon.app.GUI.views;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.JPanel;

import teaspoon.app.TeaspoonMask;
import teaspoon.app.GUI.models.TeaspoonMaskModel;

/**
 * <b>TEASPOON:<b>
 * <i>Tools for Evolutionary Analysis of Serially-sampled POpulatiONs</i>
 * Jayna Raghwani, Samir Bhatt, Joe Parker &amp; Oliver G. Pybus
 * University of Oxford, 2010-2018.
 * 
 * @author <a href="http://github.com/lonelyjoeparker">@lonelyjoeparker</a>
 * @since 11 Jul 2018
 * @version 0.1
 */
public class MaskDisplayPanel extends JPanel {

	TeaspoonMaskModel model = null;

	private void doMaskDrawing(Graphics g) {

		Graphics2D g2d = (Graphics2D) g;

		g2d.setColor(new Color(212, 212, 212));
		int xpos, ypos, width, height, maskCount, maxLength;
		xpos = 0;
		ypos = 0;
		width = 90;
		height = 60;
		
		try{
			maskCount = model.getMasks().size();
			TeaspoonMask firstMask = model.getMasks().get(0);
			maxLength = firstMask.getLength();
		}catch (NullPointerException ex){
			maskCount = 0;
			maxLength = width;
			ex.printStackTrace();
		}
		
		// draw ruler
		height = 10;
		width = maxLength;
		g2d.drawRect(xpos, ypos, width, height);
		g2d.setColor(new Color(0, 0, 0));
		g2d.fillRect(xpos, ypos, width, height);

		// draw ticks
		g2d.setColor(new Color(255, 0, 0));
		int remainder = maxLength;
		int drawTickPos = 0;
		while(remainder > 10){
			// test if multuple of 100.
			int positionRemainder = drawTickPos % 100;
			if((positionRemainder) == 0){
				// draw a labelled tick
				g2d.drawLine(drawTickPos, ypos, drawTickPos, ypos+20);
				g2d.setColor(new Color(0, 0, 0));
				g2d.drawChars((drawTickPos+"").toCharArray(), (0), ((drawTickPos+"").toCharArray().length), drawTickPos+2, ypos+20);	
				g2d.setColor(new Color(255, 0, 0));
			}else{
				// draw a regular tick
				g2d.drawLine(drawTickPos, ypos, drawTickPos, ypos+10);
			}
			// increment position and recalculate remainder
			drawTickPos += 10;
			remainder = maxLength - (drawTickPos - 10);
		}
		
		// draw masks
		xpos = 10;
		ypos = 25;
		width = 90;
		height = 10;
		g2d.setColor(new Color(125, 167, 116));
		
		for(int i = 0; i<maskCount; i++){
			TeaspoonMask mask = model.getMasks().get(i);
			ArrayList<int[]> lineCoords = this.maskToRanges(mask);
			for(int[] line:lineCoords){
				xpos = line[0];
				width = line[1] - line[0];
				g2d.drawRect(xpos, ypos, width, height);
				g2d.setColor(new Color(125, 167, 116));
				g2d.fillRect(xpos, ypos, width, height);
			}
			/*
			xpos = 10+mask.getFirstStart();
			width = mask.getLastEnd() - mask.getFirstStart();
			*/
			
			// now name
			g2d.setColor(Color.BLACK);			
			// draw the whole name for now
			g2d.drawChars(mask.toString().toCharArray(), (0), (mask.toString().toCharArray().length), xpos, ypos+10);
			
			ypos += 15;
		}
		
		this.setPreferredSize(new Dimension(maxLength+10,ypos+height));
	}

	/**
	 * Utility method. Takes an arraylist and walks through it's positions
	 * pulling out one line for every contiguous run of mask sites
	 * @param mask
	 * @return
	 */
	private ArrayList<int[]> maskToRanges(TeaspoonMask mask){
		ArrayList<int[]> lines = new ArrayList<int[]>();
		boolean[] positions = mask.getPositions();
		boolean lastState = false;
		int[] currentLine = new int[2];
		for(int pos=0;pos<positions.length; pos++){
			if(positions[pos] && lastState){
				// this is true, last was true - do nothing
			}
			if(!positions[pos] && !lastState){
				// this is false, last was false - do nothing
			}
			if(positions[pos] && !lastState){
				// this is true, last was false - start new line
				currentLine[0] = pos;
			}
			if(!positions[pos] && lastState){
				// this is false, last was true - end new line
				currentLine[1] = pos;
				lines.add(currentLine);
				currentLine = new int[2];
			}
			lastState = positions[pos];
		}
		// tidy up the last position
		if(lastState && positions[positions.length-1]){
			currentLine[1] = positions.length-1;
			lines.add(currentLine);
		}
		return lines;
	}

	@Override
	public void paintComponent(Graphics g) {

		super.paintComponent(g);
		if(model != null){
			doMaskDrawing(g);
		}else{
			doDemoDrawing(g);
		}
	}

	public void setModel(TeaspoonMaskModel maskData){
		this.model = maskData;
	}

	private void doDemoDrawing(Graphics g) {
	
		Graphics2D g2d = (Graphics2D) g;
	
		g2d.setColor(new Color(212, 212, 212));
		int xpos, ypos, width, height;
		xpos = 10;
		ypos = 15;
		width = 90;
		height = 60;
		String initText = "Mask list empty. To see mask display: load alignments, set ancestral alignment, and create a new mask.";
		g2d.drawChars(initText.toCharArray(), (0), (initText.toCharArray().length), xpos, ypos+10);

		/*
		g2d.drawRect(xpos, ypos, width, height);
		g2d.drawRect(130, 15, 90, 60);
		g2d.drawRect(250, 15, 90, 60);
		g2d.drawRect(10, 105, 90, 60);
		g2d.drawRect(130, 105, 90, 60);
		g2d.drawRect(250, 105, 90, 60);
		g2d.drawRect(10, 195, 90, 60);
		g2d.drawRect(130, 195, 90, 60);
		g2d.drawRect(250, 195, 90, 60);
	
		g2d.setColor(new Color(125, 167, 116));
		g2d.fillRect(xpos, ypos, width, height);
	
		g2d.setColor(new Color(42, 179, 231));
		g2d.fillRect(130, 15, 90, 60);
	
		g2d.setColor(new Color(70, 67, 123));
		g2d.fillRect(250, 15, 90, 60);
	
		g2d.setColor(new Color(130, 100, 84));
		g2d.fillRect(10, 105, 90, 60);
	
		g2d.setColor(new Color(252, 211, 61));
		g2d.fillRect(130, 105, 90, 60);
	
		g2d.setColor(new Color(241, 98, 69));
		g2d.fillRect(250, 105, 90, 60);
	
		g2d.setColor(new Color(217, 146, 54));
		g2d.fillRect(10, 195, 90, 60);
	
		g2d.setColor(new Color(63, 121, 186));
		g2d.fillRect(130, 195, 90, 60);
	
		g2d.setColor(new Color(31, 21, 1));
		g2d.fillRect(250, 195, 90, 60);
		*/
	}

}
