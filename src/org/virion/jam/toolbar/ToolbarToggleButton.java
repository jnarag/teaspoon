package org.virion.jam.toolbar;

import javax.swing.*;

/**
 * @author rambaut
 *         Date: Oct 18, 2005
 *         Time: 10:09:21 PM
 */
public class ToolbarToggleButton extends JButton implements ToolbarItem {

	public ToolbarToggleButton(ToolbarAction action, boolean doubleClickToToggle) {
		super(action);

		this.doubleClickToToggle = doubleClickToToggle;

		setHorizontalTextPosition(SwingConstants.CENTER);
		setVerticalTextPosition(SwingConstants.BOTTOM);
		putClientProperty("Quaqua.Button.style", "toolBarTab");
		putClientProperty("JButton.buttonType", "toolbar");
		setBorderPainted(true);

	    setToolTipText(action.getToolTipText());

		setDisabledIcon(action.getDisabledIcon());
		setPressedIcon(action.getPressedIcon());
		setSelectedIcon(action.getPressedIcon());
	}

	public void setToolbarOptions(ToolbarOptions options) {
		switch (options.getDisplay()) {
			case ToolbarOptions.ICON_AND_TEXT:
				setText(action.getLabel());
				setIcon(action.getIcon());
				setDisabledIcon(action.getDisabledIcon());
				setPressedIcon(action.getPressedIcon());
				break;
			case ToolbarOptions.ICON_ONLY:
				setText(null);
				setIcon(action.getIcon());
				setDisabledIcon(action.getDisabledIcon());
				setPressedIcon(action.getPressedIcon());
				break;
			case ToolbarOptions.TEXT_ONLY:
				setText(action.getLabel());
				setIcon(null);
				setDisabledIcon(null);
				setPressedIcon(null);
				break;
		}
	}

	public void doClick() {
		super.doClick();
	}

	public void doClick(int i) {
		super.doClick(i);
	}

	public void setAction(Action action) {
		super.setAction(action);
		if (action instanceof ToolbarAction) {
			this.action = (ToolbarAction)action;
		}
	}

	private ToolbarAction action;
	private boolean doubleClickToToggle;
}