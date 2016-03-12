package org.virion.jam.disclosure;

import org.virion.jam.util.IconUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id: DisclosurePanel.java 948 2008-10-02 00:19:47Z rambaut $
 */
public class DisclosurePanel extends JPanel {

	/**
	 *
	 * @param title The title of the panel
	 * @param panel The contents of the panel
	 * @param isOpen Whether the panel should start open
	 */
	public DisclosurePanel(final String title, final JPanel panel, boolean isOpen) {
        this(new JLabel(title), -1, panel, isOpen, 50);
    }

    /**
     *
     * @param titleComponent The title component of the panel
     * @param panel The contents of the panel
     * @param isOpen Whether the panel should start open
     */
    public DisclosurePanel(final JComponent titleComponent, final JPanel panel, boolean isOpen, int openSpeed) {
        this(titleComponent, -1, panel, isOpen, openSpeed);
    }

	/**
	 *
	 * @param titleComponent The component to use as the title of the panel
	 * @param contentPanel The contents of the panel
	 * @param isOpen Whether the panel should start open
	 * @param openSpeed The opening speed in milliseconds
	 */
    public DisclosurePanel(final JComponent titleComponent, int preferredTitleHeight,
                           final JPanel contentPanel,
	                       boolean isOpen, int openSpeed) {

        setOpaque(false);

        panel = contentPanel;
        panel.setOpaque(false);

		button = new DisclosureButton(openSpeed);

        this.titleComponent = titleComponent;
        titleComponent.setOpaque(false);

		JPanel panel1 = new JPanel(new BorderLayout(6, 0)) {
			public void paint(Graphics graphics) {
				graphics.drawImage(backgroundImage, 0, 0, getWidth(), getHeight(), null);
                super.paint(graphics);
            }
		};
        panel1.setOpaque(false);
        panel1.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                button.doClick();
            }
        });

        JPanel componentPanel = new JPanel(new BorderLayout(6, 0));
        componentPanel.setOpaque(false);

        componentPanel.add(button, BorderLayout.WEST);
		componentPanel.add(titleComponent, BorderLayout.CENTER);
        panel1.add(componentPanel, BorderLayout.CENTER);

        if (preferredTitleHeight > 0) {
            panel1.setPreferredSize(new Dimension(0, preferredTitleHeight));
        }
        
        setLayout(new BorderLayout(0, 0));
        add(panel1, BorderLayout.NORTH);
		add(panel, BorderLayout.CENTER);

        setOpen(isOpen);

        button.addDisclosureListener(new DisclosureListener() {
			public void opening(Component component) {
				fireOpening();
			}

			public void opened(Component component) {
				panel.setVisible(true);
				fireOpened();
			}

			public void closing(Component component) {
				fireClosing();
			}

			public void closed(Component component) {
				panel.setVisible(false);
				fireClosed();
			}
		});
    }

    public Dimension getMinimumSize() {
        Dimension size = super.getPreferredSize();
        return size;
    }

    public void setOpen(boolean isOpen) {
		button.setSelected(isOpen);
		panel.setVisible(isOpen);
    }

	public void addDisclosureListener(DisclosureListener listener) {
		listeners.add(listener);
	}

	public void removeDisclosureListener(DisclosureListener listener) {
		listeners.remove(listener);
	}

	private void fireOpening() {
        for (DisclosureListener listener : listeners) {
            (listener).opening(this);
        }
    }

	private void fireOpened() {
        for (DisclosureListener listener : listeners) {
            (listener).opened(this);
        }
    }

	private void fireClosing() {
        for (DisclosureListener listener : listeners) {
            (listener).closing(this);
        }
    }

	private void fireClosed() {
        for (DisclosureListener listener : listeners) {
            (listener).closed(this);
        }
    }

	public DisclosureButton getDisclosureButton() {
		return button;
	}

	public Component getTitleComponent() {
		return titleComponent;
	}

	public JPanel getContentsPanel() {
		return panel;
	}

    private final DisclosureButton button;
    private final Component titleComponent;
    private final JPanel panel;
	private final java.util.List<DisclosureListener> listeners = new ArrayList<DisclosureListener>();

	private static BufferedImage backgroundImage = null;

	static {
		try {
			backgroundImage = IconUtils.getBufferedImage(DisclosurePanel.class, "images/titleBackground.png");

		} catch (Exception e) {
			// no icons...
		}
	}
}
