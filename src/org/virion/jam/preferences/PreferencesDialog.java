package org.virion.jam.preferences;

import org.virion.jam.toolbar.*;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.*;
import java.util.List;

/**
 * PreferencesDialog.java
 *
 * @author			Andrew Rambaut
 * @version			$Id: PreferencesDialog.java 751 2007-08-06 22:41:46Z rambaut $
 */
public class PreferencesDialog {

    private JFrame frame;
    private CardLayout cardLayout;
    private JPanel sectionsPanel;

    public PreferencesDialog(JFrame frame) {
        this.frame = frame;
    }

    public void showDialog() {

        JPanel panel = new JPanel(new BorderLayout());
        Toolbar toolbar = new Toolbar(null);
        toolbar.setFloatable(false);

        cardLayout = new CardLayout();
        sectionsPanel = new JPanel(cardLayout);
	    sectionsPanel.setBorder(new EmptyBorder(12,12,12,12));

        panel.add(toolbar, BorderLayout.NORTH);
        panel.add(sectionsPanel, BorderLayout.CENTER);

        JOptionPane optionPane = new JOptionPane(panel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.DEFAULT_OPTION,
                null,
                new String[] { "Done" },
                null);
        optionPane.setBorder(new EmptyBorder(0,0,12,0));

        final JDialog dialog = optionPane.createDialog(frame, currentSection);

        for (PreferencesSection section : sections) {
            final String title = section.getTitle();
            if (currentSection == null) {
                currentSection = title;
            }
            final ToolbarButton button = new ToolbarButton(
                    new ToolbarAction(title, title, section.getIcon()) {
                        public void actionPerformed(ActionEvent e) {
                            showSection(title);
                            currentSection = title;
                            dialog.setTitle(currentSection);
                        }
                    }
            );
            JPanel buttonPanel = new JPanel(new BorderLayout());
            buttonPanel.setBorder(BorderFactory.createEmptyBorder(0,1,0,1));
            buttonPanel.add(button, BorderLayout.CENTER);

            toolbar.addComponent(buttonPanel);
            sectionsPanel.add(section.getPanel(), title);
            buttons.put(title, buttonPanel);

            section.retrievePreferences();
        }
        toolbar.addFlexibleSpace();

        showSection(currentSection);

	    for (PreferencesSection section : sections) {
	        section.retrievePreferences();
	    }

        dialog.pack();
        dialog.setVisible(true);
	   
    }

    public void showSection(String title) {
        cardLayout.show(sectionsPanel, title);

        JPanel buttonPanel = buttons.get(currentSection);
        buttonPanel.setBorder(BorderFactory.createEmptyBorder(0,1,0,1));
        buttonPanel.setOpaque(false);

        buttonPanel = buttons.get(title);
        buttonPanel.setBorder(BorderFactory.createMatteBorder(0,1,0,1,Color.gray));
        buttonPanel.setBackground(new Color(0.85F, 0.85F, 0.85F, 0.5F));
        buttonPanel.setOpaque(true);
        buttonPanel.repaint();
    }


    public void addSection(PreferencesSection section) {
        sections.add(section);
    }

    String currentSection = null;

    private List<PreferencesSection> sections = new ArrayList<PreferencesSection>();
    private Map<String, JPanel> buttons = new HashMap<String, JPanel>();
}