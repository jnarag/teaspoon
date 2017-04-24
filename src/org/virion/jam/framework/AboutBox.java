package org.virion.jam.framework;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.StringTokenizer;

import javax.swing.BorderFactory;
import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JEditorPane;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.virion.jam.html.SimpleLinkListener;
import org.virion.jam.util.IconUtils;
import org.virion.jam.util.Utils;

public class AboutBox extends AbstractFrame {

    /**
     * Creates an AboutBox with a given title, message and icon
     * and centers it over the parent component.
     */
    public AboutBox(String title, String message, Icon icon) {
        super();

        if (icon != null) {
            setIconImage(IconUtils.getBufferedImageFromIcon(icon));
        }

        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent event) {
                close();
            }
        });
        addKeyListener(new KeyAdapter() {
            public void keyPressed(KeyEvent event) {
                if (event.getKeyCode() == KeyEvent.VK_ESCAPE
                        ||
                        event.getKeyCode() == KeyEvent.VK_ENTER) {
                    close();
                }
            }
        });

        JPanel contentsPanel = new JPanel(new GridBagLayout());
        contentsPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        contentsPanel.setBackground(Color.white);

        JLabel iconLabel = new JLabel(icon, JLabel.CENTER);
        JLabel titleLabel = new JLabel(title, JLabel.CENTER);

        Font font = titleLabel.getFont();
        titleLabel.setFont(font.deriveFont(16.0f).deriveFont(Font.BOLD));

        GridBagConstraints c = new GridBagConstraints();
        c.gridwidth = GridBagConstraints.REMAINDER;
        c.insets = new Insets(5,5,5,5);
        contentsPanel.add(iconLabel, c);
        contentsPanel.add(titleLabel, c);

        font = font.deriveFont(11.0f);

        if (message.startsWith("<html>")) {
            JEditorPane editorPane = new JEditorPane("text/html", message);
            editorPane.setOpaque(false);
            editorPane.setFont(font);
            editorPane.setEditable(false);
            editorPane.addHyperlinkListener(new SimpleLinkListener());
            contentsPanel.add(editorPane, c);
        } else {
            StringTokenizer tokens = new StringTokenizer(message, "\n");
            while (tokens.hasMoreElements()) {
                String text = tokens.nextToken();
                JLabel messageLabel = new JLabel(text, JLabel.CENTER);
                messageLabel.setFont(font);
                contentsPanel.add(messageLabel, c);
            }
        }


        getSaveAction().setEnabled(false);
        getSaveAsAction().setEnabled(false);
        getPrintAction().setEnabled(false);
        getPageSetupAction().setEnabled(false);

        getCutAction().setEnabled(false);
        getCopyAction().setEnabled(false);
        getPasteAction().setEnabled(false);
        getDeleteAction().setEnabled(false);
        getSelectAllAction().setEnabled(false);
        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);
        getMinimizeWindowAction().setEnabled(true);
        getCloseWindowAction().setEnabled(true);

        getContentPane().add(contentsPanel);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setResizable(false);
        pack();
        Utils.centerComponent(this, null);
    }

    /**
     * Sets the visibility to false and disposes the frame.
     */
    public void close() {
        setVisible(false);
        dispose();
    }

    protected void initializeComponents() {
    }

    public boolean requestClose() {
        return false;
    }

    public JComponent getExportableComponent() {
        return null;
    }
}


