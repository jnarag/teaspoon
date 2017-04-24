package org.virion.jam.html;

import java.awt.BorderLayout;

import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JScrollPane;

/**
 * General-purpose class to display HTML in a standalone frame.
 */
public class HTMLViewer extends JFrame {

    private JEditorPane editorPane;

    public HTMLViewer(String title, String html) {
        super(title);
        getContentPane().setLayout(new BorderLayout());
        editorPane = new JEditorPane("text/html", html);
        editorPane.setEditable(false);
        getContentPane().add(new JScrollPane(editorPane), BorderLayout.CENTER);
    }
}





