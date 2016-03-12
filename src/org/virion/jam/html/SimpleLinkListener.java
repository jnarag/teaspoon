package org.virion.jam.html;


import org.virion.jam.util.BrowserLauncher;

import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

/**
 * iSeek prototype. Codename seekquence.
 *
 * This class listens to Hyperlink Events, and opens the url in a browser window.
 *
 * Open a browser from a Java application on Windows, Unix, or Macintosh.
 * see  http://ostermiller.org/utils/Browser.html  for more information
 *
 * @author Nasser Giacaman
 * @version $Id: SimpleLinkListener.java 844 2007-11-21 23:17:11Z twobeers $
 *          Date: 26/01/2005
 *          Time: 11:54:50
 */
public class SimpleLinkListener implements HyperlinkListener {

    public void hyperlinkUpdate(HyperlinkEvent he) {

        if (he.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
            try{
                BrowserLauncher.openURL(he.getDescription());
            }catch(Exception ioe){
                ioe.printStackTrace();
            }
        }
    }
}
