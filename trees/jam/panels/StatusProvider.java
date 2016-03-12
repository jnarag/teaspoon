package virion.jam.panels;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author rambaut
 *         Date: Jul 27, 2004
 *         Time: 9:32:24 AM
 */
public interface StatusProvider {

    /**
     * Status providers must be able to store a list of StatusListeners. They should
     * then call the appropriate methods on all of these to update the status.
     *
     * @param statusListener the StatusListener to be added
     */
    void addStatusListener(StatusListener statusListener);

    /**
     * Remove the given StatusListener from the provider's list.
     *
     * @param statusListener the StatusListener to be removed
     */
    void removeStatusListener(StatusListener statusListener);

    void fireStatusChanged(int status, String statusText);

    /**
     * The status bar has been pressed. This method should not really be
     * here. You should instead call {@link #fireStatusButtonPressed()}
     */
    void statusButtonPressed();

    String getStatusText();
    int getStatus();

    /**
     * Fire a status bar event to anything interested.
     */
    void fireStatusButtonPressed();

    public void addOverrideProvider(StatusProvider provider);
    public void removeOverrideProvider(StatusProvider provider);


    public class Helper implements StatusProvider {
        private int lastStatus = StatusPanel.NORMAL;
        private String lastStatusText = "";
        // Modifications to this collection are guarded by a lock on this Object,
        // but the listeners shouldn't be called when holding any lock.
        private List listeners = new ArrayList();

        public void addStatusListener(StatusListener statusListener) {
            synchronized (this) {
                listeners.add(statusListener);
            }
            statusListener.statusChanged(lastStatus, lastStatusText);
        }

        public synchronized void removeStatusListener(StatusListener statusListener) {
            synchronized (this) {
                listeners.remove(statusListener);
            }
            statusListener.statusChanged(StatusPanel.NORMAL, "");
        }

        /**
         * This method should not be called when holding a lock on this Object,
         * because it invokes listeners which can do arbitrary things.
         */
        private void realFireStatusChanged() {
            int status;
            String statusText;
            StatusProvider override = null;
            ArrayList listeners;
            synchronized (this) {
                status = lastStatus;
                statusText = lastStatusText;
                int index = overrideProviders.size()- 1;
                if (index >= 0) {
                    override =(StatusProvider) overrideProviders.get(index);
                }
                listeners = new ArrayList(this.listeners);
            }
            if (override != null) {
                status = override.getStatus();
                statusText = override.getStatusText();
            }
            Iterator i = listeners.iterator();
            while (i.hasNext()) {
                ((StatusListener) i.next()).statusChanged(status, statusText);
            }

        }

        public void fireStatusChanged(int status, String statusText) {
            synchronized (this) {
                lastStatus = status;
                lastStatusText = statusText;
            }
            realFireStatusChanged();
        }

        // Modifications to these lists are guarded by a lock on this Object,
        // but the listeners shouldn't be called when holding any lock.
        private List overrideListeners = new ArrayList();
        private List overrideProviders = new ArrayList();

        public void addOverrideProvider(final StatusProvider provider) {
            StatusListener listener = new StatusListener(){
                public void statusChanged(int status, String statusText) {
                    boolean needToFireStatusChange;
                    synchronized(StatusProvider.Helper.this) {
                        needToFireStatusChange = (overrideProviders.size() > 0 && overrideProviders.get(overrideProviders.size() - 1)== provider);
                    }
                    if (needToFireStatusChange) {
                        realFireStatusChanged();
                    }
                }
            };
            synchronized (this) {
                provider.addStatusListener(listener);
                overrideProviders.add(provider);
                overrideListeners.add(listener);
            }
        }

        public void removeOverrideProvider(StatusProvider provider) {
            synchronized (this) {
                for (int i = 0; i < overrideProviders.size(); i++) {
                    if(overrideProviders.get(i) == provider) {
                        overrideProviders.remove(i);
                        StatusListener listener =(StatusListener) overrideListeners.get(i);
                        provider.removeStatusListener(listener);
                        overrideListeners.remove(i);
                        break;
                    }
                }
            }
            realFireStatusChanged();
        }

        public void fireStatusButtonPressed(){
            StatusProvider override = null;
            synchronized (this) {
                int index = overrideProviders.size ();
                if (index > 0) {
                    override = ((StatusProvider)overrideProviders.get(index - 1));
                }
            }

            if (override != null) {
                override.fireStatusButtonPressed();
            } else {
                statusButtonPressed();
            }
        }

        public void statusButtonPressed() {
            // do nothing... override
        }

        public synchronized int getStatus() {
            return lastStatus;
        }

        public synchronized String getStatusText() {
            return lastStatusText;
        }

        public synchronized boolean hasOverrideProvider() {
            return !overrideProviders.isEmpty();
        }
    }
}
