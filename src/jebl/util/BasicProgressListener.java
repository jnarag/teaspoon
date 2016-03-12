package jebl.util;

/**
 * A basic {@link jebl.util.ProgressListener} implementation that allows the caller to set the canceled
 * status. Typically this class is used when you don't care about providing immediate feedback to
 * the user, but still need to be able to cancel an operation.
 * @author Matt Kearse
 * @version $Id$
 */

public class BasicProgressListener extends ProgressListener {
    private volatile boolean canceled=false;
    private volatile String message="";
    private volatile double fractionCompleted = 0;
    private volatile boolean indeterminate = false;
    protected void _setProgress(double fractionCompleted) {
        this.indeterminate = false;
        this.fractionCompleted = fractionCompleted;
    }

    protected void _setIndeterminateProgress() {
        this.indeterminate = true;
    }

    protected void _setMessage(String message) {
        this.message = message;
    }

    public boolean isCanceled() {
        return canceled;
    }

    /**
     * Sets this progress listener as cancel, so that {@link #isCanceled()} will return true.
     */
    public void cancel() {
        this.canceled = true;
    }

    /**
     *
     * @return the most recent message set on this progress listener.
     */
    public String getMessage() {
        return message;
    }

    /**
     * @return the current fraction complete.
     */
    public double getFractionCompleted() {
        return fractionCompleted;
    }

    /**
     * @return true if {@link #setIndeterminateProgress()} has been called, and {@link #setProgress(double)} has not been called since
     */
    public boolean isIndeterminate() {
        return indeterminate;
    }
}
