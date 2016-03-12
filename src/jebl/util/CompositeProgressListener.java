package jebl.util;

import jebl.math.MachineAccuracy;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * A {@link jebl.util.ProgressListener} that is suitable for a task that consists of several subtasks.
 * You specify the relative duration of each subtask, and then the subtasks' setProgress()
 * calls with values between 0 and 1 are translated to reflect the overall progress on the whole
 * (combined) task. In other words, each subtask reports progress as if it were the whole task,
 * and the CompositeProgressListener translates this into overall progress. This also implies
 * that calling {@link jebl.util.ProgressListener#setComplete()} or {@link ProgressListener#setProgress(double) setProgress(1.0)}
 * marks the current subtask rather than the entire task completed.
 * <p/>
 * As the combined progress listener cannot know which subtask it is currently being called from,
 * you have to explicitely let it know when a new subtask (not the first) starts, by calling
 * {@link #beginNextSubtask()}. Thus when the constructor is passed an array of N doubles as its second
 * argument, {@link #beginNextSubtask()} should be called precisely N-1 times.
 * <p/>
 * Alternatively, instead of calling {@link #beginNextSubtask()} after each subtask (except the last),
 * you can instead call {@link #beginSubtask()} before each subtask (including the first).
 * <p/>
 *
 * @author Tobias Thierer
 * @version $Id: CompositeProgressListener.java 944 2008-09-26 05:08:15Z matthew_cheung $
 */
public final class CompositeProgressListener extends ProgressListener {
    protected int numOperations;
    protected ProgressListener listener;
    protected int currentOperationNum = 0;
    protected double[] time;
    protected double baseTime = 0.0; // overall progress (0..1) at the start of the current sub-operation
    protected double currentOperationProgress = 0.0;
    private boolean beganFirstSubTask=false;
    private String currentSubTaskMessage="";

    /**
     * construct a new composite ProgressListener.
     *
     * @param listener the ProgressListener that all progress reports are forwarded to after adjusting them for the currently active sub-task
     * @param operationDuration a list of relative weightings to give each sub task.
     */
    public CompositeProgressListener(ProgressListener listener, double ... operationDuration) {
        numOperations = operationDuration.length;
        if (numOperations == 0) {
            // Give a slightly more helpful message in this special case (would
            // otherwise be caught by the totalTime != 0.0 test below as well)
            throw new IllegalArgumentException("Composite operation must have > 0 subtasks");
        }
        if (listener == null) {
            this.listener = ProgressListener.EMPTY;
        } else {
            this.listener = listener;
        }
        this.time = operationDuration.clone();

        // scale times to a sum of 1
        double totalTime = 0.0;
        for (double d : operationDuration) {
            if (d < 0.0) {
                throw new IllegalArgumentException("Operation cannot take negative time: " + d);
            }
            totalTime += d;
        }
        if (MachineAccuracy.same(totalTime, 0.0)) { // will always be the case if numOperations == 0
            throw new IllegalArgumentException("There must be at least one subtask that takes > 0 time");
        }
        for (int i = 0; i < numOperations; i++) {
            this.time[i] = (operationDuration[i] / totalTime);
        }
    }

    private static double[] getWeights(int numberOfEvenlyWeightedSubTasks) {
        double[] results=new double[numberOfEvenlyWeightedSubTasks];
        Arrays.fill(results,1);
        return results;
    }

    /**
     * Construct a CompositeProgressListener with a number of evenly weighted subtasks.
     * @param listener the ProgressListener that all progress reports are forwarded to after adjusting them for the currently active sub-task
     * @param numberOfEvenlyWeightedSubTasks the number of evenly weighted sub-tasks.
     */
    public CompositeProgressListener(ProgressListener listener, int numberOfEvenlyWeightedSubTasks) {
        //noinspection PrimitiveArrayArgumentToVariableArgMethod
        this(listener,getWeights(numberOfEvenlyWeightedSubTasks));
    }

    public static CompositeProgressListener forFiles(ProgressListener listener, List<File> files) {
        int n = files.size();
        double[] lengths = new double[n];
        int i =0;
        for (File file : files) {
            lengths[i++] = (double) file.length();
        }
        return new CompositeProgressListener(listener, lengths);
    }

    /**
     * Used as an alternative to {@link #beginNextSubtask()}.
     * Instead of calling {@link #beginNextSubtask()} once after each subtask
     * (except the last), you can instead call beginSubTask at the beginning
     * of every subtask including the first.
     */
    public void beginSubtask() {
        currentSubTaskMessage = "";
        if (!beganFirstSubTask) {
            beganFirstSubTask = true;
        } else {
            beginNextSubtask();
        }
    }

    /**
     * Used as an alternative to {@link #beginNextSubtask()}.
     * Instead of calling {@link #beginNextSubtask()} once after each subtask
     * (except the last), you can instead call beginSubTask at the beginning
     * of every subtask including the first.
     * @param message a message to be displayed to the user as a prefix in the progress message
     */
    public void beginSubtask(String message) {
        currentSubTaskMessage = message;
        setMessage("");
        beginSubtask();
        currentSubTaskMessage = message;
    }

    protected void _setProgress(double fractionCompleted) {
        currentOperationProgress = fractionCompleted;
        listener._setProgress(baseTime + fractionCompleted * time[currentOperationNum]);
    }

    protected void _setIndeterminateProgress() {
        listener._setIndeterminateProgress();
    }

    protected void _setMessage(String message) {
        if (currentSubTaskMessage.length()>0) {
             message=currentSubTaskMessage+(message.length()>0?": "+message:"");
            // concatentate the parent message and the sub-task messages. Previous behaviour was
            // just to overwrite the parent sub-task message, but I think this is just wrong.
        }
        listener._setMessage(message);
    }

    public boolean isCanceled() {
        return listener.isCanceled();
    }

    public boolean addProgress(double fractionCompletedDiff) {
        return setProgress(currentOperationProgress + fractionCompletedDiff);
    }

    /**
     * @return true if there is another subtask available after the current one
     */
    public boolean hasNextSubtask() {
        return (currentOperationNum < (numOperations - 1));
    }

    /**
     * Clear all progress, including that of previous subtasks.
     * Note: if the task has already been canceled, this does not reset its status to non-canceled.
     */
    public void clearAllProgress () {
        currentOperationNum = 0;
        baseTime = 0.0;
        setProgress(0.0);
        //listener.setProgress(0);
    }

    /**
     * Convenience method to start the next operation AND set a new message.
     * @param message message to set (will be passed to setMessage()
     */
    public void beginNextSubtask(String message) {
        beginNextSubtask();
        currentSubTaskMessage=message;
        setMessage("");
    }

    /**
     * begins the next subtask. Should not be called on the first subtask, but should only be called
     * to start tasks after the first one. If you wish to call a begin subtask method
     * for each task including the first, use {@link #beginSubtask()} instead.
     */
    public void beginNextSubtask() {
        currentSubTaskMessage = "";
        setComplete();
        if (!hasNextSubtask()) {
            throw new IllegalStateException(currentOperationNum + " " + numOperations);
        }
        baseTime += time[currentOperationNum];
        currentOperationNum++;
        currentOperationProgress = 0.0;
    }

    public final boolean setComplete() {
        return setProgress(1.0);
    }

    /**
     * @return The root {@link jebl.util.ProgressListener} that this forwards adjusted progress reports to.
     */
    public ProgressListener getRootProgressListener() {
        if(listener instanceof CompositeProgressListener) {
            return ((CompositeProgressListener)listener).getRootProgressListener();
        } else {
            return listener;
        }
    }
}
