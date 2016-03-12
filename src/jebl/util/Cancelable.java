package jebl.util;

/**
 * Checks whether some operation has been cancelled; This will typically
 * be used as a callback because it can only query whether a task has
 * been cancelled, rather than cancel it (there is no cancel() method).
 *
 * (The name Cancelable kinda seems to imply that there is only a
 * void cancel() method, so maybe this interface should really have a
 * different name - any ideas?)
 *
 * @author Tobias Thierer
 * @version $Id: Cancelable.java 755 2007-08-10 04:47:14Z twobeers $
 *          <p/>
 *          Created on 10/08/2007 16:25:35
 */
public interface Cancelable {
    boolean isCanceled();
}
