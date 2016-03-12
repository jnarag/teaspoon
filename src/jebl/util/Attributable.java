/*
 * Attributable.java
 *
 * (c) 2005 JEBL Development Team
 *
 * This package is distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package jebl.util;

import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Interface for associating attributeNames with an object.
 *
 * @version $Id: Attributable.java 849 2007-12-06 00:10:14Z twobeers $
 *
 * @author Andrew Rambaut
 */
public interface Attributable {

    /**
     * Sets an named attribute for this object.
     * @param name the name of the attribute.
     * @param value the new value of the attribute.
     */
    void setAttribute(String name, Object value);

    /**
     * @return an object representing the named attributed for this object.
     * @param name the name of the attribute of interest, or null if the attribute doesn't exist.
     */
    Object getAttribute(String name);

    /**
     * @param name name of attribute to remove
     */
    void removeAttribute(String name);

    /**
     * @return an array of the attributeNames that this object has.
     */
    Set<String> getAttributeNames();

    /**
     * Gets the entire attribute map.
     * @return an unmodifiable map
     */
    Map<String, Object> getAttributeMap();

    public static class Utils {
        private Utils() { }  // make class uninstantiable

        Set<String> getAttributeNames(Collection<Attributable> attributables) {
            Set<String> names = new HashSet<String>();

            for (Attributable attributable : attributables) {
                names.addAll(attributable.getAttributeNames());
            }

            return names;
        }
    }
}


