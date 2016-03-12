/**
 * Taxon.java
 *
 * (c) 2005 JEBL Development Team
 *
 * This package is distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package jebl.evolution.taxa;

import jebl.util.AttributableHelper;

import java.lang.ref.WeakReference;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: Taxon.java 865 2007-12-21 03:59:11Z matt_kearse $
 */
public final class Taxon implements jebl.util.Attributable, Comparable {

    /**
     * A private constructor. Taxon objects can only be created by the static Taxon.getTaxon()
     * factory method.
     * @param name the name of the taxon
     */
    public Taxon(String name) {
        this(name, null);
    }


    /**
     * A private constructor. Taxon objects can only be created by the static Taxon.getTaxon()
     * factory method.
     * @param name the name of the taxon
     */
    private Taxon(String name, TaxonomicLevel taxonomicLevel) {
        this.name = name.intern(); // Intern to make the String isn't a substring
        // which references the original string's characters 
        this.taxonomicLevel = taxonomicLevel;
    }

    /**
     * get the name of the taxon
     * @return the name
     */
    public String getName() {
        return name;
    }


    /**
     * get the taxonomic level of the taxon
     * @return the taxonomic level
     */
    public TaxonomicLevel getTaxonomicLevel() {
        return taxonomicLevel;
    }

	// Attributable IMPLEMENTATION

	public void setAttribute(String name, Object value) {
		if (helper == null) {
			helper = new AttributableHelper();
		}
		helper.setAttribute(name, value);
	}

	public Object getAttribute(String name) {
		if (helper == null) {
			return null;
		}
		return helper.getAttribute(name);
	}

    public void removeAttribute(String name) {
        if( helper != null ) {
            helper.removeAttribute(name);
        }
    }

    public Set<String> getAttributeNames() {
        if (helper == null) {
            return Collections.emptySet();
        }
        return helper.getAttributeNames();
    }

	public Map<String, Object> getAttributeMap() {
		if (helper == null) {
			return Collections.emptyMap();
		}
		return helper.getAttributeMap();
	}

	private AttributableHelper helper = null;

    private static void purgeGarbageCollectedTaxa() {
        synchronized(taxa) {
            List<String> taxonNamesToPurge = new ArrayList<String>();
            for (Map.Entry<String,WeakReference<Taxon>> entry : taxa.entrySet()) {
                Taxon taxon = entry.getValue().get();
                if (taxon == null) {
                    taxonNamesToPurge.add(entry.getKey());
                }
            }
            for (String taxonName : taxonNamesToPurge) {
                taxa.remove(taxonName);
            }
            taxaCreatedSinceLastPurge.set(0);
            lastPurgeTime.set(System.currentTimeMillis());
        }
    }

    // Static factory methods

    /**
     * @return a Set containing all the currently created Taxon objects.
     */
    public static Set<Taxon> getAllTaxa() {
        Set<Taxon> result = new HashSet<Taxon>();
        synchronized (taxa) {
            purgeGarbageCollectedTaxa();
            for (Map.Entry<String,WeakReference<Taxon>> entry : taxa.entrySet()) {
                Taxon taxon = entry.getValue().get();
                if (taxon != null) { // might have been garbage collected just now
                    result.add(taxon);
                }
            }
        }
        return Collections.unmodifiableSet(result);
    }

    private static AtomicInteger taxaCreatedSinceLastPurge = new AtomicInteger(0);
    private static AtomicLong lastPurgeTime = new AtomicLong(0);

    /**
     * A static method that returns a Taxon object with the given name. If this has
     * already been created then the same instance will be returned.
     *
     * @param name
     * @return the taxon
     */
    public static Taxon getTaxon(String name) {
        if (name == null) {
            throw new IllegalArgumentException("Illegal null string for taxon name");
        }
        if (name.length() == 0) {
            throw new IllegalArgumentException("Illegal empty string for taxon name");
        }
        Taxon taxon;
        synchronized(taxa) {
            // TT: Can we somehow get around having to manually purge WeakReferences whose values have been
            // garbage collected? E.g. using a WeakHashMap<Taxon,Taxon>() which maps each taxon to is canonicalized form
            // and could be accessed as map.get(new Taxon(...))? Or would a WeakHashMap entry where key==value not get
            // garbage collected?
            if (taxaCreatedSinceLastPurge.get() > 10000 || (lastPurgeTime.get() + 10000 < System.currentTimeMillis())) {
                purgeGarbageCollectedTaxa();
            }
            WeakReference<Taxon> taxonReference = taxa.get(name);
            taxon = (taxonReference == null ? null : taxonReference.get());
            if (taxon == null) {
                taxaCreatedSinceLastPurge.incrementAndGet();
                taxon = new Taxon(name);
                taxa.put(taxon.getName(), new WeakReference<Taxon>(taxon));
            }
        }
        return taxon;
    }

	// private members


    /**
     * The name of this taxon.
     */
    private final String name;

    /**
     * Maps taxon name to the taxon of that name. We use WeakReferences so that if no other references
     * to a taxon exist, it no longer wastes memoryt.
     */
    private static final Map<String, WeakReference<Taxon>> taxa = new HashMap<String, WeakReference<Taxon>>();

    /**
     * the taxonomic level of this taxon.
     */
    private final TaxonomicLevel taxonomicLevel;

    public String toString() {
        return name;
    }

	public int compareTo(Object o) {
		return name.compareTo(((Taxon)o).getName());
	}


    public boolean equals(Taxon t) {
	    return name.equals(t.getName());
    }

    public int hashCode() {
        return name.hashCode();
    }
}
