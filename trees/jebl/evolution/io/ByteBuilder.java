package jebl.evolution.io;

import java.io.IOException;
import java.util.logging.Logger;

/**
 * Similar to a StringBuilder, but its internal buffer is a byte[] with
 * one entry for each character, so it can only correctly append single-byte
 * characters.
 * 
 * @author Joseph Heled
 * @version $Id: ByteBuilder.java 879 2008-01-31 04:23:43Z twobeers $
 */
public class ByteBuilder implements CharSequence, Appendable {
    private static final Logger logger = Logger.getLogger(ByteBuilder.class.getName());
    final int maxCapacity;
    int current = 0;
    byte[] data;

    void ensureCapacity(final int cap) {
        if (cap > maxCapacity) {
            throw new IllegalArgumentException("requested capacity " + cap + " is > the allowed maximum capacity " + maxCapacity + " for this ByteBuilder");
        }
        if( cap > data.length ) {
            int newLen = 2*(cap+1);
            if( newLen <= 0 ) {
                newLen += 256;
            }
            if( newLen > maxCapacity) newLen = maxCapacity;
            if( newLen < cap ) newLen = cap;
            byte[] d  = new byte[newLen];
            System.arraycopy(data, 0, d, 0, data.length);
            data = d;
        }
    }

    /**
     * Constructs a ByteBuilder that will never grow beyond <code>maxCapacity</code>
     * bytes in length. If you don't want to limit the size this ByteBuilder can
     * grow to, you should pass in Integer.MAX_VALUE here
     * @param maxCapacity The maximum, NOT the initial capacity of this ByteBuilder
     */
    public ByteBuilder(int maxCapacity) {
        this.maxCapacity = maxCapacity;
        if (maxCapacity < 0) {
            throw new IllegalArgumentException("maxCapacity must be positive, but got " + maxCapacity);
        }
        data = new byte[16];
        if (maxCapacity < data.length) {
            logger.warning("ByteBuilder with a very small maxCapacity constructed: " + maxCapacity);
        }
    }

    public static boolean isCharacterAscii(final char c) {
        return ((int) c) < 128;
    }

    /**
     * Appends an ASCII character (see {@link #isCharacterAscii(char)}) to this ByteBuilder.
     * @param c ASCII character to append
     * @return this ByteBuilder
     * @throws IllegalArgumentException if c is not an ASCII character
     */
    public ByteBuilder append(char c) {
        if( current + 1 > data.length ) {
            ensureCapacity(current + 1); // will throw an exception if insufficient capacity (maxCapacity reached)
        }
        if (!isCharacterAscii(c)) {
            throw new IllegalArgumentException("Can't append multi-byte character #" + (int) c + " to a ByteBuilder: " + c);
        }
        data[current] = (byte)c;
        ++current;
        return this;
    }

    /**
     * Appends an ASCII CharSequence to this ByteBuilder.
     * @param charSequence ASCII CharSequence to append
     * @return this ByteBuilder
     * @throws IllegalArgumentException if charSequence contains non-ASCII characters
     */
    public ByteBuilder append(CharSequence charSequence) throws IOException {
        return append(charSequence, 0, charSequence.length());
    }

    public ByteBuilder append(CharSequence csq, int start, int end) throws IOException {
        for (int i = start; i < end; i++) {
            append(csq.charAt(i));
        }
        return this;
    }

    public int length() {
        return current;
    }

    public char charAt(int index) {
        return (char)data[index];
    }

    public CharSequence subSequence(int start, int end) {
        return new String(data, start, end - start);
    }

    public String toString() {
        return new String(data, 0, current);
    }
}
