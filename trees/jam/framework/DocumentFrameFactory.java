/**
 * DocumentFrameFactory.java
 */

package virion.jam.framework;



public interface DocumentFrameFactory {

    DocumentFrame createDocumentFrame(Application app, MenuBarFactory menuBarFactory);
}