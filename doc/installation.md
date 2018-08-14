## Installation

### Check prerequisites

 - Java v1.7+ (1.8+ recommended)
 - 2Gb RAM (8Gb recommended)
 - Quad-core 1GHz clock (2.5Ghz recommended)
 
### Getting the latest version

The latest version can be downloaded from the [GitHub](https://github.com/jnarag/teaspoon/releases) repository.

### Installation (general)

Installation is easy. Unzip the installation file. Open a command-line terminal.

### Installation (OSX)

An installer will be made available. Note that depending on your OSX settings/version you may be asked for confirmation before opening TEASPOON for the first time.

### Test the example data

To verify your system has TEASPOON installed correctly, run the example data. Navigate to the TEASPOON folder, and run:

```java -jar TEASPOON.jar example_data/example.fasta example_data/outgroup.fasta example_data/example.mask```

### Note on RAM requirements

TEASPOON is fast but quite memory hungry because sequence alignments are loaded in a fairly simplistic way. This means you may need to allocate more RAM than usual to the Java virtual machine (JVM) which TEASPOON actually runs inside.

The default and maximum JVM sizes vary from system to system, but while 2Gb RAM may be adequate for small (hundreds of sequences/nucleotides) datasets, bigger datasets (thousands of sequences/positions) may well need 10Gb, 16Gb or more. If TEASPOON *does* run out of memory you'll usually get a `java.io.OutOfMemoryException`-type error letting you know.

It's easiest to specify the initial, minimum, and maximum memory directly. In brief see below (but also see [full notes from Oracle/Sun here](http://www.oracle.com/technetwork/java/javase/tech/vmoptions-jsp-140102.html)

```java -Xms16G -Xmx16G -jar TEASPOON.jar example_data/example.fasta example_data/outgroup.fasta example_data/example.mask```

In this example, the `-Xms` and `-Xms` flags set the initial and maximum RAM allocated to the JVM; the `16G` specifies the amount (16Gb). Note that these are positional arguments, e.g. they have to come before the `-jar` argument, and all the following arguments for TEASPOON itself.
