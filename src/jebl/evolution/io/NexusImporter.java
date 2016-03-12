/*
 * NexusImporter.java
 *
 * (c) 2002-2005 JEBL development team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */
package jebl.evolution.io;

import jebl.evolution.alignments.Alignment;
import jebl.evolution.alignments.BasicAlignment;
import jebl.evolution.distances.BasicDistanceMatrix;
import jebl.evolution.distances.DistanceMatrix;
import jebl.evolution.graphs.Node;
import jebl.evolution.sequences.BasicSequence;
import jebl.evolution.sequences.Sequence;
import jebl.evolution.sequences.SequenceType;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.CompactRootedTree;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.trees.SimpleRootedTree;
import jebl.evolution.trees.Tree;
import jebl.util.Attributable;

import java.awt.*;
import java.io.EOFException;
import java.io.IOException;
import java.io.Reader;
import java.util.*;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class for importing NEXUS file format.
 *
 * This is a good starting point for documentation about the nexus file format:
 *
 *    https://www.nescent.org/wg_phyloinformatics/NEXUS_Specification
 *
 * For a quick check if something is probably valid or not, consider this attempt at a
 * nexus grammar (the above page says "Don't treat this as gospel, its just an attempt
 * to get the syntax rules"):
 *
 *    http://www.cs.nmsu.edu/~epontell/nexus/nexus_grammar
 *
 * @version $Id: NexusImporter.java 935 2008-07-22 16:52:04Z rambaut $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class NexusImporter implements AlignmentImporter, SequenceImporter, TreeImporter, DistanceMatrixImporter {

    /**
     * Represents the block types used in the nexus format
     */
    public enum NexusBlock {
		UNKNOWN,
		TAXA,
		CHARACTERS,
		DATA,
		UNALIGNED,
		DISTANCES,
		TREES
	}

    private boolean compactTrees = false;

    // NEXUS specific ImportException classes
    /**
     * Thrown when a block is missing that is required for importing a particular type of data from the nexus input.
     */
    public static class MissingBlockException extends ImportException {
		public MissingBlockException() { super(); }
		public MissingBlockException(String message) { super(message); }
	}

	/**
     * @param reader
     * @param expectedLength Expected length of the input in bytes, or 0 if unknown. Used for optimization and tracking
     *                       progress.
	 */
	public NexusImporter(Reader reader, long expectedLength) {
		helper = new ImportHelper(reader);
        helper.setExpectedInputLength(expectedLength);
        initHelper();
	}

	public NexusImporter(Reader reader) {
		this(reader, 0);
	}

    /**
     *
     * @param reader
     * @param compactTrees true to import trees as {@link jebl.evolution.trees.CompactRootedTree}, false to use
     * {@link jebl.evolution.trees.SimpleRootedTree}.
     * @param expectedInputLength Expected length of the input in bytes, or 0 if unknown. Used for optimization and tracking
     *                            progress.
     */
    public NexusImporter(Reader reader, boolean compactTrees, long expectedInputLength) {
        this(reader, expectedInputLength);
        this.compactTrees = compactTrees;
    }

    /**
     * @param reader
     * @param compactTrees
     * @deprecated Use NexusImporter(Reader reader, boolean compactTrees, long expectedInputLength)
     */
    @Deprecated
    public NexusImporter(Reader reader, boolean compactTrees) {
        // a wild guess on the low side
        this(reader, compactTrees, 4*1024);
    }

    private void initHelper() {
        // ! defines a comment to be written out to a log file
        // & defines a meta comment
        helper.setCommentDelimiters('[', ']', '\0', '!', '&');
    }

    /**
     * Read ahead to the next block in the input.
     *
	 * This should be overridden to provide support for other blocks. Once
	 * the block is read in, nextBlock is automatically set to UNKNOWN by
	 * findEndBlock.
     *
     * @return the type of the next block
	 */
	public NexusBlock findNextBlock() throws IOException
	{
		findToken("BEGIN", true);
		nextBlockName = helper.readToken(";").toUpperCase();
		return findBlockName(nextBlockName);
	}

	/**
	 * This function returns an enum class to specify what the
	 * block given by blockName is.
	 */
	private NexusBlock findBlockName(String blockName)
	{
		try {
			nextBlock = NexusBlock.valueOf(blockName);
		} catch( IllegalArgumentException e ) {
			// handle unknown blocks. java 1.5 throws an exception in valueOf
			nextBlock = null;
		}

		if (nextBlock == null) {
			nextBlock = NexusBlock.UNKNOWN;
		}

		return nextBlock;
	}

    /**
     *
     * @return the name of the next block which will be read. eg "TREES"
     * @see {@link jebl.evolution.io.NexusImporter.NexusBlock}
     */
    public String getNextBlockName() {
		return nextBlockName;
	}

	/**
	 * @return an Iterator over the trees in the nexus input.
	 */
	public Iterator<Tree> iterator() {
		return new Iterator<Tree>() {

			public boolean hasNext() {
				boolean hasNext = false;
				try {
					hasNext = hasTree();
				} catch (IOException e) {
					// deal with errors by stopping the iteration
				} catch (ImportException e) {
					// deal with errors by stopping the iteration
				}
				return hasNext;
			}

			public Tree next() {
				Tree tree = null;
				try {
					tree = importNextTree();
				} catch (IOException e) {
					// deal with errors by stopping the iteration
				} catch (ImportException e) {
					// deal with errors by stopping the iteration
				}
				if (tree == null) throw new NoSuchElementException("No more trees in this file");
				return tree;
			}

			public void remove() {
				throw new UnsupportedOperationException("operation is not supported by this Iterator");
			}
		};
	}

	/**
	 * Parse the next 'TAXA' block encountered in the input.
	 */
	public List<Taxon> parseTaxaBlock() throws ImportException, IOException
	{
		return readTaxaBlock();
	}

	/**
	 * Parse the next 'CHARACTERS' block encountered in the input.
	 */
	public List<Sequence> parseCharactersBlock(List<Taxon> taxonList) throws ImportException, IOException
	{
		return readCharactersBlock(taxonList);
	}

	/**
	 * Parse the next 'DATA' block encountered in the input.
	 */
	public List<Sequence> parseDataBlock(List<Taxon> taxonList) throws ImportException, IOException
	{
		return readDataBlock(taxonList);
	}

	/**
	 * Parse the next 'TREES' block encountered in the input.
	 */
	public List<Tree> parseTreesBlock(List<Taxon> taxonList) throws ImportException, IOException
	{
		return readTreesBlock(taxonList);
	}

    /**
	 * Parse the next 'DISTANCES' block encountered in the input.
	 */
    public DistanceMatrix parseDistancesBlock(List<Taxon> taxonList) throws ImportException, IOException
	{
		return readDistancesBlock(taxonList);
	}

	// **************************************************************
	// SequenceImporter IMPLEMENTATION
	// **************************************************************

	/**
	 * Import all alignments in the input from the current position.
	 */
	public List<Alignment> importAlignments() throws IOException, ImportException
	{
		boolean done = false;

		List<Taxon> taxonList = null;
		List<Alignment> alignments = new ArrayList<Alignment>();

		while (!done) {
			try {

				NexusBlock block = findNextBlock();

				if (block == NexusBlock.TAXA) {

					taxonList = readTaxaBlock();

				} else if (block == NexusBlock.CHARACTERS) {

					if (taxonList == null) {
						throw new MissingBlockException("TAXA block is missing");
					}

					List<Sequence> sequences = readCharactersBlock(taxonList);
					alignments.add(new BasicAlignment(sequences));

				} else if (block == NexusBlock.DATA) {

					// A data block doesn't need a taxon block before it
					// but if one exists then it will use it.
					List<Sequence> sequences = readDataBlock(taxonList);
					alignments.add(new BasicAlignment(sequences));

				} else {
					// Ignore the block..
				}

			} catch (EOFException ex) {
				done = true;
			}
		}

		if (alignments.size() == 0) {
			throw new MissingBlockException("DATA or CHARACTERS block is missing");
		}

		return alignments;
	}

	/**
	 * Import all sequences in the input from the current position
	 */
	public List<Sequence> importSequences() throws IOException, ImportException {
		boolean done = false;

		List<Taxon> taxonList = null;
		List<Sequence> sequences = null;

		while (!done) {
			try {
				NexusBlock block = findNextBlock();
				if (block == NexusBlock.TAXA) {
					taxonList = readTaxaBlock();
				} else if (block == NexusBlock.CHARACTERS) {
					if (taxonList == null) {
						throw new MissingBlockException("TAXA block is missing");
					}
					sequences = readCharactersBlock(taxonList);
					done = true;
				} else if (block == NexusBlock.DATA) {
					// A data block doesn't need a taxon block before it
					// but if one exists then it will use it.
					sequences = readDataBlock(taxonList);
					done = true;
				} else {
					// Ignore the block..
				}
			} catch (EOFException ex) {
				done = true;
			}
		}
		if (sequences == null) {
			throw new MissingBlockException("DATA or CHARACTERS block is missing");
		}
		return sequences;
	}

	// **************************************************************
	// TreeImporter IMPLEMENTATION
	// **************************************************************

	private boolean isReadingTreesBlock = false;
	private List<Taxon> treeTaxonList = null;
	private Map<String, Taxon> translationMap = Collections.emptyMap();
	private Tree nextTree = null;
	private String[] lastToken = new String[1];

	/**
     * If not currently reading a TREES block then read ahead to the next TREES block, parsing TRANSLATE and TAXA blocks
     * in the process if necessary.
     * <p/>
     * Then determine if the current (or next) TREES block contains any more trees. If true then importNextTree will return a
     * non-null value.
     *
	 * @return true if another tree is available, false otherwise.
	 */
	public boolean hasTree() throws IOException, ImportException {
		if (!isReadingTreesBlock) {
			isReadingTreesBlock = startReadingTrees();
			translationMap = readTranslationMap(treeTaxonList, lastToken);
		}

		if (!isReadingTreesBlock) return false;

		if (nextTree == null) {
			nextTree = readNextTree(lastToken);
		}

		return (nextTree != null);
	}


	/**
     * If not currently reading a TREES block then read ahead to the next TREES block, parsing TRANSLATE and TAXA blocks
     * in the process if necessary.
     * <p/>
	 * Then parse the next available tree.
     *
	 * @return the next available tree or null if no more trees are available
	 */
	public Tree importNextTree() throws IOException, ImportException
	{
		// call hasTree to do the hard work...
		if (!hasTree()) {
			isReadingTreesBlock = false;
			return null;
		}

		Tree tree = nextTree;
		nextTree = null;

		return tree;
	}

    /**
     * Import all trees in the file from the current position. Will read ahead to the next TREES block if necessary
     *
     * @return list of trees
     * @throws IOException
     * @throws ImportException
     */
    public List<Tree> importTrees() throws IOException, ImportException {
        // We can't call startReadingTrees() here because if hasTree() was called before
        // then this importer will already have read into the trees block. However
        // is hasTree() was called then the following is still guaranteed to work as
        // per the TreeImporter.hasTree() javadoc.
        List<Tree> result = new ArrayList<Tree>();
        while (hasTree()) {
            result.add(importNextTree());
        }
        if (result.isEmpty()) {
            throw new MissingBlockException("TREES block is missing");
        } else {
            return Collections.unmodifiableList(result);
        }
	}

    /**
     * Read nexus blocks until the next TREES block (or the end of the input) is encountered.
     *
     * @return true if a TREES block was found, false otherwise
     * @throws IOException
     * @throws ImportException
     */
    public boolean startReadingTrees() throws IOException, ImportException
	{
		treeTaxonList = null;

		while (true) {
			try {
				NexusBlock block = findNextBlock();
				switch( block )  {
					case TAXA: treeTaxonList = readTaxaBlock(); break;
					case TREES: return true;
						// Ignore the block..
					default: break;
				}
			} catch (EOFException ex) {
				break;
			}
		}

		return false;
	}

	/**
	 * Import all distance matrices from all DISTANCES blocks in the input form the current position.
	 */
	public List<DistanceMatrix> importDistanceMatrices() throws IOException, ImportException {
		boolean done = false;
		List<Taxon> taxonList = null;
		List<DistanceMatrix> distanceMatrices = new ArrayList<DistanceMatrix>();
		while (!done) {
			try {
				NexusBlock block = findNextBlock();
				if (block == NexusBlock.TAXA) {
					taxonList = readTaxaBlock();
				} else if (block == NexusBlock.DISTANCES) {
					if (taxonList == null) {
						throw new MissingBlockException("TAXA block is missing");
					}
					DistanceMatrix distanceMatrix = readDistancesBlock(taxonList);
					distanceMatrices.add(distanceMatrix);
				} else {
					// Ignore the block..
				}
			} catch (EOFException ex) {
				done = true;
			}
		}
		return distanceMatrices;
	}

	// **************************************************************
	// PRIVATE Methods
	// **************************************************************

	/**
	 * Finds the end of the current block.
	 */
	private void findToken(String query, boolean ignoreCase) throws IOException
	{
		String token;
		boolean found = false;

		do {
			token = helper.readToken();

			if ( (ignoreCase && token.equalsIgnoreCase(query)) || token.equals(query) ) {
				found = true;
			}
		} while (!found);
	}

	/**
	 * Read ahead to the end of the current block.
	 */
	public void findEndBlock() throws IOException
	{
		try {
			String token;

			do {
				token = helper.readToken(";");
			} while ( !token.equalsIgnoreCase("END") && !token.equalsIgnoreCase("ENDBLOCK") );
		} catch (EOFException e) {
			// Doesn't matter if the End is missing
		}

		nextBlock = NexusBlock.UNKNOWN;
	}

	/**
	 * Reads the header information (DIMENSIONS, TITLE and FORMAT) for the current 'DATA', 'CHARACTERS' or 'TAXA' block.
     * Requires that the importer is already at the block's position (this will not read ahead to the specified block).
     *
     * @param tokenToLookFor the first token which is expected to occur after the header fields. once this token has is
     * reached, the method will return.
     * @param block the type of the block currently being parsed.
	 */
	private void readDataBlockHeader(String tokenToLookFor, NexusBlock block) throws ImportException, IOException
	{

		boolean foundDimensions = false, foundTitle = false, foundFormat = false;
		String token;

		do {
			token = helper.readToken();

			if ( token.equalsIgnoreCase("TITLE") ) {
				if (foundTitle) {
					throw new ImportException.DuplicateFieldException("TITLE");
				}

				foundTitle = true;
			} else if ( token.equalsIgnoreCase("DIMENSIONS") ) {

				if (foundDimensions) {
					throw new ImportException.DuplicateFieldException("DIMENSIONS");
				}

				boolean nchar = (block == NexusBlock.TAXA);
				boolean ntax = (block == NexusBlock.CHARACTERS);

				do {
					String token2 = helper.readToken( "=;" );

					if (helper.getLastDelimiter() != '=') {
						throw new ImportException.BadFormatException("Unknown subcommand, '" + token2 + "', or missing '=' in DIMENSIONS command");
					}

					if ( token2.equalsIgnoreCase("NTAX") ) {

						if (block == NexusBlock.CHARACTERS) {
							throw new ImportException.BadFormatException("NTAX subcommand in CHARACTERS block");
						}

						taxonCount = helper.readInteger( ";" );
						ntax = true;

					} else if ( token2.equalsIgnoreCase("NCHAR") ) {

						if (block == NexusBlock.TAXA) {
							throw new ImportException.BadFormatException("NCHAR subcommand in TAXA block");
						}

						siteCount = helper.readInteger( ";" );
						nchar = true;

					} else {
						throw new ImportException.BadFormatException("Unknown subcommand, '" + token2 + "', in DIMENSIONS command");
					}

				} while (helper.getLastDelimiter() != ';');

				if (!ntax) {
					throw new ImportException.BadFormatException("NTAX subcommand missing from DIMENSIONS command");
				}
				if (!nchar) {
					throw new ImportException.BadFormatException("NCHAR subcommand missing from DIMENSIONS command");
				}
				foundDimensions = true;

			} else if ( token.equalsIgnoreCase("FORMAT") ) {

				if (foundFormat) {
					throw new ImportException.DuplicateFieldException("FORMAT");
				}

				sequenceType = null;

				do {
					String token2 = helper.readToken("=;");

					if (token2.equalsIgnoreCase("GAP")) {

						if (helper.getLastDelimiter() != '=') {
							throw new ImportException.BadFormatException("Expecting '=' after GAP subcommand in FORMAT command");
						}

						gapCharacters = helper.readToken(";");

					} else if (token2.equalsIgnoreCase("MISSING")) {

						if (helper.getLastDelimiter() != '=') {
							throw new ImportException.BadFormatException("Expecting '=' after MISSING subcommand in FORMAT command");
						}

						missingCharacters = helper.readToken(";");

					} else if (token2.equalsIgnoreCase("MATCHCHAR")) {

						if (helper.getLastDelimiter() != '=') {
							throw new ImportException.BadFormatException("Expecting '=' after MATCHCHAR subcommand in FORMAT command");
						}

						matchCharacters = helper.readToken(";");

					} else if (token2.equalsIgnoreCase("DATATYPE")) {

						if (helper.getLastDelimiter() != '=') {
							throw new ImportException.BadFormatException("Expecting '=' after DATATYPE subcommand in FORMAT command");
						}

						String token3 = helper.readToken(";");
						if (token3.equalsIgnoreCase("NUCLEOTIDE") ||
								token3.equalsIgnoreCase("DNA") ||
								token3.equalsIgnoreCase("RNA")) {

							sequenceType = SequenceType.NUCLEOTIDE;

						} else if (token3.equalsIgnoreCase("PROTEIN")) {

							sequenceType = SequenceType.AMINO_ACID;

						} else if (token3.equalsIgnoreCase("CONTINUOUS")) {

							throw new ImportException.UnparsableDataException("Continuous data cannot be parsed at present");

						}
					} else if (token2.equalsIgnoreCase("INTERLEAVE")) {
						isInterleaved = true;
					}

				} while (helper.getLastDelimiter() != ';');

				foundFormat = true;
			}
		} while ( !token.equalsIgnoreCase(tokenToLookFor) );

		if ( !foundDimensions ) {
			throw new ImportException.MissingFieldException("DIMENSIONS");
		}
		if ( block != NexusBlock.TAXA && sequenceType == null ) {
			throw new ImportException.MissingFieldException("DATATYPE. Only Nucleotide or Protein sequences are supported.");
		}
	}

	/**
	 * Reads sequences in a 'DATA' or 'CHARACTERS' block.
	 */
	private List<Sequence> readSequenceData(List<Taxon> taxonList) throws ImportException, IOException
	{
        boolean sequencherStyle = false;
        String firstSequence = null;
		List<Sequence> sequences = new ArrayList<Sequence>();

		if (isInterleaved) {
			List<String> sequencesData = new ArrayList<String>(taxonCount);
			List<Taxon> taxons =  new ArrayList<Taxon>();
			List<Taxon> taxList =  (taxonList != null) ? taxonList : taxons;

			int[] charsRead = new int[taxonCount];
			for (int i = 0; i < taxonCount; i++) {
				sequencesData.add("");
				charsRead[i] = 0;
			}
			//throw new ImportException.UnparsableDataException("At present, interleaved data is not parsable");
			boolean firstLoop = true;

			int readCount = 0;
			while (readCount < siteCount * taxonCount) {

				for (int i = 0; i < taxonCount; i++) {

					String token = helper.readToken();

					int sequenceIndex;
					Taxon taxon = Taxon.getTaxon(token);
					if (firstLoop) {
						if (taxonList != null ) {
							sequenceIndex = taxonList.indexOf(taxon);
						} else {
							sequenceIndex = taxons.size();
							taxons.add(taxon);
						}
					} else {
						sequenceIndex = taxList.indexOf(taxon);
					}

					if( sequenceIndex < 0 ) {
						// taxon not found in taxon list...
						// ...perhaps it is a numerical taxon reference?
						throw new ImportException.UnknownTaxonException("Unexpected taxon:" + token
								+ " (expecting " + taxList.get(i).getName() + ")");
					}

					StringBuffer buffer = new StringBuffer();

					helper.readSequenceLine(buffer, sequenceType, ";", gapCharacters, missingCharacters,
							matchCharacters, firstSequence);

					String seqString = buffer.toString();

                    // We now check if this file is in Sequencher* style NEXUS, this style has the taxon and site counts
                    // before the sequence data.
                    try{
                        if(firstLoop && Integer.parseInt(taxon.toString()) == taxonCount &&
                                Integer.parseInt(seqString) == siteCount){
                            i--;
                            taxons.remove(taxon);
                            sequencherStyle = true;
                            continue;
                        }
                    } catch(NumberFormatException e) {
                        // Do nothing, this just means that this is the NEXUS format we usually expect rather than sequencher
                    }

                    readCount += seqString.length();
					charsRead[sequenceIndex] += seqString.length();

					sequencesData.set(sequenceIndex, sequencesData.get(sequenceIndex).concat(seqString));
					if (i == 0) {
						firstSequence = seqString;
					}

					if (helper.getLastDelimiter() == ';') {
						if (i < taxonCount - 1) {
							throw new ImportException.TooFewTaxaException();
						}
						for (int k = 0; k < taxonCount; k++) {
							if (charsRead[k] != siteCount) {
								throw new ImportException.ShortSequenceException(taxList.get(k).getName()
                                        + " has length " + charsRead[k] + ", expecting " + siteCount);
							}
						}
					}
                }

				firstLoop = false;
			}

            // Sequencher style apparently doesnt use a ';' after the sequence data.
            if (!sequencherStyle && helper.getLastDelimiter() != ';') {
				throw new ImportException.BadFormatException("Expecting ';' after sequences data");
			}

			for (int k = 0; k < taxonCount; k++) {
				Sequence sequence = new BasicSequence(sequenceType, taxList.get(k), sequencesData.get(k));
				sequences.add(sequence);
			}

		} else {

			for (int i = 0; i < taxonCount; i++) {
				String token = helper.readToken();

				Taxon taxon = Taxon.getTaxon(token);

				if (taxonList != null && !taxonList.contains(taxon)) {
					// taxon not found in taxon list...
					// ...perhaps it is a numerical taxon reference?
                    StringBuilder message = new StringBuilder("Expected: ").append(token).append("\nActual taxa:\n");
                    for (Taxon taxon1 : taxonList) {
                        message.append(taxon1).append("\n");
                    }
                    throw new ImportException.UnknownTaxonException(message.toString());
				}

				StringBuilder buffer = new StringBuilder() ;
				helper.readSequence(buffer, sequenceType, ";", siteCount, gapCharacters,
						missingCharacters, matchCharacters, firstSequence, true);
				String seqString = buffer.toString();

				if (seqString.length() != siteCount) {
					throw new ImportException.ShortSequenceException(taxon.getName()
                            + " has length " + seqString.length() + ", expecting " + siteCount);
				}

				if (i == 0) {
					firstSequence = seqString;
				}

                if (helper.getLastDelimiter() == ';' && i < taxonCount - 1) {
					throw new ImportException.TooFewTaxaException();
				}

				Sequence sequence = new BasicSequence(sequenceType, taxon, seqString);
				sequences.add(sequence);
			}

			if (helper.getLastDelimiter() != ';') {
				throw new ImportException.BadFormatException("Expecting ';' after sequences data");
			}

		}

		return sequences;
	}


	/**
	 * Reads a 'TAXA' block.
	 */
	private List<Taxon> readTaxaBlock() throws ImportException, IOException {
		taxonCount = 0;

		readDataBlockHeader("TAXLABELS", NexusBlock.TAXA);

		if (taxonCount == 0) {
			throw new ImportException.MissingFieldException("NTAXA");
		}

		List<Taxon> taxa = new ArrayList<Taxon>();

		do {
            String name = helper.readToken(";");
            if (name.equals("")) {
                throw new ImportException.UnknownTaxonException("Expected nonempty taxon name, got empty string");
            }
            Taxon taxon = Taxon.getTaxon(name);
			taxa.add(taxon);

            parseAndClearMetaComments(taxon, helper);
		} while (helper.getLastDelimiter() != ';');

		if (taxa.size() != taxonCount) {
			throw new ImportException.BadFormatException("Number of taxa doesn't match NTAXA field");
		}

		findEndBlock();

		return taxa;
	}

    /**
     * Parse any meta comment pairs that have been read and clear them from the ImportHelper. Comment pairs will be attached
     * to the given Attributable.
     *
     * @param item Attributable to attach comments too
     * @param importHelper ImportHelper which may have read meta comments.
     * @throws ImportException.BadFormatException
     */
    static void parseAndClearMetaComments(Attributable item, ImportHelper importHelper) throws ImportException.BadFormatException {
        for (String meta : importHelper.getMetaComments()) {
            // A meta-comment which should be in the form:
            // \[&label[=value][,label[=value]>[,/..]]\]
            parseMetaCommentPairs(meta, item);

        }
        importHelper.clearLastMetaComment();
    }

    /**
	 * Reads a 'CHARACTERS' block.
	 */
	private List<Sequence> readCharactersBlock(List<Taxon> taxonList) throws ImportException, IOException
	{

		siteCount = 0;
		sequenceType = null;

		readDataBlockHeader("MATRIX", NexusBlock.CHARACTERS);

		List<Sequence> sequences = readSequenceData(taxonList);

		findEndBlock();

		return sequences;
	}

	/**
	 * Reads a 'DATA' block.
	 */
	private List<Sequence> readDataBlock(List<Taxon> taxonList) throws ImportException, IOException
	{

		taxonCount = 0;
		siteCount = 0;
		sequenceType = null;

		readDataBlockHeader("MATRIX", NexusBlock.DATA);

		List<Sequence> sequences = readSequenceData(taxonList);

		findEndBlock();

		return sequences;
	}

	/**
	 * Reads a 'DISTANCES' block.
	 */
	private DistanceMatrix readDistancesBlock(List<Taxon> taxonList) throws ImportException, IOException
	{
		if (taxonList == null) {
			throw new ImportException.BadFormatException("Missing Taxa for reading distances");
		}

		Triangle triangle = Triangle.LOWER;
		boolean diagonal = true;
		boolean labels = false;

		boolean ttl = false, fmt = false;

		String token = helper.readToken();
		while ( !token.equalsIgnoreCase("MATRIX") ) {

			if ( token.equalsIgnoreCase("TITLE") ) {
				if (ttl) {
					throw new ImportException.DuplicateFieldException("TITLE");
				}

				ttl = true;
			} else if ( token.equalsIgnoreCase("FORMAT") ) {

				if (fmt) {
					throw new ImportException.DuplicateFieldException("FORMAT");
				}

				sequenceType = null;

				do {
					String token2 = helper.readToken("=;");

					if (token2.equalsIgnoreCase("TRIANGLE")) {

						if (helper.getLastDelimiter() != '=') {
							throw new ImportException.BadFormatException("Expecting '=' after TRIANGLE subcommand in FORMAT command");
						}

						String token3 = helper.readToken(";");
						if (token3.equalsIgnoreCase("LOWER")) {
							triangle = Triangle.LOWER;
						} else if (token3.equalsIgnoreCase("UPPER")) {
							triangle = Triangle.UPPER;
						} else if (token3.equalsIgnoreCase("BOTH")) {
							triangle = Triangle.BOTH;
						}
					} else
						// based on the example in the PAUP manual
						if (token2.equalsIgnoreCase("NODIAGONAL")) {
							// in a valid file, triangle != both
							diagonal = false;
						} else if (token2.equalsIgnoreCase("LABELS")) {
							labels = true;
						}

				} while (helper.getLastDelimiter() != ';');

				fmt = true;
			}
			token = helper.readToken();
		}

		double[][] distances = new double[taxonList.size()][taxonList.size()];

		for (int i = 0; i < taxonList.size(); i++) {
			token = helper.readToken();

			Taxon taxon = Taxon.getTaxon(token);

			int index = taxonList.indexOf(taxon);

			if (index < 0) {
				// taxon not found in taxon list...
				// ...perhaps it is a numerical taxon reference?
                StringBuilder message = new StringBuilder("Expected: ").append(token).append("\nActual taxa:\n");
                for (Taxon taxon1 : taxonList) {
                    message.append(taxon1).append("\n");
                }
                throw new ImportException.UnknownTaxonException(message.toString());
			}

			if (index != i) {
				throw new ImportException.BadFormatException("The taxon labels are in a different order to those in the TAXA block");
			}

			if (triangle == Triangle.LOWER) {
				for (int j = 0; j < i + 1; j++) {
					if (i != j) {
						distances[i][j] = helper.readDouble();
						distances[j][i] = distances[i][j];
					} else {
						if (diagonal) {
							distances[i][j] = helper.readDouble();
						}
					}
				}
			} else if (triangle == Triangle.UPPER) {
				for (int j = i; j < taxonList.size(); j++) {
					if (i != j) {
						distances[i][j] = helper.readDouble();
						distances[j][i] = distances[i][j];
					} else {
						if (diagonal) {
							distances[i][j] = helper.readDouble();
						}
					}
				}
			} else {
				for (int j = 0; j < taxonList.size(); j++) {
					if (i != j || diagonal) {
						distances[i][j] = helper.readDouble();
					} else {
						distances[i][j] = 0.0;
					}
				}
			}


			if (helper.getLastDelimiter() == ';' && i < taxonList.size() - 1) {
				throw new ImportException.TooFewTaxaException();
			}
		}

		if (helper.nextCharacter() != ';') {
			throw new ImportException.BadFormatException("Expecting ';' after sequences data");
		}

		findEndBlock();

		return new BasicDistanceMatrix(taxonList, distances);
	}


	/**
	 * Reads a 'TREES' block.
	 */
	private List<Tree> readTreesBlock(List<Taxon> taxonList) throws ImportException, IOException {
		List<Tree> trees = new ArrayList<Tree>();
		String[] lastToken = new String[1];
		translationMap = readTranslationMap(taxonList, lastToken);
        while( true ) {
			RootedTree tree = readNextTree(lastToken);
			if (tree == null) {
				break;
			}
			trees.add(tree);
		}

		if (trees.size() == 0) {
			throw new ImportException.BadFormatException("No trees defined in TREES block");
		}

		nextBlock = NexusBlock.UNKNOWN;
		return trees;
	}

	private Map<String, Taxon> readTranslationMap(List<Taxon> taxonList, String[] lastToken) throws ImportException, IOException {
		Map<String, Taxon> translationList = new HashMap<String, Taxon>();

		String token = helper.readToken(";");

		if ( token.equalsIgnoreCase("TRANSLATE") ) {

			do {
				String token2 = helper.readToken(",;");

				if (helper.getLastDelimiter() == ',' || helper.getLastDelimiter() == ';') {
                    if( token2.length() == 0 && (char)helper.getLastDelimiter() == ';') {
                        //assume an extra comma at end of list
                        break;
                    }
                    throw new ImportException.BadFormatException("Missing taxon label in TRANSLATE command of TREES block");
				}

				String token3 = helper.readToken(",;");

				if (helper.getLastDelimiter() != ',' && helper.getLastDelimiter() != ';') {
					throw new ImportException.BadFormatException("Expecting ',' or ';' after taxon label in TRANSLATE command of TREES block");
				}

				Taxon taxon = Taxon.getTaxon(token3);

				if (taxonList != null && !taxonList.contains(taxon)) {
					// taxon not found in taxon list...
					// ...perhaps it is a numerical taxon reference?
                    StringBuilder message = new StringBuilder("Expected: ").append(token).append("\nActual taxa:\n");
                    for (Taxon taxon1 : taxonList) {
                        message.append(taxon1).append("\n");
                    }
					throw new ImportException.UnknownTaxonException(message.toString());
				}
				translationList.put(token2, taxon);
			} while (helper.getLastDelimiter() != ';');
			token = helper.readToken(";");
		} else if (taxonList != null) {
			for (Taxon taxon : taxonList) {
				translationList.put(taxon.getName(), taxon);
			}
		}
		lastToken[0] = token;
		return translationList;
	}

	private RootedTree readNextTree(String[] lastToken) throws ImportException, IOException
	{
		try {
			SimpleRootedTree tree = null;
			String token = lastToken[0];

            boolean isUnrooted = token.equalsIgnoreCase("UTREE");
            if ( isUnrooted || token.equalsIgnoreCase("TREE")) {
				if (helper.nextCharacter() == '*') {
					// Star is used to specify a default tree - ignore it
					helper.readCharacter();
				}

                {
                    // According to the Nexus specification at http://www.cs.nmsu.edu/~epontell/nexus/nexus_grammar
                    // and all Nexus files I (TT) have seen, the [&U] unrooted meta comment must actually occur not
                    // here but only after the ' = ' behind the name. However until 2008-05-05 JEBL produced Nexus
                    // files with the [&U] in this wrong location so we need to continue supporting such broken
                    // files; this has caused bug 5150.
                    for (String meta : helper.getMetaComments()) {
                        // Look for the unrooted meta comment [&U]
                        if (meta.equalsIgnoreCase("U")) {
                            isUnrooted = true;
                        }
                    }
                    helper.clearLastMetaComment();
                }

				final String treeName = helper.readToken( "=;" );

				if (helper.getLastDelimiter() != '=') {
					throw new ImportException.BadFormatException("Missing label for tree '" + treeName + "' or missing '=' in TREE command of TREES block");
				}

				try {
					if (helper.nextCharacter() != '(') {
						throw new ImportException.BadFormatException("Missing tree definition in TREE command of TREES block");
					}

					// Save tree comment and attach it later
					final List<String> comments = helper.getMetaComments();
					helper.clearLastMetaComment();

					tree = new SimpleRootedTree();
					readInternalNode(tree);

					// save name as attribute
                    if( ! NexusExporter.isGeneratedTreeName(treeName) )  {
                        tree.setAttribute(NexusExporter.treeNameAttributeKey, treeName);
                    }

                    int last = helper.getLastDelimiter();
					if( last == ':' ) {
						// root length - discard for now
						/*double rootLength = */ helper.readDouble(";");
						last = helper.getLastDelimiter();
					}

					if (last != ';') {
						throw new ImportException.BadFormatException("Expecting ';' after tree, '" + treeName + "', TREE command of TREES block");
					}

                    for (String comment : comments) {
                        if (Character.toUpperCase(comment.charAt(0)) == 'U') { // [&U] unrooted meta comment, see tree_rest, root in http://www.cs.nmsu.edu/~epontell/nexus/nexus_grammar
                            isUnrooted = true;
                        } else if (comment.matches("^W\\s+[\\+\\-]?[\\d\\.]+")) { // if '[W number]' (MrBayes), set weight attribute
                            tree.setAttribute("weight", Float.valueOf(comment.substring(2)));
						} else {
							try {
								parseMetaCommentPairs(comment, tree);
							} catch(ImportException.BadFormatException e) {
								// set generic comment attribute
								tree.setAttribute("comment", comment);
							}
						}
					}

					tree.setConceptuallyUnrooted(isUnrooted);

				} catch (EOFException e) {
					// If we reach EOF we may as well return what we have?
					return tree;
				}

				token = helper.readToken(";");
			} else if ( token.equalsIgnoreCase("ENDBLOCK") || token.equalsIgnoreCase("END") ) {
				return null;
			} else {
				throw new ImportException.BadFormatException("Unknown command '" + token + "' in TREES block");
			}

			//added this to escape readNextTree loop correctly -- AJD
			lastToken[0] = token;

			return compactTrees ? new CompactRootedTree(tree) : tree;

		} catch (EOFException e) {
			return null;
		}
	}

	/**
	 * Reads a branch in. This could be a node or a tip (calls readNode or readTip
	 * accordingly). It then reads the branch length and SimpleNode that will
	 * point at the new node or tip.
	 */
	private Node readBranch(SimpleRootedTree tree) throws IOException, ImportException
	{
		Node branch;

		helper.clearLastMetaComment();
		if (helper.nextCharacter() == '(') {
			// is an internal node
			branch = readInternalNode(tree);

		} else {
			// is an external node
			branch = readExternalNode(tree);
		}

		if (helper.getLastDelimiter() == ':') {
			final double length = helper.readDouble(",():;");
			tree.setLength(branch, length);
		}

		// If there is a metacomment after the branch length indicator (:), then it is a branch attribute
		// however, in the present implementation, this simply gets added to the node attributes.
        parseAndClearMetaComments(tree.getParentEdge(branch), helper);

		return branch;
	}

	/**
	 * Reads a node in. This could be a polytomy. Calls readBranch on each branch
	 * in the node.
     * @param tree
     * @return
     */
	private Node readInternalNode(SimpleRootedTree tree) throws IOException, ImportException
	{
		List<Node> children = new ArrayList<Node>();

		// read the opening '('
		helper.readCharacter();

		// read the first child
		children.add( readBranch(tree) );

		if (helper.getLastDelimiter() != ',') {
			//throw new ImportException.BadFormatException("Missing ',' in tree");
		}
        // MK: previously, an internal node must have at least 2 children.
        // MK: We we now allow trees with a single child so that we can create proper taxonomy
        // MK: trees with only a single child at a taxonomy level.

        // read subsequent children

        while(helper.getLastDelimiter()==',') {
			children.add( readBranch(tree) );
		}

        //System.out.println("kids="+children.size());
        // should have had a closing ')'
		if (helper.getLastDelimiter() != ')') {
			throw new ImportException.BadFormatException("Missing closing ')' in tree");
		}

		Node node = tree.createInternalNode(children);

		// find the next delimiter
		String token = helper.readToken(":(),;").trim();

		// if there is a token before the branch length, treat it as a node label
		// and store it as an attribute of the node...
		if (token.length() > 0) {
			node.setAttribute("label", parseValue(token));
		}

		// If there is a metacomment before the branch length indicator (:), then it is a node attribute
        parseAndClearMetaComments(node, helper);

		return node;
	}

	/**
	 * Reads an external node in.
	 */
	private Node readExternalNode(SimpleRootedTree tree) throws ImportException, IOException
	{
		String label = helper.readToken(":(),;");

        Taxon taxon;
        try {
            taxon = Taxon.getTaxon(label);
        } catch (IllegalArgumentException e) {
            throw new ImportException.UnknownTaxonException(e.getMessage());
        }

        if (translationMap.size() > 0) {
			taxon = translationMap.get(label);

			if (taxon == null) {
				// taxon not found in taxon list...
				throw new ImportException.UnknownTaxonException("Taxon in tree, '" + label + "' is unknown");
			}
		}

        try {
            final Node node = tree.createExternalNode(taxon);

	        // Attempt to parse external node attributes
	        // If there is a metacomment before the branch length indicator (:), then it is a node attribute
            parseAndClearMetaComments(node, helper);

            return node;
        } catch (IllegalArgumentException e) {
           throw new ImportException.DuplicateTaxaException(e.getMessage());
        }
    }

	static void parseMetaCommentPairs(String meta, Attributable item) throws ImportException.BadFormatException {
		// This regex should match key=value pairs, separated by commas
		// This can match the following types of meta comment pairs:
		// value=number, value="string", value={item1, item2, item3}
        // (label must be quoted if it contains spaces (i.e. "my label"=label)

        Pattern pattern = Pattern.compile("(\"[^\"]*\"+|[^,=\\s]+)\\s*(=\\s*(\\{[^=}]*\\}|\"[^\"]*\"+|[^,]+))?");
		Matcher matcher = pattern.matcher(meta);

		while (matcher.find()) {
			String label = matcher.group(1);
            if( label.charAt(0) == '\"' ) {
                label = label.substring(1, label.length() - 1);
            }
            if (label == null || label.trim().length() == 0) {
				throw new ImportException.BadFormatException("Badly formatted attribute: '"+ matcher.group()+"'");
			}
			final String value = matcher.group(2);
			if (value != null && value.trim().length() > 0) {
				// there is a specified value so try to parse it
				item.setAttribute(label, parseValue(value.substring(1)));
			} else {
				item.setAttribute(label, Boolean.TRUE);
			}
		}
	}

	/**
	 * This method takes a string and tries to decode it returning the object
	 * that best fits the data. It will recognize command delimited lists enclosed
	 * in {..} and call parseValue() on each element. It will also recognize Boolean,
	 * Integer and Double. If the value starts with a # then it will attempt to decode
	 * the following integer as an RGB colour - see Color.decode(). If nothing else fits
	 * then the value will be returned as a string but trimmed of leading and trailing
	 * white space.
	 * @param value the string
	 * @return the object
	 */
	static Object parseValue(String value) {

		value = value.trim();

		if (value.startsWith("{")) {
			// the value is a list so recursively parse the elements
			// and return an array
			String[] elements = value.substring(1, value.length() - 1).split(",");
			Object[] values = new Object[elements.length];
			for (int i = 0; i < elements.length; i++) {
				values[i] = parseValue(elements[i]);
			}
			return values;
		}

		if (value.startsWith("#")) {
			// I am not sure whether this is a good idea but
			// I am going to assume that a # denotes an RGB colour
			try {
				return Color.decode(value.substring(1));
			} catch (NumberFormatException nfe1) {
				// not a colour
			}
		}

        // A string qouted by the nexus exporter and such
        if( value.startsWith("\"") && value.endsWith("\"") ) {
            return value.subSequence(1, value.length() - 1);
        }

        if (value.equalsIgnoreCase("TRUE") || value.equalsIgnoreCase("FALSE")) {
			return Boolean.valueOf(value);
		}

		// Attempt to format the value as an integer
		try {
			return Integer.parseInt(value);
		} catch (NumberFormatException nfe1) {
			// not an integer
		}

		// Attempt to format the value as a double
		try {
			return Double.parseDouble(value);
		} catch (NumberFormatException nfe2) {
			// not a double
		}

		// return the trimmed string
		return value;
	}

	// private stuff
	private NexusBlock nextBlock = null;
	private String nextBlockName = null;

	private int taxonCount = 0, siteCount = 0;
	private SequenceType sequenceType = null;
	private String gapCharacters = "-";
	private String matchCharacters = ".";
	private String missingCharacters = "?";
	private boolean isInterleaved = false;

	protected final ImportHelper helper;
}