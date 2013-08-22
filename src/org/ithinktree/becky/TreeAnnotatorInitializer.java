/**
 * TreeAnnotatorInitializer.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky;

import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

import dr.app.tools.NexusExporter;
import dr.app.util.Arguments;
import dr.app.util.Arguments.ArgumentException;
import dr.app.util.Arguments.Option;
import dr.app.util.Arguments.StringOption;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.TaxonList;

/**
 * @author Arman D. Bilge
 *
 */
public class TreeAnnotatorInitializer {

	/**
	 * @throws IOException 
	 * @throws ImportException 
	 * 
	 */
	public TreeAnnotatorInitializer(final String hostTreesFile, final String symbiontTreesFile, final String hostTreesOutputFile, final String symbiontTreesOutputFile) throws IOException, ImportException {
		
		final FileReader hostFileReader = new FileReader(hostTreesFile);
		final FileReader symbiontFileReader = new FileReader(symbiontTreesFile);
		
		final TreeImporter hostTreeImporter = new NexusImporter(hostFileReader);
		final TreeImporter symbiontTreeImporter = new NexusImporter(symbiontFileReader);
		
		final PrintStream hostTreesStream = new PrintStream(new FileOutputStream(hostTreesOutputFile));
		final PrintStream symbiontTreesStream = new PrintStream(new FileOutputStream(symbiontTreesOutputFile));
		
		final NexusExporter hostNexusExporter = new NexusExporter(hostTreesStream);
		final NexusExporter symbiontNexusExporter = new NexusExporter(symbiontTreesStream);
		
		Map<String, Integer> hostNexusHeader = null;
		Map<String, Integer> symbiontNexusHeader = null;
		
		int totalTrees = 10000;
        final int stepSize = totalTrees / 60;
        
        PrintStream progressStream = System.err;
        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");
		totalTrees = 0;
		while (hostTreeImporter.hasTree() && symbiontTreeImporter.hasTree()) {
			MutableTree hostTree = (MutableTree) hostTreeImporter.importNextTree();
			MutableTree symbiontTree = (MutableTree) symbiontTreeImporter.importNextTree();
			if (hostTaxonList == null) hostTaxonList = hostTree;
			ids = new int[hostTree.getNodeCount()];
			addNewClades(hostTree, hostTree.getRoot());
			
			for (int i = 0; i < symbiontTree.getNodeCount(); ++i) {
				NodeRef n = symbiontTree.getNode(i);
				symbiontTree.setNodeAttribute(n, HOST_NODE_REF, ids[(Integer) symbiontTree.getNodeAttribute(n, HOST_NODE_REF)]);
			}
			if (hostNexusHeader == null) {
				hostNexusHeader = hostNexusExporter.writeNexusHeader(hostTree);
				symbiontNexusHeader = symbiontNexusExporter.writeNexusHeader(symbiontTree);
			}
			hostNexusExporter.writeNexusTree(hostTree, hostTree.getId(), true, hostNexusHeader);
			symbiontNexusExporter.writeNexusTree(symbiontTree, symbiontTree.getId(), true, symbiontNexusHeader);
			if (totalTrees++ % stepSize == 0) {
                progressStream.print("*");
                progressStream.flush();
            }
		}
		progressStream.println();
		progressStream.println("Read " + totalTrees + " trees.");
		hostTreesStream.println("End;");
		symbiontTreesStream.println("End;");
		progressStream.println();
		progressStream.println();
		progressStream.println("Done.");
		
		hostFileReader.close();
		symbiontFileReader.close();
		hostTreesStream.close();
		symbiontTreesStream.close();
		
	}
	
	private final Map<BitSet,Integer> clades = new HashMap<BitSet,Integer>();
	private int[] ids;
	private TaxonList hostTaxonList;
	private final String NODE_REF = "nodeRef";
	private final String HOST_NODE_REF = "host.nodeRef";
	
	private BitSet addNewClades(final MutableTree tree, final NodeRef node) {

		final BitSet bitSet = new BitSet();
		if (tree.isExternal(node)) {
			bitSet.set(hostTaxonList.getTaxonIndex(tree.getNodeTaxon(node).getId()));
		} else {
			for (int i = 0; i < tree.getChildCount(node); ++i) {
				bitSet.or(addNewClades(tree, tree.getChild(node, i)));
			}
		}
		
		if (!clades.containsKey(bitSet))
			clades.put(bitSet, clades.size());

		final int id = clades.get(bitSet);
		ids[(Integer) tree.getNodeAttribute(node, NODE_REF)] = id;
		tree.setNodeAttribute(node, NODE_REF, id);

		return bitSet;
	}

	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		final Arguments arguments = new Arguments(new Option[]{
				new StringOption("h", "filename", "host input trees file name"),
				new StringOption("i", "filename", "host output trees file name"),
				new StringOption("s", "filename", "symbiont input trees file name"),
				new StringOption("t", "filename", "symbiont output trees file name")
		}, false);
		
		try {
			arguments.parseArguments(args);
		} catch (ArgumentException e) {
			e.printStackTrace(System.err);
			arguments.printUsage("annotationinitializer", "");
			System.exit(1);
		}

		try {
			new TreeAnnotatorInitializer(arguments.getStringOption("h"), arguments.getStringOption("s"), arguments.getStringOption("i"), arguments.getStringOption("t"));
		} catch (Exception e) {
			e.printStackTrace(System.err);
			arguments.printUsage("annotationinitializer", "");
			System.exit(1);
		}
		
	}

}
