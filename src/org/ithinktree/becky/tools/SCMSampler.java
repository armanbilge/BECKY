package org.ithinktree.becky.tools;

import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Locale;

import org.ithinktree.becky.CophylogenyModel;
import org.ithinktree.becky.SimpleCophylogenyModel;
import org.ithinktree.becky.SimpleCophylogenyModel.EventType;

import dr.app.util.Arguments;
import dr.app.util.Arguments.ArgumentException;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Units;
import dr.inference.model.Parameter;

public class SCMSampler {
	
	private static boolean acceptTree(final Tree tree, final int extantSpeciesCount) {
		
		if (tree.getExternalNodeCount() != extantSpeciesCount) return false;
		
		for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
			final SimpleNode n = (SimpleNode) tree.getNode(i);
			if (n.getAttribute(CoevolutionSimulator.COEVOLUTIONARY_EVENT) == EventType.LOSS &&
					!tree.isRoot(n) && ((SimpleNode) tree.getParent(n)).getAttribute(CoevolutionSimulator.COEVOLUTIONARY_EVENT) != EventType.NO_EVENT) {
				final SimpleNode n_p = (SimpleNode) tree.getParent(n);
				if (n_p.getAttribute(CoevolutionSimulator.COEVOLUTIONARY_EVENT) != EventType.HOST_SWITCH) {
					return false;
				} else if (!tree.isRoot(n_p) && ((SimpleNode) tree.getParent(n_p)).getAttribute(CoevolutionSimulator.COEVOLUTIONARY_EVENT) != EventType.HOST_SWITCH) {
					return false;
				} else if (!((SimpleNode) CophylogenyModel.Utils.getSisters(tree, n).get(0)).getAttribute("host.nodeRef").equals(((SimpleNode) tree.getParent(n_p)).getAttribute("host.nodeRef"))) {
					return false;
				}
			}
		}
		return true;
	}
	
	public static void main(String[] args) throws IOException, ImportException {
		
		Locale.setDefault(Locale.US);
		Arguments arguments = new Arguments(new Arguments.Option[]{
				new Arguments.StringOption("h", "file name", "host tree file"),
				new Arguments.IntegerOption("t", "number of symbiont taxa"),
				new Arguments.RealArrayOption("r", 3, "coevolutionary rates"),
				new Arguments.RealOption("c", "relaxed clock stdev"),
				new Arguments.LongOption("seed", "random number generator seed"),
				new Arguments.IntegerOption("s", "sample size"),
				new Arguments.Option("n", "use newick for output")
		});
		
		try {
			arguments.parseArguments(args);
		} catch (ArgumentException e) {
			e.printStackTrace(System.err);
			arguments.printUsage("scmsampler", "");
			System.exit(1);
		}

		final FileReader hostFileReader = new FileReader(arguments.getStringOption("h"));
		final TreeImporter hostTreeImporter = new NexusImporter(hostFileReader);
		final Tree hostTree = hostTreeImporter.importNextTree();
		hostFileReader.close();
		
		final int taxonCount = arguments.getIntegerOption("t");
//		int uniqueTaxonCount = taxonCount;
		
		final CoevolutionSimulator sim = new CoevolutionSimulator();
		
		final SimpleCophylogenyModel model = new SimpleCophylogenyModel(new Parameter.Default(arguments.getRealArrayOption("r")[0]), new Parameter.Default(arguments.getRealArrayOption("r")[1]), new Parameter.Default(arguments.getRealArrayOption("r")[2]), Units.Type.YEARS);
		
//		final NexusExporter exporter = new NexusExporter(System.out);
//		Map<String,Integer> header = null;
		
		int totalTrees = arguments.getIntegerOption("s");
        int stepSize = totalTrees / 60;
        if (stepSize == 0) ++stepSize;
		PrintStream progressStream = System.err;
        progressStream.println("Simulating trees...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");
		
        PrintStream paramLog = new PrintStream(new FileOutputStream("param.log"));
        paramLog.println("STATE\tProbability");
        
        System.out.println("#NEXUS\nbegin trees;\n");
		for (int i = 0; i < arguments.getIntegerOption("s"); ++i) {
			
//			if (header == null) {
//				header = exporter.writeNexusHeader(hostTree);
//				System.out.println("\t\t;");
//			}
			
			Tree tree;
			do {
				tree = sim.simulateCoevolution(hostTree, 1.0, model, false, true);
			} while (taxonCount != -1 && !acceptTree(tree, taxonCount));
			
//			for (Taxon t : tree.asList()) {
//				if (!header.containsKey(t.toString())) header.put(t.toString(), ++uniqueTaxonCount);
//			}
//			exporter.writeNexusTree(tree, "TREE" + i+1, true, header);
			System.out.println("tree tree_" + (i+1) + " = " + Tree.Utils.newick(tree, new TreeTraitProvider[]{sim.provider}));
			paramLog.println(i+1 + "\t" + sim.getSimulationLogLikelihood());
			
			if (i % stepSize == 0) {
	            progressStream.print("*");
	            progressStream.flush();
	        }
		}
		
		System.out.println("\nend;");
	    paramLog.close();
        progressStream.println();
	}

}
