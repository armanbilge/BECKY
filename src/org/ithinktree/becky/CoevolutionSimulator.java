/**
 * CoevolutionSimulator.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import org.apache.commons.lang3.ArrayUtils;
import org.ithinktree.becky.CophylogenyModel.Utils;
import org.ithinktree.becky.CophylogenyModel.Utils.Relationship;

import dr.app.seqgen.SeqGen;
import dr.app.tools.NexusExporter;
import dr.app.util.Arguments;
import dr.app.util.Arguments.ArgumentException;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.MutableTaxonList;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * @author Arman D. Bilge
 */
public class CoevolutionSimulator {

	public CoevolutionSimulator() {} // Empty constructor, like the CoalescentSimulator
	
	/**
	 * Simulates a random coevolutionary process to generate a valid mapping for the starting tree.
	 * <p/>
	 * Truth be told, this hardly simulates coevolution, but instead does what is needed to generate a random, valid mapping.
	 * 
	 * @param hostTree
	 * @param symbiontTree
	 * @param hostAttributeName
	 */
	public void simulateCoevolution(Tree hostTree, MutableTree symbiontTree, CophylogenyLikelihood cophylogenyLikelihood, String hostAttributeName, boolean samplingNoHost) {
		
		for (int i = 0; i < symbiontTree.getExternalNodeCount(); i++) {
			NodeRef node = symbiontTree.getExternalNode(i);
			Taxon hostTaxon = (Taxon) symbiontTree.getNodeTaxon(node).getAttribute(hostAttributeName);
			NodeRef hostNode = null;
			if (hostTaxon != null) {
				hostNode = hostTree.getExternalNode(hostTree.getTaxonIndex(hostTaxon.getId()));
			}
			cophylogenyLikelihood.setStatesForNode(node, hostNode);
		}
		
		int[] postOrderList = new int[symbiontTree.getNodeCount()];
		Tree.Utils.postOrderTraversalList(symbiontTree, postOrderList);
		int[] hostNodeIds = MathUtils.shuffled(hostTree.getNodeCount() + (samplingNoHost ? 1 : 0));
		for (int i = 0; i < postOrderList.length; ++i) {
			NodeRef node = symbiontTree.getNode(postOrderList[i]);
			if (!symbiontTree.isExternal(node)) {
				
				NodeRef child1HostNode = cophylogenyLikelihood.getStatesForNode(symbiontTree.getChild(node, 0));
				NodeRef child2HostNode = cophylogenyLikelihood.getStatesForNode(symbiontTree.getChild(node, 1));

				NodeRef hostNode = null;
				EnumSet<Relationship> relationships;
				MathUtils.shuffle(hostNodeIds);
				int j = 0;
				do {
					if (j >= hostNodeIds.length) {
						i = -1; // Force the simulation to start over
						break;
					}
					relationships = EnumSet.noneOf(Relationship.class);
					int r = hostNodeIds[j++] - (samplingNoHost ? 1 : 0);
					hostNode = r == CophylogenyLikelihood.NO_HOST ? null : hostTree.getNode(r);
					relationships.add(Utils.determineRelationship(hostTree, hostNode, child1HostNode).relationship);
					relationships.add(Utils.determineRelationship(hostTree, hostNode, child2HostNode).relationship);
				} while (relationships.contains(Relationship.ANCESTOR) ||
							(relationships.contains(Relationship.SELF) && relationships.contains(Relationship.DESCENDANT)));
				cophylogenyLikelihood.setStatesForNode(node, hostNode);
				double maxHeight = hostTree.isRoot(hostNode) ? Double.MAX_VALUE : hostTree.getNodeHeight(hostTree.getParent(hostNode));
				double minHeight = Math.max(hostTree.getNodeHeight(hostNode), Math.max(symbiontTree.getNodeHeight(symbiontTree.getChild(node, 0)), symbiontTree.getNodeHeight(symbiontTree.getChild(node, 1))));
				if (minHeight > maxHeight) i = -1; // Don't know if this can happen, but restarts the simulation
				symbiontTree.setNodeHeight(node, MathUtils.nextDouble() * (maxHeight - minHeight) + minHeight);
			}
		}
	}
	
	public Tree simulateCoevolution(final Tree hostTree, final double rate, final SimpleCophylogenyModel model, final boolean isRelaxed) {
		if (!isRelaxed)
			return simulateCoevolution(hostTree, rate, model, false, 0.0);
		throw new IllegalArgumentException();
	}
	
	public Tree simulateCoevolution(final Tree hostTree, final double rate, final SimpleCophylogenyModel model, final boolean isRelaxed, final double stdev) {
		
		SimpleNode root;
		do {
			symbiontCounts = new int[hostTree.getTaxonCount()];
			associations.clear();
			root = simulateCoevolution(hostTree,
					hostTree.getRoot(),
					hostTree.getNodeHeight(hostTree.getRoot()) + MathUtils.nextExponential(hostTree.getNodeCount() / hostTree.getNodeHeight(hostTree.getRoot())),
					rate,
					model.getDuplicationRate(),
					model.getHostShiftRate(),
					model.getLossRate(),
					isRelaxed,
					stdev);
		} while (root == null);
		return new SimpleTree(root);
	}
	
	private int[] symbiontCounts;
	public final Map<String,String> associations = new HashMap<String,String>();
	private SimpleNode simulateCoevolution(final Tree hostTree, final NodeRef hostNode, final double height, final double rate, final double duplicationRate, final double hostShiftRate, final double lossRate, final boolean isRelaxed, final double stdev) {
				
		final SimpleNode node = new SimpleNode();
		
		final double relaxedRate;
		if (isRelaxed) {
			relaxedRate = Math.exp(rate + stdev * MathUtils.nextGaussian());
		} else {
			relaxedRate = rate;
		}
		
		node.setAttribute("cophylogeny.clock.rate", relaxedRate);
		node.setAttribute("host.nodeRef", hostNode.getNumber());
		
		final SimpleNode child1;
		final SimpleNode child2;
		
		final EventIndexAndTime nextEvent = simulateSimultaneousPoissonProcesses(relaxedRate * duplicationRate, relaxedRate * hostShiftRate, relaxedRate * lossRate);
		final double eventHeight = height - nextEvent.time;
		
		final double hostNodeHeight = hostTree.getNodeHeight(hostNode);
		if (hostNodeHeight > eventHeight) {
			node.setHeight(hostTree.getNodeHeight(hostNode));
			if (hostTree.isExternal(hostNode)) {
				// Cannot coevolve anymore
				int i = Integer.parseInt(hostTree.getNodeTaxon(hostNode).getId().substring(4));
				String taxonId = "symbiont" + i + "." + ++symbiontCounts[i - 1];
				node.setTaxon(new Taxon(taxonId));
				associations.put(taxonId,hostTree.getNodeTaxon(hostNode).getId());
				return node;
			}
			// Cospeciation event;
			child1 = simulateCoevolution(hostTree, hostTree.getChild(hostNode, 0), hostNodeHeight, rate, duplicationRate, hostShiftRate, lossRate, isRelaxed, stdev);
			child2 = simulateCoevolution(hostTree, hostTree.getChild(hostNode, 1), hostNodeHeight, rate, duplicationRate, hostShiftRate, lossRate, isRelaxed, stdev);
		} else {
			node.setHeight(eventHeight);
			switch(nextEvent.index) {
			case 0:
				// Duplication event
				child1 = simulateCoevolution(hostTree, hostNode, eventHeight, rate, duplicationRate, hostShiftRate, lossRate, isRelaxed, stdev);
				child2 = simulateCoevolution(hostTree, hostNode, eventHeight, rate, duplicationRate, hostShiftRate, lossRate, isRelaxed, stdev);
				break;
			case 1:
				// Host-shift event
				int nodeCount = hostTree.getNodeCount();
				NodeRef newHost;
				Relationship r;
				if (!hostTree.isRoot(hostNode)) { // Can't host-shift if at the root!
					do {
						newHost = hostTree
								.getNode(MathUtils.nextInt(nodeCount));
						r = Utils.determineRelationship(hostTree, hostNode,
								newHost).relationship;
					} while (hostTree.isRoot(newHost)
							|| (hostTree.getNodeHeight(newHost) >= eventHeight || eventHeight > hostTree // TODO Use of >= is debatable
									.getNodeHeight(hostTree.getParent(newHost))
									|| (r != Relationship.COUSIN || r != Relationship.SISTER)));
					child1 = simulateCoevolution(hostTree, newHost,
							eventHeight, rate, duplicationRate, hostShiftRate,
							lossRate, isRelaxed, stdev);
				} else {
					child1 = null; // Like a host-shift to a totally different tree that we're not following
				}
				child2 = simulateCoevolution(hostTree, hostNode, eventHeight, rate, duplicationRate, hostShiftRate, lossRate, isRelaxed, stdev);
				break;
			case 2: return null; // Loss event; null indicates the child lineage was lost
			default: throw new RuntimeException("Unknown cophylogenetic event: " + nextEvent.index); // Shouldn't be needed
			}
		}
		
		if (child1 != null && child2 != null) {
			// Both lineages survived
			node.addChild(child1);
			node.addChild(child2);
		} else if (child1 == child2) {
			// Both are null, hence this entire lineage was lost
			return null;
		} else {
			// Just one child lineage is null/lost
			return child2 == null ? child1 : child2; // Set this node to the continuing child lineage
		}
		return node;
	}
		
	private class EventIndexAndTime {
		public final int index;
		public final double time;
		public EventIndexAndTime(int index, double time) {
			this.index = index;
			this.time = time;
		}
	}
	
	private final EventIndexAndTime simulateSimultaneousPoissonProcesses(final double...lambdas) {
		final double lambda = MathUtils.getTotal(lambdas);
		final double[] p = new double[lambdas.length - 1];
		p[0] = lambdas[0] / lambda;
		for (int i = 1; i < p.length; ++i)
			p[i] = lambdas[i] / lambda + p[i - 1];
		final double time = MathUtils.nextExponential(lambda);
		final double U = 1 - MathUtils.nextDouble();
		int i;
		for (i = 0; i < p.length && p[i] < U; ++i);
		return new EventIndexAndTime(i, time);
	}
				
	public static void main(String[] args) {
		
		Locale.setDefault(Locale.US);
		Arguments arguments = new Arguments(new Arguments.Option[]{
				new Arguments.IntegerOption("t", "# taxa in host tree"),
				new Arguments.RealArrayOption("r", 3, "coevolutionary rates"),
				new Arguments.RealOption("c", "relaxed clock stdev"),
				new Arguments.LongOption("seed", "random number generator seed")
		});
		
		try {
			arguments.parseArguments(args);
		} catch (ArgumentException e) {
			e.printStackTrace(System.err);
			arguments.printUsage("coevolutionsim", "");
			System.exit(1);
		}
	
		long seed = MathUtils.getSeed();
		if (arguments.hasOption("seed")) {
            seed = arguments.getLongOption("seed");
            if (seed <= 0) {
                System.err.println("The random number seed should be > 0");
                System.exit(1);
            }
    		MathUtils.setSeed(seed);
        }
		System.out.println("Seed used: " + seed);
		System.out.println("Rates: " + ArrayUtils.toString(arguments.getRealArrayOption("r")));
		if (arguments.hasOption("c")) System.out.println(arguments.getRealOption("c"));
		System.out.println();
		
		final MutableTaxonList taxa = new Taxa();
		final int TAXA = arguments.getIntegerOption("t");
		for (int i = 1; i <= TAXA; ++i) taxa.addTaxon(new Taxon("host" + i));
		
		final Tree hostTree = new CoalescentSimulator().simulateTree(taxa, new ConstantPopulationModel(new Parameter.Default(100), Units.Type.YEARS));
		
		MutableTree mutableTree = (MutableTree) hostTree;
		for (int i = 0; i < mutableTree.getNodeCount(); ++i) {
			NodeRef node = mutableTree.getNode(i);
			mutableTree.setNodeAttribute(node, "nodeRef", node.getNumber());
		}
		
		PrintStream stream;
		File f;
		try {
			stream = new PrintStream(new FileOutputStream("host.tre"));
			new NexusExporter(stream).exportTree(hostTree);
			stream.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		f = new File("host.newick");
		f.deleteOnExit();
		try {
			stream = new PrintStream(new FileOutputStream(f));
			stream.println(Tree.Utils.newick(hostTree));
			stream.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		final double[] rates = arguments.getRealArrayOption("r");
		final SimpleCophylogenyModel scm = new SimpleCophylogenyModel(new Parameter.Default(rates[1]), new Parameter.Default(rates[2]), new Parameter.Default(1), Units.Type.YEARS);
		final CoevolutionSimulator cs = new CoevolutionSimulator();
		final Tree symbiontTree;
		if (arguments.hasOption("c")) {
			symbiontTree = cs.simulateCoevolution(hostTree, rates[0], scm, true, arguments.getRealOption("c"));
		} else {
			symbiontTree = cs.simulateCoevolution(hostTree, rates[0], scm, false);
		}
		
		try {
			stream = new PrintStream("symbiont.tre");
			new NexusExporter(stream).exportTree(symbiontTree);
			stream.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		f = new File("symbiont.newick");
		f.deleteOnExit();
		try {
			stream = new PrintStream(new FileOutputStream(f));
			stream.println(Tree.Utils.newick(symbiontTree));
			stream.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		try {
			stream = new PrintStream(new FileOutputStream("associations.map"));
			for (String key : cs.associations.keySet()) {
				stream.println(key + "\t" + cs.associations.get(key));
			}
			stream.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		SeqGen.main(new String[]{"host.newick", "host", "0.1"});
		SeqGen.main(new String[]{"symbiont.newick", "symbiont", "0.1"});
		
	}
	
	public static void debugHelper(Tree hostTree, Tree symbiontTree, CophylogenyLikelihood cophylogenyLikelihood) {
				
		System.err.println("host tree:");
		for (int i = 0; i < hostTree.getNodeCount(); i++) {
			NodeRef node = hostTree.getNode(i);
			if (!hostTree.isRoot(node)) {System.err.print(hostTree.getParent(node).getNumber() + ", ");} else {System.err.print(-1 + ", ");}
		}
		System.err.println();
		System.err.println("symbiont tree:");
		for (int i = 0; i < symbiontTree.getNodeCount(); i++) {
			NodeRef node = symbiontTree.getNode(i);
			if (!symbiontTree.isRoot(node)) System.err.print(symbiontTree.getParent(node).getNumber() + ", ");
		}
		System.err.println();
		System.err.println("associations:");
		for (int i = 0; i < symbiontTree.getNodeCount(); i++) {
			NodeRef node = symbiontTree.getNode(i);
			System.err.print(cophylogenyLikelihood.getStatesForNode(node).getNumber() + ", ");
		}
		System.err.println();
	}
	
}
