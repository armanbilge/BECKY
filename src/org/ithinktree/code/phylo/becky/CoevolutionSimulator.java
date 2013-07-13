/**
 * CoevolutionSimulator.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.code.phylo.becky;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.EnumSet;

import org.ithinktree.code.phylo.becky.CophylogenyModel.Relationship;

import dr.app.tools.NexusExporter;
import dr.app.util.Arguments;
import dr.app.util.Arguments.ArgumentException;
import dr.app.util.Arguments.IntegerOption;
import dr.app.util.Arguments.Option;
import dr.app.util.Arguments.RealArrayOption;
import dr.app.util.Arguments.StringOption;
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
 *
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
	public void simulateCoevolution(Tree hostTree, Tree symbiontTree, CophylogenyLikelihood cophylogenyLikelihood, String hostAttributeName) {
		
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
		int hostNodeCount = hostTree.getNodeCount();
		for (int i = 0; i < postOrderList.length; i++) {
			NodeRef node = symbiontTree.getNode(postOrderList[i]);
			if (!symbiontTree.isExternal(node)) {
				
				NodeRef child1HostNode = cophylogenyLikelihood.getStatesForNode(symbiontTree.getChild(node, 0));
				NodeRef child2HostNode = cophylogenyLikelihood.getStatesForNode(symbiontTree.getChild(node, 1));

				NodeRef hostNode = null;
				EnumSet<Relationship> relationships;
				int[] hostNodeIds = MathUtils.shuffled(hostNodeCount + 1);
				int j = 0;
				do {
					if (j > hostNodeCount) {
						i = -1; // Force the simulation to start over
						break;
					}
					relationships = EnumSet.noneOf(Relationship.class);
					int r = hostNodeIds[j++] - 1;
					hostNode = r == CophylogenyLikelihood.NO_HOST ? null : hostTree.getNode(r);
					relationships.add(Relationship.determineRelationship(hostTree, hostNode, child1HostNode).relationship);
					relationships.add(Relationship.determineRelationship(hostTree, hostNode, child2HostNode).relationship);
				} while (relationships.contains(Relationship.ANCESTOR) ||
							(relationships.contains(Relationship.SELF) && relationships.contains(Relationship.DESCENDANT)));
				cophylogenyLikelihood.setStatesForNode(node, hostNode);
			}
		}
	}
	
	public static Tree simulateCoevolution(Tree hostTree, SimpleCophylogenyModel model) {
		
		symbiontCounts = new int[hostTree.getTaxonCount()];
		
		SimpleNode root;
		do {
			root = simulateCoevolution(hostTree,
					hostTree.getRoot(),
					hostTree.getNodeHeight(hostTree.getRoot()),
					model.getDuplicationRate(),
					model.getHostShiftRate(),
					model.getLossRate());
		} while (root == null);
		return new SimpleTree(root);
	}
	
	private static int[] symbiontCounts;
	private static SimpleNode simulateCoevolution(Tree hostTree, NodeRef hostNode, double height, double duplicationRate, double hostShiftRate, double lossRate) {
				
		final SimpleNode node = new SimpleNode();
		node.setAttribute("host.nodeRef", hostNode.getNumber());
		
		SimpleNode child1 = null;
		SimpleNode child2 = null;
		
		final EventIndexAndTime nextEvent = simulateSimultaneousPoissonProcesses(duplicationRate, hostShiftRate, lossRate);
		final double eventHeight = height - nextEvent.time;
		
		final double hostNodeHeight = hostTree.getNodeHeight(hostNode);
		if (hostNodeHeight > eventHeight) {
			node.setHeight(hostTree.getNodeHeight(hostNode));
			if (hostTree.isExternal(hostNode)) {
				// Cannot coevolve anymore
				int i = hostNode.getNumber();
				node.setTaxon(new Taxon("symbiont" + i + "." + ++symbiontCounts[i]));
				return node;
			}
			// Cospeciation event;
			child1 = simulateCoevolution(hostTree, hostTree.getChild(hostNode, 0), hostNodeHeight, duplicationRate, hostShiftRate, lossRate);
			child2 = simulateCoevolution(hostTree, hostTree.getChild(hostNode, 1), hostNodeHeight, duplicationRate, hostShiftRate, lossRate);
		} else {
			node.setHeight(eventHeight);
			switch(nextEvent.index) {
			case 0:
				// Duplication event
				child1 = simulateCoevolution(hostTree, hostNode, eventHeight, duplicationRate, hostShiftRate, lossRate);
				child2 = simulateCoevolution(hostTree, hostNode, eventHeight, duplicationRate, hostShiftRate, lossRate);
				break;
			case 1:
				// Host-shift event
				int nodeCount = hostTree.getNodeCount();
				NodeRef newHost;
				Relationship r;
				do {
					newHost = hostTree.getNode(MathUtils.nextInt(nodeCount));
					r = Relationship.determineRelationship(hostTree, hostNode, newHost).relationship;
				} while (hostTree.getNodeHeight(newHost) >= eventHeight && eventHeight < hostTree.getNodeHeight(hostTree.getParent(newHost)) 
						&& (r != Relationship.COUSIN || r != Relationship.SISTER));
				child1 = simulateCoevolution(hostTree, newHost, eventHeight, duplicationRate, hostShiftRate, lossRate);
				child2 = simulateCoevolution(hostTree, hostNode, eventHeight, duplicationRate, hostShiftRate, lossRate);
			case 2: return null; // Loss event; null indicates the child linneage was lost
			default: throw new RuntimeException("Unknown cophylogenetic event: " + nextEvent.index); // Shouldn't be needed
			}
		}
		
		if (child1 != null && child2 != null) {
			// Both linneages survived
			node.addChild(child1);
			node.addChild(child2);
		} else if (child1 == child2) {
			// Both are null, hence this entire linneage was lost
			return null;
		} else {
			// Just one child linneage is null/lost
			return child2 == null ? child1 : child2; // Set this node to the continuing child linneage
		}
		return node;
	}
		
	private static class EventIndexAndTime {
		public final int index;
		public final double time;
		public EventIndexAndTime(int index, double time) {
			this.index = index;
			this.time = time;
		}
	}
	
	private static final double sum(double...values) {
		double sum = 0;
		for (double value : values) sum += value;
		return sum;
	}
	
	private static final EventIndexAndTime simulateSimultaneousPoissonProcesses(double...lambdas) {
		final double lambda = sum(lambdas);
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
		
		Arguments arguments = new Arguments(new Option[]{
				new StringOption("h", "filename", "symbiont tree file name"),
				new RealArrayOption("r", 3, "coevolutionary rates"),
				new StringOption("s","filename", "symbiont tree file name"),
				new IntegerOption("t", "# taxa in host tree")
		}, false);
		
		try {
			arguments.parseArguments(args);
		} catch (ArgumentException e) {
			e.printStackTrace(System.err);
			arguments.printUsage("coevolutionsim", "");
			System.exit(1);
		}
	
		MutableTaxonList taxa = new Taxa();
		final int TAXA = arguments.getIntegerOption("t");
		for (int i = 1; i <= TAXA; ++i) taxa.addTaxon(new Taxon("host" + i));
		
		final Tree hostTree = new CoalescentSimulator().simulateTree(taxa, new ConstantPopulationModel(new Parameter.Default(100), Units.Type.YEARS));
		
		MutableTree mutableTree = (MutableTree) hostTree;
		final double rootHeight = hostTree.getNodeHeight(hostTree.getRoot());
		for (int i = 0; i < mutableTree.getNodeCount(); i++) {
			NodeRef node = mutableTree.getNode(i);
			mutableTree.setNodeHeight(node, mutableTree.getNodeHeight(node) / rootHeight);
			mutableTree.setNodeAttribute(node, "nodeRef", node.getNumber());
		}
		
		PrintStream stream;
		try {
			stream = new PrintStream(new FileOutputStream(arguments.getStringOption("h")));
			new NexusExporter(stream).exportTree(hostTree);
			stream.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		final double[] rates = arguments.getRealArrayOption("r");
		final SimpleCophylogenyModel scm = new SimpleCophylogenyModel(new Parameter.Default(rates[0]), new Parameter.Default(rates[1]), new Parameter.Default(rates[2]), Units.Type.YEARS);
		final Tree symbiontTree = simulateCoevolution(hostTree, scm);
		
		try {
			stream = new PrintStream(new FileOutputStream(arguments.getStringOption("s")));
			new NexusExporter(stream).exportTree(symbiontTree);
			stream.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
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