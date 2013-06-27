/**
 * CoevolutionSimulator.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package dr.cophylogeny;

import java.util.EnumSet;

import dr.cophylogeny.CophylogenyModel.Relationship;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel.Node;
import dr.math.MathUtils;
import dr.math.distributions.ExponentialDistribution;

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
	
	public Tree simulateCoevolution(Tree hostTree, SimpleCophylogenyModel model) {
		
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
	
	private SimpleNode simulateCoevolution(Tree hostTree, NodeRef hostNode, double height, double duplicationRate, double hostShiftRate, double lossRate) {
				
		SimpleNode node = new SimpleNode();
		node.setHeight(height);
		node.setAttribute("host", hostNode);
		
		double nextDuplication = height - MathUtils.nextExponential(duplicationRate);
		double nextHostShift = height - MathUtils.nextExponential(hostShiftRate);
		double nextLoss = height - MathUtils.nextExponential(lossRate);
		
		SimpleNode child1 = null;
		SimpleNode child2 = null;
		
		double hostNodeHeight = hostTree.getNodeHeight(hostNode);
		if (hostNodeHeight > nextDuplication && hostNodeHeight > nextHostShift && hostNodeHeight > nextLoss) {
			if (hostTree.isExternal(hostNode)) {
				// Cannot coevolve anymore
				return node;
			}
			child1 = simulateCoevolution(hostTree, hostTree.getChild(hostNode, 0), hostNodeHeight, duplicationRate, hostShiftRate, lossRate);
			child2 = simulateCoevolution(hostTree, hostTree.getChild(hostNode, 1), hostNodeHeight, duplicationRate, hostShiftRate, lossRate);
		} else {
			switch(getMaxIndex(nextDuplication, nextHostShift, nextLoss)) {
			case 0:
				child1 = simulateCoevolution(hostTree, hostNode, nextDuplication, duplicationRate, hostShiftRate, lossRate);
				child2 = simulateCoevolution(hostTree, hostNode, nextDuplication, duplicationRate, hostShiftRate, lossRate);
				break;
			case 1:
				int nodeCount = hostTree.getNodeCount();
				NodeRef newHost;
				Relationship r;
				do {
					newHost = hostTree.getNode(MathUtils.nextInt(nodeCount));
					r = Relationship.determineRelationship(hostTree, hostNode, newHost).relationship;
				} while (hostTree.getNodeHeight(newHost) >= nextHostShift && nextHostShift < hostTree.getNodeHeight(hostTree.getParent(newHost)) 
						&& (r != Relationship.COUSIN || r != Relationship.SISTER));
				return simulateCoevolution(hostTree, newHost, nextHostShift, duplicationRate, hostShiftRate, lossRate);
			case 2: return null; // null indicates the child linneage was lost
			default: break;
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
	
	private static int getMaxIndex(double ... values) {
		double max = 0;
		int maxIndex = -1;
		for (int i = 0; i < values.length; ++i) {
			if (values[i] > max) maxIndex = i;
		}
		return maxIndex;
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
