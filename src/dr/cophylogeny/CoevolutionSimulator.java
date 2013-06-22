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
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
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
					if (j++ > hostNodeCount) {
						i = -1; // Force the simulation to start over
						break;
					}
					relationships = EnumSet.noneOf(Relationship.class);
					int r = hostNodeIds[j] - 1;
					hostNode = r == CophylogenyLikelihood.NO_HOST ? null : hostTree.getNode(r);
					relationships.add(Relationship.determineRelationship(hostTree, hostNode, child1HostNode).relationship);
					relationships.add(Relationship.determineRelationship(hostTree, hostNode, child2HostNode).relationship);
				} while (relationships.contains(Relationship.ANCESTOR) ||
							(relationships.contains(Relationship.SELF) && relationships.contains(Relationship.DESCENDANT)));
				cophylogenyLikelihood.setStatesForNode(node, hostNode);
			}
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
