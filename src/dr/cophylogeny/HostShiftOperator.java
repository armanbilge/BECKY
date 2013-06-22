/**
 * HostShiftOperator.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package dr.cophylogeny;

import java.util.EnumSet;

import dr.cophylogeny.CophylogenyModel.Relationship;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * @author Arman D. Bilge
 *
 */
public class HostShiftOperator extends SimpleMCMCOperator {
	
	private Tree hostTree;
	private Tree symbiontTree;
	private CophylogenyLikelihood cophylogenyLikelihood;
	private int[] hostNodeIndices;
	
	/**
	 * 
	 */
	public HostShiftOperator(Tree hostTree, Tree symbiontTree, CophylogenyLikelihood cophylogenyLikelihood, double weight) {
		this.hostTree = hostTree;
		this.symbiontTree = symbiontTree;
		this.cophylogenyLikelihood = cophylogenyLikelihood;
		setWeight(weight);
		hostNodeIndices = new int[hostTree.getNodeCount() + 1];
		for (int i = 0; i < hostNodeIndices.length; i++) {
			hostNodeIndices[i] = i - 1; // Introduces a -1, aka CophylogenyLikelihood.NO_HOST
		}
		// Note that the MathUtils.shuffled() function will also initialize the array for us.
		//     I choose not to use it because I want a -1 in my array, as well as because this array
		//     will be reused, using their function doesn't really save me a step or efficiency.
		MathUtils.shuffle(hostNodeIndices);
	}

	@Override
	public String getPerformanceSuggestion() {
		return "No performance suggestion";
	}

	@Override
	public String getOperatorName() {
		return HostShiftOperatorParser.HOST_SHIFT_OPERATOR + "(" + symbiontTree.getId() + ")";
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		
		NodeRef node = symbiontTree.getInternalNode(MathUtils.nextInt(symbiontTree.getInternalNodeCount()));
		NodeRef child1HostNode = cophylogenyLikelihood.getStatesForNode(symbiontTree.getChild(node, 0));
		NodeRef child2HostNode = cophylogenyLikelihood.getStatesForNode(symbiontTree.getChild(node, 1));
		NodeRef hostNode;
		EnumSet<Relationship> relationships;
		MathUtils.shuffle(hostNodeIndices);
		int i = 0;
		do {
			relationships = EnumSet.noneOf(Relationship.class);
			hostNode = hostNodeIndices[i] == CophylogenyLikelihood.NO_HOST ? null : hostTree.getNode(hostNodeIndices[i]);
			if (hostNode != null) {
				relationships.add(Relationship.determineRelationship(hostTree, hostNode, child1HostNode).relationship);
				relationships.add(Relationship.determineRelationship(hostTree, hostNode, child2HostNode).relationship);
			}
			if (++i >= hostNodeIndices.length) break;
		} while (relationships.contains(Relationship.ANCESTOR) ||
					(relationships.contains(Relationship.SELF) && relationships.contains(Relationship.DESCENDANT)));
		
		// Check if loop ended with a valid operation
		if (!(relationships.contains(Relationship.ANCESTOR) ||
					(relationships.contains(Relationship.SELF) && relationships.contains(Relationship.DESCENDANT))))
			cophylogenyLikelihood.setStatesForNode(node, hostNode);
//		else
//			System.err.println("No valid host-shift operations!"); // Actually only at this node
		
		return 0;
	}

}
