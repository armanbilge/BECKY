/**
 * HostShiftOperator.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky;

import java.util.EnumSet;

import org.ithinktree.becky.CophylogenyModel.Utils.Relationship;
import org.ithinktree.becky.xml.HostShiftOperatorParser;

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
	
	private final Tree hostTree;
	private final Tree symbiontTree;
	private final CophylogenyLikelihood cophylogenyLikelihood;
	private final int[] hostNodeIndices;
	
	/**
	 * 
	 */
	public HostShiftOperator(Tree hostTree, Tree symbiontTree, CophylogenyLikelihood cophylogenyLikelihood, boolean sampleNoHost, double weight) {
		this.hostTree = hostTree;
		this.symbiontTree = symbiontTree;
		this.cophylogenyLikelihood = cophylogenyLikelihood;
		setWeight(weight);
		if (sampleNoHost) {
			hostNodeIndices = new int[hostTree.getNodeCount() + 1];
			for (int i = 0; i < hostNodeIndices.length; i++) {
				hostNodeIndices[i] = i - 1; // Introduces a -1, aka CophylogenyLikelihood.NO_HOST
			}
			// Note that the MathUtils.shuffled() function will also initialize the array for us.
			//     I choose not to use it because I want a -1 in my array, as well as because this array
			//     will be reused, using their function doesn't really save me a step or efficiency.
			MathUtils.shuffle(hostNodeIndices);
		} else {
			hostNodeIndices = MathUtils.shuffled(hostTree.getNodeCount());
		}
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
		
		final NodeRef node = symbiontTree.getInternalNode(MathUtils.nextInt(symbiontTree.getInternalNodeCount()));
		final NodeRef parentHostNode = symbiontTree.isRoot(node) ? null : cophylogenyLikelihood.getStatesForNode(symbiontTree.getParent(node));
		final NodeRef child1HostNode = cophylogenyLikelihood.getStatesForNode(symbiontTree.getChild(node, 0));
		final NodeRef child2HostNode = cophylogenyLikelihood.getStatesForNode(symbiontTree.getChild(node, 1));
		final double nodeHeight = symbiontTree.getNodeHeight(node);
		NodeRef hostNode;
		EnumSet<Relationship> relationships;
		boolean temporallyValid;
		MathUtils.shuffle(hostNodeIndices);
		int i = 0;
		do {
			if (i >= hostNodeIndices.length) return 0;
			relationships = EnumSet.noneOf(Relationship.class);
			hostNode = hostNodeIndices[i] == CophylogenyLikelihood.NO_HOST ? null : hostTree.getNode(hostNodeIndices[i]);
			temporallyValid = (hostTree.isRoot(hostNode) ? true : (hostTree.getNodeHeight(hostTree.getParent(hostNode)) > nodeHeight)) && nodeHeight >= hostTree.getNodeHeight(hostNode);
			if (temporallyValid && hostNode != null) {
				relationships.add(CophylogenyModel.Utils.determineRelationship(hostTree, hostNode, child1HostNode).relationship);
				relationships.add(CophylogenyModel.Utils.determineRelationship(hostTree, hostNode, child2HostNode).relationship);
				if (parentHostNode != null)
					relationships.add(CophylogenyModel.Utils.determineRelationship(hostTree, parentHostNode, hostNode).relationship);
			}
			i++;
		} while (relationships.contains(Relationship.ANCESTOR) ||
					(relationships.contains(Relationship.SELF) && relationships.contains(Relationship.DESCENDANT)) || !temporallyValid);
		
		cophylogenyLikelihood.setStatesForNode(node, hostNode);
		
		return 0;
	}

}
