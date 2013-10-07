/**
 * HostShiftOperator.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky;

import java.util.List;

import org.ithinktree.becky.xml.HostShiftOperatorParser;

import dr.evolution.tree.MutableTree;
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
	private final MutableTree symbiontTree;
	private final CophylogenyLikelihood cophylogenyLikelihood;
	private final int internalNodeCount;
	private final boolean sampleNoHost;
	
	/**
	 * 
	 */
	public HostShiftOperator(final Tree hostTree, final MutableTree symbiontTree, final CophylogenyLikelihood cophylogenyLikelihood, final boolean sampleNoHost, final double weight) {
		this.hostTree = hostTree;
		this.symbiontTree = symbiontTree;
		this.cophylogenyLikelihood = cophylogenyLikelihood;
		this.sampleNoHost = sampleNoHost;
		internalNodeCount = symbiontTree.getInternalNodeCount();
		setWeight(weight);
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
		
		final NodeRef node = symbiontTree.getInternalNode(MathUtils.nextInt(internalNodeCount));
		final double nodeParentHeight = symbiontTree.isRoot(node) ? Double.POSITIVE_INFINITY : symbiontTree.getNodeHeight(symbiontTree.getParent(node));
		final double nodeChildHeight = Math.max(symbiontTree.getNodeHeight(symbiontTree.getChild(node, 0)), symbiontTree.getNodeHeight(symbiontTree.getChild(node, 0)));
		final List<NodeRef> hostNodes = CophylogenyModel.Utils.getLineagesInTimeRange(hostTree, nodeParentHeight, nodeChildHeight);
		final int i = sampleNoHost ? MathUtils.nextInt(hostNodes.size() + 1) - 1 : MathUtils.nextInt(hostNodes.size());
		final NodeRef hostNode = i < 0 ? null : hostNodes.get(i);
		if ((hostNode != null && hostNode.equals(cophylogenyLikelihood.getStatesForNode(node)))
				|| cophylogenyLikelihood.getStatesForNode(node) == null) throw new OperatorFailedException("No change of state");
		if (!CophylogenyModel.Utils.isContemporaneous(hostTree, hostNode, symbiontTree.getNodeHeight(node))) {
		    final double min = Math.max(nodeChildHeight, hostTree.getNodeHeight(hostNode));
		    symbiontTree.setNodeHeight(node, MathUtils.nextDouble() * (Math.min(nodeParentHeight, hostTree.isRoot(hostNode) ? Double.POSITIVE_INFINITY : hostTree.getNodeHeight(hostTree.getParent(hostNode))) - min) + min);
		}
		cophylogenyLikelihood.setStatesForNode(node, hostNode);
		
		return 1.0; // I think that this is the correct hastings ratio, b/c all nodes have equal opportunity of being selected in both directions
	}

}
