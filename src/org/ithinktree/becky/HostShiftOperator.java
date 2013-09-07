/**
 * HostShiftOperator.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky;

import java.util.List;

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
	private final int internalNodeCount;
	private final boolean sampleNoHost;
	
	/**
	 * 
	 */
	public HostShiftOperator(final Tree hostTree, final Tree symbiontTree, final CophylogenyLikelihood cophylogenyLikelihood, final boolean sampleNoHost, final double weight) {
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
		final List<NodeRef> hostNodes = CophylogenyModel.Utils.getContemporaneousLineages(hostTree, symbiontTree.getNodeHeight(node));
		final int i = sampleNoHost ? MathUtils.nextInt(hostNodes.size() + 1) - 1 : MathUtils.nextInt(hostNodes.size());
		final NodeRef hostNode = i < 0 ? null : hostNodes.get(i);
		if ((hostNode != null && hostNode.equals(cophylogenyLikelihood.getStatesForNode(node)))
				|| cophylogenyLikelihood.getStatesForNode(node) == null) throw new OperatorFailedException("No change of state");
		cophylogenyLikelihood.setStatesForNode(node, hostNode);
		
		return 1.0; // I think that this is the correct hastings ratio, b/c all nodes have equal opportunity of being selected in both directions
	}

}
