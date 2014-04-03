/**
 * HostSwitchOperator.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky;

import java.util.List;

import org.ithinktree.becky.xml.HostSwitchOperatorParser;

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
public class HostSwitchOperator extends SimpleMCMCOperator {
	
	protected final Tree hostTree;
	protected final MutableTree symbiontTree;
	protected final CophylogenyLikelihood cophylogenyLikelihood;
	protected final boolean sampleNoHost;
	
	/**
	 * 
	 */
	public HostSwitchOperator(final Tree hostTree, final MutableTree symbiontTree, final CophylogenyLikelihood cophylogenyLikelihood, final boolean sampleNoHost, final double weight) {
		this.hostTree = hostTree;
		this.symbiontTree = symbiontTree;
		this.cophylogenyLikelihood = cophylogenyLikelihood;
		this.sampleNoHost = sampleNoHost;
		setWeight(weight);
	}

	@Override
	public String getPerformanceSuggestion() {
		return "No performance suggestion";
	}

	@Override
	public String getOperatorName() {
		return HostSwitchOperatorParser.HOST_SWITCH_OPERATOR + "(" + symbiontTree.getId() + ")";
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		
		NodeRef node = symbiontTree.getInternalNode(MathUtils.nextInt(symbiontTree.getInternalNodeCount()));
		final double nodeParentHeight = symbiontTree.isRoot(node) ? Double.POSITIVE_INFINITY : symbiontTree.getNodeHeight(symbiontTree.getParent(node));
		final double nodeChildHeight = Math.max(symbiontTree.getNodeHeight(symbiontTree.getChild(node, 0)), symbiontTree.getNodeHeight(symbiontTree.getChild(node, 1)));
		final List<NodeRef> hostNodes = CophylogenyModel.Utils.getLineagesInTimeRange(hostTree, nodeParentHeight, nodeChildHeight);
		final int i = sampleNoHost ? MathUtils.nextInt(hostNodes.size() + 1) - 1 : MathUtils.nextInt(hostNodes.size());
		final NodeRef proposedHost = i < 0 ? null : hostNodes.get(i);
		final NodeRef currentHost = cophylogenyLikelihood.getStatesForNode(node);
		if ((proposedHost != null && proposedHost.equals(currentHost))
				|| currentHost == null) throw new OperatorFailedException("No change of state");
		double hastingsRatio = 1.0;
		if (!CophylogenyModel.Utils.isContemporaneous(hostTree, proposedHost, symbiontTree.getNodeHeight(node))) {
		    final double min = Math.max(nodeChildHeight, hostTree.getNodeHeight(proposedHost));
		    double max = Math.min(nodeParentHeight, hostTree.isRoot(proposedHost) ? Double.POSITIVE_INFINITY : hostTree.getNodeHeight(hostTree.getParent(proposedHost)));
		    if (Double.isInfinite(max)) max = hostTree.getNodeHeight(proposedHost);
		    symbiontTree.setNodeHeight(node, MathUtils.nextDouble() * (max - min) + min);
		    final double hastingsMin = Math.max(nodeChildHeight, hostTree.getNodeHeight(currentHost));
		    double hastingsMax = Math.min(nodeParentHeight, hostTree.isRoot(currentHost) ? Double.POSITIVE_INFINITY : hostTree.getNodeHeight(hostTree.getParent(currentHost)));
		    if (Double.isInfinite(hastingsMax)) hastingsMax = hostTree.getNodeHeight(currentHost);
		    hastingsRatio = (hastingsMax - hastingsMin) / (max - min);
		}
		cophylogenyLikelihood.setStatesForNode(node, proposedHost);
		
		return hastingsRatio;
	}

}
