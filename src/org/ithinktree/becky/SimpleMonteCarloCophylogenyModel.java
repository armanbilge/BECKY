package org.ithinktree.becky;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

@SuppressWarnings("serial")
public class SimpleMonteCarloCophylogenyModel extends SimpleCophylogenyModel {

	protected MonteCarloLikelihoodNoDescendants mc;
	
	public SimpleMonteCarloCophylogenyModel(Parameter duplicationRateParameter,
			Parameter hostSwitchRateParameter, Parameter lossRateParameter, int monteCarloIterations,
			Type units) {
		super(duplicationRateParameter, hostSwitchRateParameter, lossRateParameter,
				units);
		mc = new MonteCarloLikelihoodNoDescendants(this, monteCarloIterations);
	}

	 protected double likelihoodLineageLoss(final Tree tree, final NodeRef lineage, final double rate, boolean excludeRoot) {
		 if (excludeRoot) {
			 return mc.likelihoodNoDescendants(tree, lineage, tree.getNodeHeight(lineage), rate);
		 } else if (!tree.isRoot(lineage)) {
			 return mc.likelihoodNoDescendants(tree, lineage, tree.getNodeHeight(tree.getParent(lineage)), rate);
		 } else {
			 throw new RuntimeException("Not implemented");
		 }
	 }
	
}
