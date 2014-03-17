package org.ithinktree.becky;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class TipHostSwitchOperator extends HostSwitchOperator {
	
	public TipHostSwitchOperator(Tree hostTree, MutableTree symbiontTree,
			CophylogenyLikelihood cophylogenyLikelihood, boolean sampleNoHost,
			double weight) {
		super(hostTree, symbiontTree, cophylogenyLikelihood, sampleNoHost,
				weight);
	}

	public double doOperation() throws OperatorFailedException {
		
		final NodeRef node = symbiontTree.getExternalNode(MathUtils.nextInt(symbiontTree.getExternalNodeCount()));
		cophylogenyLikelihood.setStatesForNode(node, hostTree.getExternalNode(MathUtils.nextInt(hostTree.getExternalNodeCount())));
		return 1.0;
		
	}
	
}
