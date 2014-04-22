package org.ithinktree.becky;

import org.ithinktree.becky.CophylogenyModel.Utils.NodalRelationship;
import org.ithinktree.becky.CophylogenyModel.Utils.Relationship;
import org.ithinktree.becky.xml.TugOperatorParser;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public class TugOperator extends SimpleMCMCOperator {

	private final MutableTree symbiontTree;
	private final Tree hostTree;
	private final CophylogenyLikelihood cophylogenyLikelihood;
	
	public TugOperator(final MutableTree symbiontTree, final Tree hostTree, final CophylogenyLikelihood cophylogenyLikelihood, final double weight) {
		this.symbiontTree = symbiontTree;
		this.hostTree = hostTree;
		this.cophylogenyLikelihood = cophylogenyLikelihood;
		setWeight(weight);	
	}
	
	@Override
	public String getPerformanceSuggestion() {
		return "No performance suggestion";
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		
//		System.err.println("BEFORE: " + Tree.Utils.newick(symbiontTree, new TreeTraitProvider[]{cophylogenyLikelihood}));
		final NodeRef self = symbiontTree.getInternalNode(MathUtils.nextInt(symbiontTree.getInternalNodeCount()));
		final NodeRef child1 = symbiontTree.getChild(self, 0);
		final NodeRef child2 = symbiontTree.getChild(self, 1);
		final NodeRef selfHost = cophylogenyLikelihood.getStatesForNode(self);
		final NodeRef child1Host = cophylogenyLikelihood.getStatesForNode(child1);
		final NodeRef child2Host = cophylogenyLikelihood.getStatesForNode(child2);
		if (self == null || child1 == null || child2 == null || symbiontTree.isExternal(child1) || symbiontTree.isExternal(child2)) throw new OperatorFailedException("No change in state");
		final NodalRelationship rel1 = CophylogenyModel.Utils.determineRelationship(hostTree, selfHost, child1Host);
		final NodalRelationship rel2 = CophylogenyModel.Utils.determineRelationship(hostTree, selfHost, child2Host);
		
		int i = MathUtils.nextInt(2);
		if (i == 0 && rel1.relationship == Relationship.SELF && rel2.relationship == Relationship.SELF) {
			final NodeRef host = cophylogenyLikelihood.getStatesForNode(self);
			final double hostHeight = hostTree.getNodeHeight(host);

			final NodeRef hostChild1 = hostTree.getChild(host, 0);
			final NodeRef hostChild2 = hostTree.getChild(host, 1);
			cophylogenyLikelihood.setStatesForNode(child1, hostChild1);
			cophylogenyLikelihood.setStatesForNode(child2, hostChild2);
			
	        double minHeight = Math.max(Math.max(symbiontTree.getNodeHeight(symbiontTree.getChild(child1, 0)), symbiontTree.getNodeHeight(symbiontTree.getChild(child1, 1))), hostTree.isExternal(hostChild1) ? 0 : hostTree.getNodeHeight(hostChild1));
	        final double range1 = hostHeight - minHeight;
			symbiontTree.setNodeHeight(child1, MathUtils.nextDouble() * range1 + minHeight);

	        minHeight = Math.max(Math.max(symbiontTree.getNodeHeight(symbiontTree.getChild(child2, 0)), symbiontTree.getNodeHeight(symbiontTree.getChild(child2, 1))), hostTree.isExternal(hostChild2) ? 0 : hostTree.getNodeHeight(hostChild2));
	        final double range2 = hostHeight - minHeight;
			symbiontTree.setNodeHeight(child2, MathUtils.nextDouble() * range2 + minHeight);
			
	        final double maxHeight = Math.min(symbiontTree.isRoot(self) ? cophylogenyLikelihood.getOriginHeight() : symbiontTree.getNodeHeight(symbiontTree.getParent(self)), hostTree.isRoot(host) ? Double.POSITIVE_INFINITY : hostTree.getNodeHeight(hostTree.getParent(host)));
	        final double range3 = maxHeight - hostHeight;
			
			symbiontTree.setNodeHeight(self, hostHeight);

			return Math.log(range1 * range2 / range3);
			
		} else if (i == 1 && rel1.relationship == Relationship.DESCENDANT && rel2.relationship == Relationship.DESCENDANT) {
			final NodeRef host = cophylogenyLikelihood.getStatesForNode(self);
			final double hostHeight = hostTree.getNodeHeight(host);

			final NodeRef hostChild1 = cophylogenyLikelihood.getStatesForNode(child1);
			final NodeRef hostChild2 = cophylogenyLikelihood.getStatesForNode(child2);
			cophylogenyLikelihood.setStatesForNode(child1, host);
			cophylogenyLikelihood.setStatesForNode(child2, host);
			symbiontTree.setNodeHeight(child1, hostHeight);
			symbiontTree.setNodeHeight(child2, hostHeight);
	        final double maxHeight = Math.min(symbiontTree.isRoot(self) ? cophylogenyLikelihood.getOriginHeight() : symbiontTree.getNodeHeight(symbiontTree.getParent(self)), hostTree.isRoot(host) ? Double.POSITIVE_INFINITY : hostTree.getNodeHeight(hostTree.getParent(host)));
	        final double range3 = maxHeight - hostHeight;
			symbiontTree.setNodeHeight(self, MathUtils.nextDouble() * range3 + hostHeight);
			
	        double minHeight = Math.max(Math.max(symbiontTree.getNodeHeight(symbiontTree.getChild(child1, 0)), symbiontTree.getNodeHeight(symbiontTree.getChild(child1, 1))), hostTree.isExternal(hostChild1) ? 0 : hostTree.getNodeHeight(hostChild1));
	        final double range1 = hostHeight - minHeight;
			symbiontTree.setNodeHeight(child1, hostHeight);

	        minHeight = Math.max(Math.max(symbiontTree.getNodeHeight(symbiontTree.getChild(child2, 0)), symbiontTree.getNodeHeight(symbiontTree.getChild(child2, 1))), hostTree.isExternal(hostChild2) ? 0 : hostTree.getNodeHeight(hostChild2));
	        final double range2 = hostHeight - minHeight;
			symbiontTree.setNodeHeight(child2, hostHeight);

			return Math.log(range3 / (range1 * range2));
		} else {
			throw new OperatorFailedException("No change in state");
		}
		
	}

	@Override
	public String getOperatorName() {
		return TugOperatorParser.TUG_OPERATOR + "(" + symbiontTree.getId() + ")";
	}

}
