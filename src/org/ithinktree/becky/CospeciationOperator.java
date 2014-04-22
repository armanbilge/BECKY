package org.ithinktree.becky;

import org.ithinktree.becky.xml.CospeciationOperatorParser;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;


public class CospeciationOperator extends SimpleMCMCOperator {
    
    private final Tree hostTree;
    private final MutableTree symbiontTree;
    private final CophylogenyLikelihood cophylogenyLikelihood;
    
    public CospeciationOperator(final Tree hostTree, final MutableTree symbiontTree, final CophylogenyLikelihood cophylogenyLikelihood, final double weight) {
        this.hostTree = hostTree;
        this.symbiontTree = symbiontTree;
        this.cophylogenyLikelihood = cophylogenyLikelihood;
        setWeight(weight);
    }
    
    @Override
    public String getPerformanceSuggestion() {
        return "No performance suggestion";
    }
    
    @Override
    public double doOperation() throws OperatorFailedException {
        
        final NodeRef node = symbiontTree.getInternalNode(MathUtils.nextInt(symbiontTree.getInternalNodeCount()));
        final NodeRef host = cophylogenyLikelihood.getStatesForNode(node);
        if (host == null || hostTree.isExternal(host)) throw new OperatorFailedException("No change in state");
        final double hostHeight = hostTree.getNodeHeight(host);
        if (!symbiontTree.isExternal(node) && hostHeight < Math.max(symbiontTree.getNodeHeight(symbiontTree.getChild(node, 0)), symbiontTree.getNodeHeight(symbiontTree.getChild(node, 1)))) throw new OperatorFailedException("No change in state");
//        final double nodeHeight = symbiontTree.getNodeHeight(node);
        final double maxHeight = Math.min(symbiontTree.isRoot(node) ? cophylogenyLikelihood.getOriginHeight() : symbiontTree.getNodeHeight(symbiontTree.getParent(node)), hostTree.isRoot(host) ? Double.POSITIVE_INFINITY : hostTree.getNodeHeight(hostTree.getParent(host)));
        final double range = maxHeight - hostHeight;
//        if (MachineAccuracy.same(nodeHeight, hostHeight)) {
//        	symbiontTree.setNodeHeight(node, MathUtils.nextDouble() * range + hostHeight);return Double.MAX_VALUE;
////        	return range / Double.MIN_NORMAL;
//        } else {
//        	symbiontTree.setNodeHeight(node, hostTree.getNodeHeight(host));
//        	return Double.MIN_NORMAL;// / range;
//        }
        if (MathUtils.nextInt(2) == 0) {
        	symbiontTree.setNodeHeight(node, hostHeight);
        } else {
        	symbiontTree.setNodeHeight(node, MathUtils.nextDouble() * range + hostHeight);
        }
        return 0.0;
    }
    
    @Override
    public String getOperatorName() {
        return CospeciationOperatorParser.COSPECIATION_OPERATOR + "(" + symbiontTree.getId() + ")";
    }
    
}
