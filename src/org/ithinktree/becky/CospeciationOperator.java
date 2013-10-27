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
    private final int internalNodeCount;
    
    public CospeciationOperator(final Tree hostTree, final MutableTree symbiontTree, final CophylogenyLikelihood cophylogenyLikelihood, final double weight) {
        this.hostTree = hostTree;
        this.symbiontTree = symbiontTree;
        this.cophylogenyLikelihood = cophylogenyLikelihood;
        internalNodeCount = symbiontTree.getInternalNodeCount();
        setWeight(weight);
    }
    
    @Override
    public String getPerformanceSuggestion() {
        return "No performance suggestion";
    }
    
    @Override
    public double doOperation() throws OperatorFailedException {
        
        final NodeRef node = symbiontTree.getNode(MathUtils.nextInt(internalNodeCount));
        final NodeRef host = cophylogenyLikelihood.getStatesForNode(node);
        if (host == null) throw new OperatorFailedException("No change in state");
        symbiontTree.setNodeHeight(node, hostTree.getNodeHeight(host));
        // TODO Need to determine correct Hastings ratio
        return 1.0;
    }
    
    @Override
    public String getOperatorName() {
        return CospeciationOperatorParser.COSPECIATION_OPERATOR + "(" + symbiontTree.getId() + ")";
    }
    
}
