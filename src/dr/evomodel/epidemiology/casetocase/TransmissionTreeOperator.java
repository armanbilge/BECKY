package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.NodeRef;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * An operator that wraps a (phylogenetic) tree operator and adjusts the transmission tree accordingly. Only works
 * for the simpler (TT=TMRCA) version of the transmission tree. TT moves are deterministic, so the Hastings ratios
 * reported by the inner operators are not affected. Needs pretty serious documentation written at some point.
 *
 * @author Matthew Hall
 */

public class TransmissionTreeOperator extends AbstractCoercableOperator {

    private final CaseToCaseTransmissionLikelihood c2cLikelihood;
    private final AbstractTreeOperator innerOperator;
    public static final String TRANSMISSION_TREE_OPERATOR = "transmissionTreeOperator";

    public TransmissionTreeOperator(CaseToCaseTransmissionLikelihood c2cLikelihood, AbstractTreeOperator operator,
                                    CoercionMode mode) {
        super(mode);
        this.c2cLikelihood = c2cLikelihood;
        this.innerOperator = operator;
        this.setWeight(innerOperator.getWeight());
        if(c2cLikelihood.isExtended()){
            throw new RuntimeException("TransmissionTreeOperator only works on non-extended tree paintings");
        }
    }

    public TransmissionTreeOperator(CaseToCaseTransmissionLikelihood c2cLikelihood, AbstractTreeOperator operator) {
        this(c2cLikelihood,operator,CoercionMode.COERCION_OFF);
    }

    public double doOperation() throws OperatorFailedException {
        TreeModel tree = c2cLikelihood.getTree();
        AbstractCase[] branchMap = c2cLikelihood.getBranchMap();
        AbstractCase[] newBranchMap = Arrays.copyOf(branchMap, branchMap.length);
        int[] oldParents = getParentsArray(tree);
        double[] oldHeights = getHeightsArray(tree);
        double hr = innerOperator.doOperation();
        int[] newParents = getParentsArray(tree);
        ArrayList<Integer> changedNodes = new ArrayList<Integer>();
        for(int i=0; i<tree.getNodeCount(); i++){
            if(oldParents[i]!=newParents[i]){
                changedNodes.add(i);
            }
        }
        if(changedNodes.size()!=0){
            if(innerOperator instanceof ExchangeOperator){
                //this is a node swap operator
                AbstractCase[] nodePaintings = new AbstractCase[2];
                AbstractCase[] parentPaintings = new AbstractCase[2];
                for(int i=0; i<2; i++){
                    nodePaintings[i] = branchMap[changedNodes.get(i)];
                    parentPaintings[i] = branchMap[oldParents[changedNodes.get(i)]];
                }
                if(nodePaintings[0]==parentPaintings[0] || nodePaintings[1]==parentPaintings[1]){
                    //If this is not true there is nothing to do - the result is already a valid painting
                    for(int i=0; i<2; i++){
                        paintUp(tree, nodePaintings[1-i], nodePaintings[i], branchMap, newBranchMap,
                                tree.getNode(changedNodes.get(i)).getNumber(),newParents);
                    }
                }

            } else if(innerOperator instanceof SubtreeSlideOperator || innerOperator instanceof WilsonBalding){
                //this is a node transplantation operator
                int movedNode = -1;
                int oldChild = -1;
                int newChild = -1;
                for(int i=0; i<tree.getNodeCount(); i++){
                    //the chance of it having the same height as before is minuscule. If so, it will still throw an
                    // exception later, though.
                    if(tree.getNodeHeight(tree.getNode(i))!=oldHeights[i]){
                        movedNode = i;
                        break;
                    }
                }
                for(int j:changedNodes){
                    if(j!=movedNode){
                        if(tree.getParent(tree.getNode(j))==tree.getNode(movedNode)){
                            newChild = j;
                        } else {
                            oldChild = j;
                        }
                    }
                }
                if(movedNode == -1 || oldChild == -1 || newChild == -1){
                    throw new OperatorFailedException("Failed to establish relationship between relocated node and" +
                            " others");
                }
                NodeRef movedNodeObject = tree.getNode(movedNode);
                NodeRef oldChildObject = tree.getNode(oldChild);
                NodeRef newChildObject = tree.getNode(newChild);
                int otherChild = -1;
                NodeRef otherChildObject;
                //Find the other child of the moved node (the root of the transplanted subtree)
                for(int i=0; i<tree.getChildCount(movedNodeObject);i++){
                    if(tree.getChild(movedNodeObject,i)!=newChildObject){
                        otherChildObject = tree.getChild(movedNodeObject,i);
                        otherChild = otherChildObject.getNumber();
                    }
                }

                //If the child of the moved node is the earliest node with its painting:
                if(branchMap[otherChild]!=branchMap[movedNode]){
                    newBranchMap[movedNode]=branchMap[newChild];
                } else {
                    //Change all paintings up the tree from the old child that used to match the moved node to match
                    //the old child
                    paintUp(tree, branchMap[movedNode],branchMap[oldChild],branchMap,newBranchMap,oldChild,oldParents);
                    //This may have resulted in the moved node being recoloured wrong.
                    newBranchMap[movedNode]=branchMap[movedNode];
                    branchMap = newBranchMap;
                    //Change all paintings up the tree from the moved node that used to match the new child to match
                    //the moved node
                    paintUp(tree, branchMap[newChild],branchMap[movedNode],branchMap,newBranchMap,movedNode,newParents);
                }
            } else {
                //I don't know what this is
                throw new OperatorFailedException("Operator class "+innerOperator.getOperatorName()+" not yet " +
                        "supported");
            }
        }
        c2cLikelihood.setBranchMap(newBranchMap);
        c2cLikelihood.makeDirty();
        return hr;
    }




    private void paintUp(TreeModel tree, AbstractCase oldCase, AbstractCase newCase, AbstractCase[] oldbranchMap,
                         AbstractCase[] newBranchMap, int nodeNo, int[] parents){
        if(parents[nodeNo]==-1){
            return;
        }
        NodeRef newParent = tree.getNode(parents[nodeNo]);
        while(newParent!=null && oldbranchMap[newParent.getNumber()]==oldCase){
            newBranchMap[newParent.getNumber()]=newCase;
            newParent = parents[newParent.getNumber()]== -1 ? null : tree.getNode(parents[newParent.getNumber()]);
        }
    }

    private int[] getParentsArray(TreeModel tree){
        int[] out = new int[tree.getNodeCount()];
        for(int i=0; i<tree.getNodeCount(); i++){
            if(tree.getNode(i)==tree.getRoot()){
                out[i]=-1;
            } else {
                out[i]=tree.getParent(tree.getNode(i)).getNumber();
            }
        }
        return out;
    }

    private double[] getHeightsArray(TreeModel tree){
        double[] out = new double[tree.getNodeCount()];
        for(int i=0; i<tree.getNodeCount(); i++){
            out[i]=tree.getNodeHeight(tree.getNode(i));
        }
        return out;
    }

    public double getCoercableParameter() {
        if(innerOperator instanceof CoercableMCMCOperator){
            return ((CoercableMCMCOperator) innerOperator).getCoercableParameter();
        }
        throw new IllegalArgumentException();
    }

    public void setCoercableParameter(double value) {
        if(innerOperator instanceof CoercableMCMCOperator){
            ((CoercableMCMCOperator) innerOperator).setCoercableParameter(value);
            return;
        }
        throw new IllegalArgumentException();
    }

    public double getRawParameter() {
        if(innerOperator instanceof CoercableMCMCOperator){
            return ((CoercableMCMCOperator) innerOperator).getRawParameter();
        }
        throw new IllegalArgumentException();
    }

    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    public String getOperatorName() {
        return TRANSMISSION_TREE_OPERATOR;
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            CaseToCaseTransmissionLikelihood c2cL
                    = (CaseToCaseTransmissionLikelihood)xo.getChild(CaseToCaseTransmissionLikelihood.class);
            AbstractTreeOperator treeOp = (AbstractTreeOperator)xo.getChild(AbstractTreeOperator.class);

            CoercionMode mode = CoercionMode.COERCION_OFF;

            if(treeOp instanceof CoercableMCMCOperator){
                mode = ((CoercableMCMCOperator) treeOp).getMode();
            }

            return new TransmissionTreeOperator(c2cL,treeOp,mode);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        public String getParserDescription() {
            return "Performs a tree move then readjusts the transmission network in order to maintain its integrity.";
        }

        public Class getReturnType() {
            return TransmissionTreeOperator.class;
        }

        public String getParserName() {
            return TRANSMISSION_TREE_OPERATOR;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(CaseToCaseTransmissionLikelihood.class, "The transmission network likelihood element"),
                new ElementRule(AbstractTreeOperator.class, "A phylogenetic tree operator.")
        };
    };
}
