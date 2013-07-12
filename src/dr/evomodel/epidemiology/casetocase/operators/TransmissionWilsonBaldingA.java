package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.epidemiology.casetocase.AbstractCase;
import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.HashSet;

/**
 * Implements the Wilson-Balding branch swapping move if it does not change the transmission tree.
 *
 * @author Matthew Hall
 */

public class TransmissionWilsonBaldingA extends AbstractTreeOperator {

    private final CaseToCaseTreeLikelihood c2cLikelihood;
    public static final String TRANSMISSION_WILSON_BALDING_A = "transmissionWilsonBaldingA";
    private double logq;
    boolean debug = true;
    private final int tipCount;

    public TransmissionWilsonBaldingA(CaseToCaseTreeLikelihood c2cLikelihood, double weight) {
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
        tipCount = c2cLikelihood.getTree().getExternalNodeCount();
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        proposeTree();
        if (c2cLikelihood.getTree().getExternalNodeCount() != tipCount) {
            int newCount = c2cLikelihood.getTree().getExternalNodeCount();
            throw new RuntimeException("Lost some tips in modified SPR! (" +
                    tipCount + "-> " + newCount + ")");
        }

        return logq;
    }

    public void proposeTree() throws OperatorFailedException {
        TreeModel tree = c2cLikelihood.getTree();
        AbstractCase[] branchMap = c2cLikelihood.getBranchMap();
        NodeRef i;
        double oldMinAge, newMinAge, newRange, oldRange, newAge, q;
        // choose a random node avoiding root, and nodes that are ineligible for this move because they have nowhere to
        // go
        final int nodeCount = tree.getNodeCount();
        do {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
        } while (tree.getRoot() == i && !eligibleForMove(i, tree, branchMap));
        final NodeRef iP = tree.getParent(i);
        Integer[] samePaintings = c2cLikelihood.samePartition(iP, false);
        HashSet<Integer> possibleDestinations = new HashSet<Integer>();
        // we can insert the node above OR BELOW any node in the same partition
        for(int count=0; count<samePaintings.length; count++){
            possibleDestinations.add(samePaintings[count]);
            possibleDestinations.add(tree.getChild(tree.getNode(count),0).getNumber());
            possibleDestinations.add(tree.getChild(tree.getNode(count),1).getNumber());
        }
        Integer[] pd = possibleDestinations.toArray(new Integer[possibleDestinations.size()]);

        NodeRef j = tree.getNode(pd[MathUtils.nextInt(pd.length)]);
        NodeRef k = tree.getParent(j);

        while (k != null && (tree.getNodeHeight(k) <= tree.getNodeHeight(i)) || (i == j)) {
            j = tree.getNode(pd[MathUtils.nextInt(pd.length)]);
            k = tree.getParent(j);
        }

        if (iP == tree.getRoot()) {
            throw new OperatorFailedException("Root changes not allowed!");
        }

        if (k == iP || j == iP || k == i) throw new OperatorFailedException("move failed");

        final NodeRef CiP = getOtherChild(tree, iP, i);
        NodeRef PiP = tree.getParent(iP);

        newMinAge = Math.max(tree.getNodeHeight(i), tree.getNodeHeight(j));
        newRange = tree.getNodeHeight(k) - newMinAge;
        newAge = newMinAge + (MathUtils.nextDouble() * newRange);
        oldMinAge = Math.max(tree.getNodeHeight(i), tree.getNodeHeight(CiP));
        oldRange = tree.getNodeHeight(PiP) - oldMinAge;
        q = newRange / Math.abs(oldRange);

        if (j == tree.getRoot()) {

            // 1. remove edges <iP, CiP>
            tree.removeChild(iP, CiP);
            tree.removeChild(PiP, iP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(PiP, CiP);

            // iP is the new root
            tree.setRoot(iP);

        } else if (iP == tree.getRoot()) {

            // 1. remove edges <k, j>, <iP, CiP>, <PiP, iP>
            tree.removeChild(k, j);
            tree.removeChild(iP, CiP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(k, iP);

            //CiP is the new root
            tree.setRoot(CiP);

        } else {
            // 1. remove edges <k, j>, <iP, CiP>, <PiP, iP>
            tree.removeChild(k, j);
            tree.removeChild(iP, CiP);
            tree.removeChild(PiP, iP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(k, iP);
            tree.addChild(PiP, CiP);
        }

        if(debug){
            c2cLikelihood.checkPartitions();
        }
        logq = Math.log(q);

    }

    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    private boolean eligibleForMove(NodeRef node, TreeModel tree, AbstractCase[] branchMap){
        // to be eligible for this move, the node's parent and grandparent, or parent and other child, must be in the
        // same partition (so removing the parent has no effect on the transmission tree)

        return branchMap[tree.getParent(node).getNumber()]==branchMap[tree.getParent(tree.getParent(node)).getNumber()]
                || branchMap[tree.getParent(node).getNumber()]==branchMap[getOtherChild(tree,
                tree.getParent(node), node).getNumber()];
    }

    @Override
    public String getOperatorName() {
        return TRANSMISSION_WILSON_BALDING_A + " (" + c2cLikelihood.getTree().getId() +")";
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRANSMISSION_WILSON_BALDING_A;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final CaseToCaseTreeLikelihood c2cL
                    = (CaseToCaseTreeLikelihood) xo.getChild(CaseToCaseTreeLikelihood.class);
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            return new TransmissionWilsonBaldingA(c2cL, weight);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription(){
            return "This element represents a Wilson-Balding move operator, such that the transplantation of the " +
                    "subtree does not affect the topology of the transmission tree.";
        }

        public Class getReturnType(){
            return TransmissionWilsonBaldingA.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                    new ElementRule(CaseToCaseTreeLikelihood.class)
            };
        }
    };



}
