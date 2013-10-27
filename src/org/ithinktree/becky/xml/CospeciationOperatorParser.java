/**
 * HostSwitchOperator.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky.xml;

import org.ithinktree.becky.CophylogenyLikelihood;
import org.ithinktree.becky.CospeciationOperator;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.Tree;
import dr.inference.operators.MCMCOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Arman D. Bilge
 *
 */
public class CospeciationOperatorParser extends AbstractXMLObjectParser {
    
    public static final String COSPECIATION_OPERATOR = "cospeciationOperator";
    public static final String HOST_TREE = "hostTree";
    public static final String SYMBIONT_TREE = "symbiontTree";

    @Override
    public String getParserName() {
        return COSPECIATION_OPERATOR;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        
         XMLObject cxo = xo.getChild(HOST_TREE);
        final Tree hostTree = (Tree) cxo.getChild(Tree.class);
        
        cxo = xo.getChild(SYMBIONT_TREE);
        final MutableTree symbiontTree = (MutableTree) cxo.getChild(MutableTree.class);
        
        final CophylogenyLikelihood cophylogenyLikelihood = (CophylogenyLikelihood) xo.getChild(CophylogenyLikelihood.class);
        
        return new CospeciationOperator(hostTree, symbiontTree, cophylogenyLikelihood, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "This operator creates true cospeciation events.";
    }

    @SuppressWarnings("rawtypes")
    @Override
    public Class getReturnType() {
        return CospeciationOperator.class;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(HOST_TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }),
            new ElementRule(SYMBIONT_TREE, new XMLSyntaxRule[]{
                    new ElementRule(MutableTree.class)
            }),
            new ElementRule(CophylogenyLikelihood.class)
    };
    
}
