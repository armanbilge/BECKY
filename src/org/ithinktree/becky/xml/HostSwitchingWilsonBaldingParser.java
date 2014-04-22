package org.ithinktree.becky.xml;

import org.ithinktree.becky.CophylogenyLikelihood;
import org.ithinktree.becky.HostSwitchingWilsonBalding;

import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 */
public class HostSwitchingWilsonBaldingParser extends AbstractXMLObjectParser {

    public static final String HOST_SWITCHING_WILSON_BALDING = "hostSwitchingWilsonBalding";

    public String getParserName() {
        return HOST_SWITCHING_WILSON_BALDING;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final CophylogenyLikelihood cl = (CophylogenyLikelihood) xo.getChild(CophylogenyLikelihood.class);

        return new HostSwitchingWilsonBalding(treeModel, cl, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(TreeModel.class),
            new ElementRule(CophylogenyLikelihood.class)
    };

    public String getParserDescription() {
        return "An operator which performs the host-switching Wilson-Balding move on a tree";
    }

    public Class<?> getReturnType() {
        return HostSwitchingWilsonBalding.class;
    }
}
