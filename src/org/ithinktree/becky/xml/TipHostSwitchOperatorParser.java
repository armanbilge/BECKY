/**
 * HostSwitchOperator.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky.xml;

import org.ithinktree.becky.CophylogenyLikelihood;
import org.ithinktree.becky.HostSwitchOperator;
import org.ithinktree.becky.TipHostSwitchOperator;

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
public class TipHostSwitchOperatorParser extends AbstractXMLObjectParser {
	
	public static final String TIP_HOST_SWITCH_OPERATOR = "tipHostSwitchOperator";
	public static final String SAMPLE_NO_HOST = "sampleNoHost";
	public static final String HOST_TREE = "hostTree";
	public static final String SYMBIONT_TREE = "symbiontTree";

	@Override
	public String getParserName() {
		return TIP_HOST_SWITCH_OPERATOR;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
		final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
		
		final boolean usingNoHost = xo.hasAttribute(SAMPLE_NO_HOST) ?
							xo.getBooleanAttribute(SAMPLE_NO_HOST) : false;
		
		XMLObject cxo = xo.getChild(HOST_TREE);
		final Tree hostTree = (Tree) cxo.getChild(Tree.class);
		
		cxo = xo.getChild(SYMBIONT_TREE);
		final MutableTree symbiontTree = (MutableTree) cxo.getChild(MutableTree.class);
		
		final CophylogenyLikelihood cophylogenyLikelihood = (CophylogenyLikelihood) xo.getChild(CophylogenyLikelihood.class);
		
		return new TipHostSwitchOperator(hostTree, symbiontTree, cophylogenyLikelihood, usingNoHost, weight);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public String getParserDescription() {
		return "This basic operator switchs hosts on the symbiont tree while maintaining biological validity.";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return HostSwitchOperator.class;
	}

	private final XMLSyntaxRule[] rules = {
			AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
			AttributeRule.newBooleanRule(SAMPLE_NO_HOST, true),
			new ElementRule(HOST_TREE, new XMLSyntaxRule[]{
					new ElementRule(Tree.class)
			}),
			new ElementRule(SYMBIONT_TREE, new XMLSyntaxRule[]{
					new ElementRule(MutableTree.class)
			}),
			new ElementRule(CophylogenyLikelihood.class)
	};
	
}
