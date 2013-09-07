/**
 * CoevolutionSimulatorParser.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky.xml;

import org.ithinktree.becky.CophylogenyLikelihood;
import org.ithinktree.becky.tools.CoevolutionSimulator;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.Tree;
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
public class CoevolutionSimulatorParser extends AbstractXMLObjectParser {

	public static final String COEVOLUTION_SIMULATOR = "coevolutionSimulator";
	public static final String HOST_TREE = "hostTree";
	public static final String SYMBIONT_TREE = "symbiontTree";
	public static final String HOST_ATTRIBUTE_NAME = "hostAttributeName";
	public static final String SIMULATE_NO_HOST = "simulateNoHost";
	
	@Override
	public String getParserName() {
		return COEVOLUTION_SIMULATOR;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		final String hostAttributeName = xo.getStringAttribute(HOST_ATTRIBUTE_NAME);
		
		final boolean usingNoHost = xo.hasAttribute(SIMULATE_NO_HOST) ?
				xo.getBooleanAttribute(SIMULATE_NO_HOST) : false;
		
		XMLObject cxo = xo.getChild(HOST_TREE);
		final Tree hostTree = (Tree) cxo.getChild(Tree.class);
		
		cxo = xo.getChild(SYMBIONT_TREE);
		final MutableTree symbiontTree = (MutableTree) cxo.getChild(Tree.class);
		
		final CophylogenyLikelihood cophylogenyLikelihood = (CophylogenyLikelihood) xo.getChild(CophylogenyLikelihood.class);
		
		new CoevolutionSimulator().simulateCoevolution(hostTree, symbiontTree, cophylogenyLikelihood, hostAttributeName, usingNoHost);
		return null;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public String getParserDescription() {
		return "Simulates coevolution on the symbiont starting tree.";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return CoevolutionSimulator.class;
	}

	private final XMLSyntaxRule[] rules = {
			AttributeRule.newStringRule(HOST_ATTRIBUTE_NAME),
			AttributeRule.newBooleanRule(SIMULATE_NO_HOST, true),
			new ElementRule(HOST_TREE, new XMLSyntaxRule[] {
					new ElementRule(Tree.class)
			}),
			new ElementRule(SYMBIONT_TREE, new XMLSyntaxRule[] {
					new ElementRule(MutableTree.class)
			}),
			new ElementRule(CophylogenyLikelihood.class)
	};
}
