/**
 * CoevolutionSimulatorParser.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.code.phylo.becky;

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
	
	@Override
	public String getParserName() {
		return COEVOLUTION_SIMULATOR;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		CoevolutionSimulator simulator = new CoevolutionSimulator();

		final String hostAttributeName = xo.getStringAttribute(HOST_ATTRIBUTE_NAME);
		
		XMLObject cxo = xo.getChild(HOST_TREE);
		final Tree hostTree = (Tree) cxo.getChild(Tree.class);
		
		cxo = xo.getChild(SYMBIONT_TREE);
		final Tree symbiontTree = (Tree) cxo.getChild(Tree.class);
		
		final CophylogenyLikelihood cophylogenyLikelihood = (CophylogenyLikelihood) xo.getChild(CophylogenyLikelihood.class);
				
		simulator.simulateCoevolution(hostTree, symbiontTree, cophylogenyLikelihood, hostAttributeName);
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
			new ElementRule(HOST_TREE, new XMLSyntaxRule[] {
					new ElementRule(Tree.class)
			}),
			new ElementRule(SYMBIONT_TREE, new XMLSyntaxRule[] {
					new ElementRule(Tree.class)
			}),
			new ElementRule(CophylogenyLikelihood.class)
	};
}
