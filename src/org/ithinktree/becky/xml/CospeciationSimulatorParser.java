/**
 * CospeciationSimulatorParser.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky.xml;

import org.ithinktree.becky.SimpleCophylogenyModel;
import org.ithinktree.becky.tools.CoevolutionSimulator;

import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.inference.model.Parameter;
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
public class CospeciationSimulatorParser extends AbstractXMLObjectParser {

	public static final String COSPECIATION_SIMULATOR = "cospeciationSimulator";
	public static final String HOST_TREE = "hostTree";
	public static final String HOST_ATTRIBUTE_NAME = "hostAttributeName";
	public static final String SIMULATE_NO_HOST = "simulateNoHost";
	
	@Override
	public String getParserName() {
		return COSPECIATION_SIMULATOR;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
		XMLObject cxo = xo.getChild(HOST_TREE);
		final Tree hostTree = (Tree) cxo.getChild(Tree.class);
				
		return new CoevolutionSimulator().simulateCoevolution(hostTree, 0.0, new SimpleCophylogenyModel(new Parameter.Default(0.0),  new Parameter.Default(0.0), new Parameter.Default(0.0), Units.Type.YEARS), false);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public String getParserDescription() {
		return "Simulates a symbiont starting tree by cospeciating on the host tree.";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return CoevolutionSimulator.class;
	}

	private final XMLSyntaxRule[] rules = {
			AttributeRule.newBooleanRule(SIMULATE_NO_HOST, true),
			new ElementRule(HOST_TREE, new XMLSyntaxRule[] {
					new ElementRule(Tree.class)
			}),
	};
}
