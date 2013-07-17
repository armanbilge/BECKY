/**
 * CophylogenyLikelihoodParser.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */

package org.ithinktree.becky;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser for the CophylogenyLikelihood Class.
 * 
 * @author Arman D. Bilge
 *
 */
public class CophylogenyLikelihoodParser extends AbstractXMLObjectParser {

	public static final String COPHYLOGENY_LIKELIHOOD = "cophylogenyLikelihood";
	public static final String MODEL = "model";
	public static final String HOST_TREE = "hostTree";
	public static final String SYMBIONT_TREE = "symbiontTree";
	public static final String CLOCK_MODEL = "clockModel";
	public static final String RECONSTRUCTION_TAG_NAME = "stateTagName";
	public static final String HOST_ATTRIBUTE_NAME = "hostAttributeName";

	
	@Override
	public String getParserName() {
		return COPHYLOGENY_LIKELIHOOD;
	}


	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
		final String reconstructionTagName = xo.getStringAttribute(RECONSTRUCTION_TAG_NAME);
		
		final String hostAttributeName = xo.getStringAttribute(HOST_ATTRIBUTE_NAME);
		
		final CophylogenyModel cophylogenyModel = (CophylogenyModel) xo.getChild(CophylogenyModel.class);
				
		XMLObject cxo = xo.getChild(HOST_TREE);
		final Tree hostTree = (Tree) cxo.getChild(Tree.class);
		
		cxo = xo.getChild(SYMBIONT_TREE);
		final Tree symbiontTree = (Tree) cxo.getChild(Tree.class);
		
		final BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
		
		return new CophylogenyLikelihood(hostTree, symbiontTree, cophylogenyModel, branchRateModel, reconstructionTagName, hostAttributeName, xo.getId());
	}


	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}


	@Override
	public String getParserDescription() {
		return "This element represents the likelihood of the cophylogenetic mapping given the host and symbiont trees.";
	}


	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return CophylogenyLikelihood.class;
	}

	private final XMLSyntaxRule[] rules = {
			AttributeRule.newStringRule(RECONSTRUCTION_TAG_NAME),
			AttributeRule.newStringRule(HOST_ATTRIBUTE_NAME),
			new ElementRule(CophylogenyModel.class),
			new ElementRule(HOST_TREE, new XMLSyntaxRule[]{
					new ElementRule(Tree.class)
			}),
			new ElementRule(SYMBIONT_TREE, new XMLSyntaxRule[]{
					new ElementRule(Tree.class)
			}),
			new ElementRule(BranchRateModel.class)
	};
	
}
