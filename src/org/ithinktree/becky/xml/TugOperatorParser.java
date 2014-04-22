package org.ithinktree.becky.xml;

import org.ithinktree.becky.CophylogenyLikelihood;
import org.ithinktree.becky.TugOperator;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.Tree;
import dr.inference.operators.MCMCOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class TugOperatorParser extends AbstractXMLObjectParser {

	public static final String TUG_OPERATOR = "tugOperator";
	public static final String HOST_TREE = "hostTree";
	public static final String SYMBIONT_TREE = "symbiontTree";

	
	@Override
	public String getParserName() {
		return TUG_OPERATOR;
	}

	@Override
	public String getParserDescription() {
		return null;
	}

	@Override
	public Class<?> getReturnType() {
		return TugOperator.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
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


	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
	
		final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
				
		XMLObject cxo = xo.getChild(HOST_TREE);
		final Tree hostTree = (Tree) cxo.getChild(Tree.class);
		
		cxo = xo.getChild(SYMBIONT_TREE);
		final MutableTree symbiontTree = (MutableTree) cxo.getChild(MutableTree.class);
		
		final CophylogenyLikelihood cophylogenyLikelihood = (CophylogenyLikelihood) xo.getChild(CophylogenyLikelihood.class);
		
		return new TugOperator(symbiontTree, hostTree, cophylogenyLikelihood, weight);
	
	}

}
