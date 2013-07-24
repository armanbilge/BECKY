package org.ithinktree.becky;

import se.cbb.jprime.topology.RBTree;
import dr.evolution.tree.Tree;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class JPRIMERBTreeWrapperForBEASTTreeParser extends
		AbstractXMLObjectParser {

	public static final String JPRIME_RBTREE = "jprimeRBTree";
	
	@Override
	public String getParserName() {
		return JPRIME_RBTREE;
	}

	@Override
	public String getParserDescription() {
		return "Wrapper for the JPRIME_RBTREE";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return RBTree.class;
	}

	
	private final XMLSyntaxRule[] rules = {
			new ElementRule(Tree.class)
	};
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		return new JPRIMERBTreeWrapperForBEASTTree((Tree) xo.getChild(Tree.class));
	}

}
