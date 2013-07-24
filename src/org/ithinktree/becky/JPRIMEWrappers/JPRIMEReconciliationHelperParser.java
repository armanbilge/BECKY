package org.ithinktree.becky.JPRIMEWrappers;

import se.cbb.jprime.apps.dlrs.ReconciliationHelper;
import se.cbb.jprime.topology.LeafLeafMap;
import se.cbb.jprime.topology.RBTree;
import se.cbb.jprime.topology.RBTreeEpochDiscretiser;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class JPRIMEReconciliationHelperParser extends AbstractXMLObjectParser {

	public static final String JPRIME_RECONCILIATION_HELPER = "jprimeReconciliationHelper";
	public static final String GUEST = "guest";
	public static final String HOST = "host";
	
	@Override
	public String getParserName() {
		return JPRIME_RECONCILIATION_HELPER;
	}

	@Override
	public String getParserDescription() {
		return "JPRIMEReconciliationHelper";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return ReconciliationHelper.class;
	}

	private final XMLSyntaxRule[] rules = {
			new ElementRule(GUEST, new XMLSyntaxRule[]{new ElementRule(RBTree.class)}),
			new ElementRule(HOST, new XMLSyntaxRule[]{new ElementRule(RBTree.class)}),
			new ElementRule(RBTreeEpochDiscretiser.class),
			new ElementRule(LeafLeafMap.class)
	};
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		XMLObject cxo = xo.getChild(GUEST);
		final RBTree guest = (RBTree) cxo.getChild(RBTree.class);
		
		cxo = xo.getChild(HOST);
		final RBTree host = (RBTree) cxo.getChild(RBTree.class);
		
		final RBTreeEpochDiscretiser rbTreeEpochDiscretiser;

		return null;
	}

}
