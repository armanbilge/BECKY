package org.ithinktree.code.phylo.becky;

import dr.evolution.util.Taxon;
import dr.evomodel.speciation.SpeciesBindings;
import dr.util.Attribute;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class AttributableSPinfoParser extends AbstractXMLObjectParser {

	public static final String ATTR_SP = "attrSp";
	
	@Override
	public String getParserName() {
		return ATTR_SP;
	}

	@Override
	public String getParserDescription() {
		return "Taxon in a species tree";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return SpeciesBindings.SPinfo.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public Object parseXMLObject(XMLObject arg0) throws XMLParseException {
		// TODO Auto-generated method stub
		return null;
	}
	
	private final XMLSyntaxRule[] rules = {
			new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
			new ElementRule(Attribute.Default.class, true)
	};

}
