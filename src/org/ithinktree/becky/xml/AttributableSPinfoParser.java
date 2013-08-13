package org.ithinktree.becky.xml;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.util.Date;
import dr.evolution.util.Location;
import dr.evolution.util.Taxon;
import dr.evomodel.speciation.SpeciesBindings;
import dr.util.Attribute;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
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

	@SuppressWarnings("rawtypes")
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
		if (dr.xml.XMLParser.ID.contains("\'") && dr.xml.XMLParser.ID.contains("\"")) {
            // unable to handle taxon names that contain both single and double quotes
            // as it won't be possible to wrap it in either.
            throw new XMLParseException("Illegal taxon name, " + dr.xml.XMLParser.ID + ", - contains both single and double quotes");
        }
		
		List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0; i < xo.getChildCount(); ++i) {
        	Object cxo = xo.getChild(i);
        	if (cxo instanceof Taxon)
        		taxa.add((Taxon) cxo);
        }
        SpeciesBindings.SPinfo spInfo = new SpeciesBindings.SPinfo(xo.getId(), (Taxon[]) taxa.toArray());
        
        for (int i = 0; i < xo.getChildCount(); i++) {
            Object cxo = xo.getChild(i);

            if (cxo instanceof Date) {
                spInfo.setDate((Date) cxo);
            } else if (cxo instanceof Location) {
                spInfo.setLocation((Location) cxo);
            } else if (cxo instanceof Attribute) {
                final Attribute attr = (Attribute) cxo;
                spInfo.setAttribute(attr.getAttributeName(), attr.getAttributeValue());
            } else if (cxo instanceof Attribute[]) {
                Attribute[] attrs = (Attribute[]) cxo;
                for (Attribute attr : attrs) {
                    spInfo.setAttribute(attr.getAttributeName(), attr.getAttributeValue());
                }
            } else if (cxo instanceof Taxon) {
            	// Do nothing, already handled
            } else {
                throw new XMLParseException("Unrecognized element found in taxon element");
            }
        }
        
        return spInfo;
        
	}
	
	private final XMLSyntaxRule[] rules = {
			new StringAttributeRule(XMLParser.ID, "A unique identifier for this taxon"),
			new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
			new ElementRule(Attribute.Default.class, true),
            new ElementRule(Date.class, true),
            new ElementRule(Location.class, true)
	};

}
