/**
 * 
 */
package org.ithinktree.becky.xml;

import java.io.IOException;

import org.ithinktree.becky.PreAnnotator;

import dr.evolution.io.Importer.ImportException;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author armanbilge
 *
 */
public class PreAnnotatorParser extends AbstractXMLObjectParser {

	private static final String PRE_ANNOTATOR = "preAnnotator";
	private static final String SYMBIONT_IN = "symbiontIn";
	private static final String SYMBIONT_OUT = "symbiontOut";
	private static final String HOST_IN = "hostIn";
	private static final String HOST_OUT = "hostOut";
	
	/* (non-Javadoc)
	 * @see dr.xml.XMLObjectParser#getParserName()
	 */
	@Override
	public String getParserName() {
		return PRE_ANNOTATOR;
	}

	/* (non-Javadoc)
	 * @see dr.xml.AbstractXMLObjectParser#getParserDescription()
	 */
	@Override
	public String getParserDescription() {
		return "A post-mcmc, pre-annotation tool.";
	}

	/* (non-Javadoc)
	 * @see dr.xml.AbstractXMLObjectParser#getReturnType()
	 */
	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return PreAnnotator.class;
	}

	private final XMLSyntaxRule[] rules = {AttributeRule.newStringRule(SYMBIONT_IN),
											AttributeRule.newStringRule(SYMBIONT_OUT),
											AttributeRule.newStringRule(HOST_IN),
											AttributeRule.newStringRule(HOST_OUT)};
	/* (non-Javadoc)
	 * @see dr.xml.AbstractXMLObjectParser#getSyntaxRules()
	 */
	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	/* (non-Javadoc)
	 * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
	 */
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		System.out.println("Initializing Cophylogeny PreAnnotator...");
		try {
			return new PreAnnotator(xo.getStringAttribute(HOST_IN), xo.getStringAttribute(SYMBIONT_IN), xo.getStringAttribute(HOST_OUT), xo.getStringAttribute(SYMBIONT_OUT));
		} catch (IOException e) {
			throw new XMLParseException(e.getMessage());
		} catch (ImportException e) {
			throw new XMLParseException(e.getMessage());
		}
	}

}
