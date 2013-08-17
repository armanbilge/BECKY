/**
 * 
 */
package org.ithinktree.becky.jprimewrappers.xml;

import org.ithinktree.becky.jprimewrappers.CophylogenyLikelihoodWrapperForJPrIMEDLTRSModel;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
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
public class CophylogenyLikelihoodWrapperForJPrIMEDLTRSModelParser extends
		AbstractXMLObjectParser {

	private static final String DLTRS_MODEL = "dltrsModel";
	private static final String HOST_TREE = "hostTree";
	private static final String GUEST_TREE = "guestTree";
	private static final String DUPLICATION_RATE = "duplicationRate";
	private static final String HOST_SHIFT_RATE = "hostShiftRate";
	private static final String LOSS_RATE = "lossRate";
	private static final String NORMALIZE = "normalize";
	private static final String HOST_ATTRIBUTE_NAME = "hostAttributeName";
	private static final String MEAN = "mean";
	private static final String STDEV = "stdev";

	@Override
	public String getParserName() {
		return DLTRS_MODEL;
	}

	@Override
	public String getParserDescription() {
		return "Wrapper for the JPrIME DLTRS Model";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return CophylogenyLikelihoodWrapperForJPrIMEDLTRSModel.class;
	}

	private final XMLSyntaxRule[] rules = {
			new ElementRule(Taxa.class),
			new ElementRule(HOST_TREE, new XMLSyntaxRule[]{new ElementRule(Tree.class)}),
			new ElementRule(GUEST_TREE, new XMLSyntaxRule[]{new ElementRule(Tree.class)}),
			new ElementRule(DUPLICATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
			new ElementRule(HOST_SHIFT_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
			new ElementRule(LOSS_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
			new ElementRule(MEAN, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
			new ElementRule(STDEV, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
			AttributeRule.newBooleanRule(NORMALIZE),
			AttributeRule.newStringRule(HOST_ATTRIBUTE_NAME),
	};
	
	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		final Tree hostTree = (Tree) xo.getChild(HOST_TREE).getChild(Tree.class);
		final Tree guestTree = (Tree) xo.getChild(GUEST_TREE).getChild(Tree.class);
		
		final Taxa taxa = (Taxa) xo.getChild(Taxa.class);
		
		final Parameter duplicationRate = (Parameter) xo.getChild(DUPLICATION_RATE).getChild(Parameter.class);
		final Parameter hostShiftRate = (Parameter) xo.getChild(HOST_SHIFT_RATE).getChild(Parameter.class);
		final Parameter lossRate = (Parameter) xo.getChild(LOSS_RATE).getChild(Parameter.class);
		final Parameter mean = (Parameter) xo.getChild(MEAN).getChild(Parameter.class);
		final Parameter stdev = (Parameter) xo.getChild(STDEV).getChild(Parameter.class);
		
		final boolean normalize = xo.getBooleanAttribute(NORMALIZE);
		final String hostAttributeName = xo.getStringAttribute(HOST_ATTRIBUTE_NAME);
		
		return new CophylogenyLikelihoodWrapperForJPrIMEDLTRSModel(xo.getId(), hostTree, guestTree, duplicationRate, hostShiftRate, lossRate, normalize, taxa, mean, stdev, hostAttributeName);
	}

}
