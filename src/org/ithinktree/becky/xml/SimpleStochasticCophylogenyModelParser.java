package org.ithinktree.becky.xml;

import org.ithinktree.becky.SimpleStochasticCophylogenyModel;

import dr.evolution.util.Units.Type;
import dr.evoxml.util.XMLUnits;
import dr.evoxml.util.XMLUnits.Utils;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class SimpleStochasticCophylogenyModelParser extends
		AbstractXMLObjectParser {

	public static final String SIMPLE_STOCHASTIC_COPHYLOGENY_MODEL = "simpleStochasticCophylogenyModel";
	public static final String DUPLICATION_RATE = "duplicationRate";
	public static final String HOST_SWITCH_RATE = "hostSwitchRate";
	public static final String LOSS_RATE = "lossRate";
	public static final String ITERATIONS = "iterations";

	
	@Override
	public String getParserName() {
		return SIMPLE_STOCHASTIC_COPHYLOGENY_MODEL;
	}

	@Override
	public String getParserDescription() {
		return "Parses a SimpleCophylogenyModel that employs a stochastic simulator.";
	}

	@Override
	public Class<?> getReturnType() {
		return SimpleStochasticCophylogenyModel.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	private final XMLSyntaxRule[] rules = {
			new ElementRule(DUPLICATION_RATE, new XMLSyntaxRule[]{
					new ElementRule(Parameter.class)
			}),
			new ElementRule(HOST_SWITCH_RATE, new XMLSyntaxRule[]{
					new ElementRule(Parameter.class)
			}),
			new ElementRule(LOSS_RATE, new XMLSyntaxRule[]{
					new ElementRule(Parameter.class)
			}),
			XMLUnits.SYNTAX_RULES[0],
			AttributeRule.newIntegerRule(ITERATIONS)
	};
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		final Type units = Utils.getUnitsAttr(xo);
		final int iterations = xo.getIntegerAttribute(ITERATIONS);
		
		XMLObject cxo = xo.getChild(DUPLICATION_RATE);
		final Parameter drParameter = (Parameter) cxo.getChild(Parameter.class);
		
		cxo = xo.getChild(HOST_SWITCH_RATE);
		final Parameter hsrParameter = (Parameter) cxo.getChild(Parameter.class);
		
		cxo = xo.getChild(LOSS_RATE);
		final Parameter lrParameter = (Parameter) cxo.getChild(Parameter.class);
		
		return new SimpleStochasticCophylogenyModel(drParameter, hsrParameter, lrParameter, iterations, units);
	}

}
