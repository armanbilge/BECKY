/**
 * SimpleCophylogenyModelParser.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky.xml;

import org.ithinktree.becky.SimpleCophylogenyModel;

import dr.evolution.util.Units.Type;
import dr.evoxml.util.XMLUnits;
import dr.evoxml.util.XMLUnits.Utils;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Arman D. Bilge
 *
 */
public class SimpleCophylogenyModelParser extends AbstractXMLObjectParser {

	public static final String SIMPLE_COPHYLOGENY_MODEL = "simpleCophylogenyModel";
	public static final String DUPLICATION_RATE = "duplicationRate";
	public static final String HOST_SWITCH_RATE = "hostSwitchRate";
	public static final String LOSS_RATE = "lossRate";

	@Override
	public String getParserName() {
		return SIMPLE_COPHYLOGENY_MODEL;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
		final Type units = Utils.getUnitsAttr(xo);
				
		XMLObject cxo = xo.getChild(DUPLICATION_RATE);
		Parameter drParameter = (Parameter) cxo.getChild(Parameter.class);
		
		cxo = xo.getChild(HOST_SWITCH_RATE);
		Parameter hsrParameter = (Parameter) cxo.getChild(Parameter.class);
		
		cxo = xo.getChild(LOSS_RATE);
		Parameter lrParameter = (Parameter) cxo.getChild(Parameter.class);
		
		return new SimpleCophylogenyModel(drParameter, hsrParameter, lrParameter, units);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public String getParserDescription() {
		return "A simple model for cophylogenies.";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return SimpleCophylogenyModel.class;
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
			XMLUnits.SYNTAX_RULES[0]
	};
	
}
