/**
 * NodeRefProviderParser.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package dr.cophylogeny;

import dr.evolution.tree.Tree;
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
public class NodeRefProviderParser extends AbstractXMLObjectParser {

	public static final String NODE_REF_PROVIDER = "nodeRefProvider";
	public static final String TAG_NAME = "tagName";

	@Override
	public String getParserName() {
		return NODE_REF_PROVIDER;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		final String tag = xo.getStringAttribute(TAG_NAME);
		final Tree tree = (Tree) xo.getChild(Tree.class);
		return new NodeRefProvider(tree, tag);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public String getParserDescription() {
		return "A utility element that logs nodes' integer values for post-inference reference.";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return NodeRefProvider.class;
	}

	private final XMLSyntaxRule[] rules = {
			AttributeRule.newStringRule(TAG_NAME),
			new ElementRule(Tree.class)
	};
	
}
