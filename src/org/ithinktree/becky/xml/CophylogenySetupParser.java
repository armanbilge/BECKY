/**
 * CophylogenySetupParser.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */
package org.ithinktree.becky.xml;

import org.ithinktree.becky.CophylogenyLikelihood;
import org.ithinktree.becky.tools.CoevolutionSimulator;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Arman D. Bilge
 *
 */
public class CophylogenySetupParser extends AbstractXMLObjectParser {

	public static final String COPHYLOGENY_SETUP = "cophylogenySetup";
	public static final String HOST_TREE = "hostTree";
	public static final String SYMBIONT_TREE = "symbiontTree";

	
	@Override
	public String getParserName() {
		return COPHYLOGENY_SETUP;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
		XMLObject cxo = xo.getChild(HOST_TREE);
		final Tree hostTree = (Tree) cxo.getChild(Tree.class);
		
		cxo = xo.getChild(SYMBIONT_TREE);
		final MutableTree symbiontTree = (MutableTree) cxo.getChild(Tree.class);

		final CophylogenyLikelihood cl = (CophylogenyLikelihood) xo.getChild(CophylogenyLikelihood.class);
		
		for (int i = 0; i < symbiontTree.getNodeCount(); ++i) {
			NodeRef n = symbiontTree.getNode(i);
			cl.setStatesForNode(n, hostTree.getExternalNode(0));
		}
		cl.setStatesForNode(symbiontTree.getExternalNode(0), hostTree.getExternalNode(1));
		NodeRef n = symbiontTree.getRoot();
		do {
			n = Tree.Utils.postorderSuccessor(symbiontTree, n);
			if (!symbiontTree.isExternal(n) && !symbiontTree.isRoot(n)) {
				double d = Math.max(symbiontTree.getNodeHeight(symbiontTree.getChild(n, 0)), symbiontTree.getNodeHeight(symbiontTree.getChild(n, 1)));
				symbiontTree.setNodeHeight(n, (1.0 - MathUtils.nextDouble()) * (1.0 - d) + d);
			}
		} while (!symbiontTree.isRoot(n));
		cl.setStatesForNode(symbiontTree.getRoot(), hostTree.getRoot());
		symbiontTree.setNodeHeight(symbiontTree.getRoot(), 2.0);
		return null;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public String getParserDescription() {
		return "Simulates a symbiont starting tree by cospeciating on the host tree.";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return CoevolutionSimulator.class;
	}

	private final XMLSyntaxRule[] rules = {
			new ElementRule(HOST_TREE, new XMLSyntaxRule[]{
					new ElementRule(Tree.class)
			}),
			new ElementRule(SYMBIONT_TREE, new XMLSyntaxRule[]{
					new ElementRule(Tree.class)
			}),
			new ElementRule(CophylogenyLikelihood.class)
	};
}
