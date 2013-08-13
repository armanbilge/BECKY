/**
 * 
 */
package org.ithinktree.becky.test;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import junit.framework.TestCase;

/**
 * @author armanbilge
 *
 */
public abstract class CophylogenyLikelihoodTest extends TestCase {

	/**
	 * 
	 */
	public CophylogenyLikelihoodTest() {
	}

	/**
	 * @param name
	 */
	public CophylogenyLikelihoodTest(String name) {
		super(name);
	}

	protected Tree hostTree;
	protected Tree symbiontTree;
	
	public void setUp() throws Exception {
		hostTree = new NewickImporter("tree;").importTree(null);
		symbiontTree = new NewickImporter("tree;").importTree(null);
	}
	
}
