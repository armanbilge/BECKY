/**
 * CophylogenyModelTest.java
 *
 * BECKY
 */
package test.org.ithinktree.becky;

import java.util.Arrays;
import java.util.List;

import org.ithinktree.becky.CophylogenyModel;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * @author Arman D. Bilge
 *
 */
@RunWith(JUnit4.class)
public class CophylogenyModelUtilsTest {

	private Tree tree = TestUtils.DEFAULT_TREE;
	
	@Test
	public void testGetContemporaneousLineages() {
		boolean test = true;
		
		List<NodeRef> lineages = CophylogenyModel.Utils.getContemporaneousLineages(tree, tree.getNodeHeight(TestUtils.CDE));
		
		test &= lineages.contains(TestUtils.CDE);
		test &= lineages.contains(TestUtils.B);
		test &= lineages.contains(TestUtils.A);
		
		Assert.assertTrue(test);
	}

	@Test
	public void testGetContemporaneousLineageCount() {
		Assert.assertTrue(CophylogenyModel.Utils.getContemporaneousLineageCount(tree, tree.getNodeHeight(TestUtils.CDE)) ==
				CophylogenyModel.Utils.getContemporaneousLineages(tree, tree.getNodeHeight(TestUtils.CDE)).size());
	}
	
	@Test
	public void testLostLineages() {
	    
	    List<NodeRef> determineRelationshipLostLineages = Arrays.asList(CophylogenyModel.Utils.determineRelationship(tree, TestUtils.DE, TestUtils.BCDE).lostLineages);
	    List<NodeRef> lostLineagesToTime = Arrays.asList(CophylogenyModel.Utils.lostLineagesToTime(tree, TestUtils.DE, tree.getNodeHeight(TestUtils.BCDE)));
	    
	    NodeRef[] knownLostLineages = new NodeRef[]{TestUtils.B, TestUtils.C};
	    Assert.assertTrue(knownLostLineages.length == determineRelationshipLostLineages.size());
	    Assert.assertTrue(knownLostLineages.length == lostLineagesToTime.size());
	    for (NodeRef n : knownLostLineages) {
	        Assert.assertTrue(determineRelationshipLostLineages.contains(n));
	        Assert.assertTrue(lostLineagesToTime.contains(n));
	    }
	    
	}
	
}
