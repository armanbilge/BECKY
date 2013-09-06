/**
 * SimpleCophylogenyModelTest.java
 *
 * BECKY
 */
package test.org.ithinktree.becky;

import java.lang.reflect.Method;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.ithinktree.becky.CophylogenyLikelihood;
import org.ithinktree.becky.SimpleCophylogenyModel;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.inference.model.Parameter;
import dr.math.MachineAccuracy;

/**
 * @author Arman D. Bilge
 *
 */
@RunWith(JUnit4.class)
public class SimpleCophylogenyModelTest {

	private Parameter duplicationRate;
	private Parameter hostShiftRate;
	private Parameter lossRate;
	private BranchRates branchRates;
	private SimpleCophylogenyModel model;
	private Tree host;
	
	@Before
	public void before() {
		
		duplicationRate = new Parameter.Default(1.5);
		hostShiftRate = new Parameter.Default(0.5);
		lossRate = new Parameter.Default(1.0);
		model = new SimpleCophylogenyModel(duplicationRate, hostShiftRate, lossRate, Units.Type.YEARS);
		
		// Force internal variables to update
		try {
			model.calculateNodeLogLikelihood(null, null, null, null, null, null, null, null, null);
		} catch (Exception e) {
			// Expected; do nothing
		}
		
		branchRates = new BranchRates() {
			public double getBranchRate(Tree t, NodeRef n) {
				return 1.0;
			}
		};
		
		host = TestUtils.DEFAULT_TREE;
	}
		
	private CophylogenyLikelihood createCophylogenyLikelihood(final Tree symbiont) {
		return new CophylogenyLikelihood(host, (MutableTree) symbiont, model, branchRates, "host.nodeRef", "testCophylogenyLikelihood");
	}
	
	@Test
	public void testSelfSelf() {
		// Construct the symbiont tree
//		final String newickTree = "(A,B);";
//		final Tree symbiont = TestUtils.treeFromNewick(newickTree, true);
//		CophylogenyLikelihood cl = createCophylogenyLikelihood(symbiont);
//		Assert.assertEquals(0.0, cl.getLogLikelihood(), MachineAccuracy.EPSILON);
	}
	
	@Test
	public void testSelfCousin() {
		
	}
	
	@Test
	public void testCousinCousin() {
		
	}
	
	@Test
	public void testDescendantCousin() {
		
	}
	
	@Test
	public void testDescendantDescendantCousin() {
		
	}
	
	@Test
	public void testDescendantDescendantDirect() {
		
	}
	
	@Test
	public void testLikelihoodHostShiftEventAndLossInTimeIntegration() {
				
		Method m;
		try {
			m = model.getClass().getDeclaredMethod("likelihoodHostShiftEventAndLossInTime", double.class, double.class, double.class, double.class, double.class, double.class);
		} catch (Exception e) {
			Assert.fail("Fatal reflection error retrieving method: " + e.toString());
			return;
		}
		m.setAccessible(true);
		double likelihood;
		try {
			likelihood = (Double) m.invoke(model, 0.1, 0.2, 0.27, 0.23, 1.5, 1.0);
		} catch (Exception e) {
			Assert.fail("Fatal reflection error invoking method: " + e.toString());
			return;
		}
		
		Assert.assertEquals(0.0023401643365404423, likelihood, 1E-19);
	}
	
	@Test
	public void testLikelihoodLineageLoss() {
	    
	    final Tree tree = TestUtils.treeFromNewick("((A:1.0,B:1.0):2.0,C:3.0);", true);
	    final double overallRate = model.getOverallRate();
	    final double lossRate = model.getLossRate();
	    final double actualLikelihood = 0.0 + Math.exp(0.0 * overallRate) * (1 - Math.exp(3 * lossRate)) * ((1 - Math.exp(2 * lossRate) + Math.exp(2 * overallRate) * (1 - Math.exp(1 * lossRate) * (1 - Math.exp(1 * lossRate)))));
	    
	    Method m;
        try {
            m = model.getClass().getDeclaredMethod("likelihoodLineageLoss", Tree.class, NodeRef.class, double.class);
        } catch (Exception e) {
            Assert.fail("Fatal reflection error retrieving method: " + e.toString());
            return;
        }
        m.setAccessible(true);
        double likelihood;
        try {
            likelihood = (Double) m.invoke(model, tree, tree.getRoot(), 1.0);
        } catch (Exception e) {
            Assert.fail("Fatal reflection error invoking method: " + ExceptionUtils.getStackTrace(e));
            return;
        }
	    
	    Assert.assertEquals(actualLikelihood, likelihood, MachineAccuracy.EPSILON);
	}
	
}
