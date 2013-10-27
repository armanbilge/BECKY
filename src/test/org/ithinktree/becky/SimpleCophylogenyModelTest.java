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
import org.junit.Ignore;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import test.org.ithinktree.becky.TestUtils.SimpleBranchRates;
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
	private Parameter hostSwitchRate;
	private Parameter lossRate;
	private BranchRates branchRates;
	private SimpleCophylogenyModel model;
	private Tree host;
	
	@Before
	public void before() {
		
		duplicationRate = new Parameter.Default(1.5);
		hostSwitchRate = new Parameter.Default(0.5);
		lossRate = new Parameter.Default(1.0);
		model = new SimpleCophylogenyModel(duplicationRate, hostSwitchRate, lossRate, Units.Type.YEARS);
		
		// Force internal variables to update
		try {
			model.calculateNodeLogLikelihood(null, null, null, null, null, null, null, null, null);
		} catch (Exception e) {
			// Expected; do nothing
		}
		
		branchRates = new SimpleBranchRates(1.0);
		
		host = TestUtils.DEFAULT_TREE;
	}
		
	@SuppressWarnings("unused")
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
	@Ignore
	public void testLikelihoodHostSwitchEventAndLossInTimeIntegration() {
				
		Method m;
		try {
			m = model.getClass().getDeclaredMethod("likelihoodHostSwitchEventAndLossInTime", double.class, double.class, double.class, double.class, double.class, double.class);
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
	@Ignore
	public void testLikelihoodLineageLoss() {
	    
//	    final Tree tree = TestUtils.treeFromNewick("((A:1.0,B:1.0):2.0,C:3.0);", true);
	    final Tree tree = TestUtils.DEFAULT_TREE;
	    final double overallRate = model.getOverallRate();
	    final double lossRate = model.getLossRate();
//	    final double actualLikelihood = 0.0 + Math.exp(0.0 * overallRate) * (1 - Math.exp(3 * lossRate)) * ((1 - Math.exp(2 * lossRate) + Math.exp(2 * overallRate) * (1 - Math.exp(1 * lossRate) * (1 - Math.exp(1 * lossRate)))));
	    final double actualLikelihood = 1 - Math.exp(-0.5 * lossRate) + Math.exp(-0.5 * overallRate) * (1 - Math.exp(-1.0 * lossRate)) * (1 - Math.exp(-1.0 * lossRate));
	    
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
            likelihood = (Double) m.invoke(model, tree, TestUtils.DE, 1.0);
        } catch (Exception e) {
            Assert.fail("Fatal reflection error invoking method: " + ExceptionUtils.getStackTrace(e));
            return;
        }
	    
	    Assert.assertEquals(actualLikelihood, likelihood, MachineAccuracy.EPSILON);
	}
	
	@Test
	public void testCorrectness() {
	    
	    final boolean testUltrametricity = false;
	    Tree host = TestUtils.treeFromNexus("#NEXUS\n\nBegin taxa;\n    Dimensions ntax=8;\n    Taxlabels\n     host5\n     host8\n     host2\n     host3\n     host7\n     host4\n     host1\n     host6\n     ;\nEnd;\n\nBegin trees;\n   Translate\n     1 host5,\n      2 host8,\n      3 host2,\n      4 host3,\n      5 host7,\n      6 host4,\n      7 host1,\n      8 host6\n       ;\ntree TREE1 = [&R] ((1[&nodeRef=0]:0.2569508196399131,2[&nodeRef=1]:0.2569508196399131)[&nodeRef=8]:0.5523743397767453,(3[&nodeRef=2]:0.38918227984662984,((4[&nodeRef=3]:0.07415946956215326,5[&nodeRef=4]:0.07415946956215326)[&nodeRef=9]:0.16377032613746662,(6[&nodeRef=5]:0.16636101720874424,(7[&nodeRef=6]:0.07184717193772008,8[&nodeRef=7]:0.07184717193772008)[&nodeRef=10]:0.09451384527102416)[&nodeRef=11]:0.07156877849087562)[&nodeRef=12]:0.15125248414700998)[&nodeRef=13]:0.42014287957002866)[&nodeRef=14];\nEnd;\n", testUltrametricity);
	    Tree symbiont = TestUtils.treeFromNexus("#NEXUS\n\nBegin taxa;\n   Dimensions ntax=8;\n    Taxlabels\n     'symbiont5.1'\n     'symbiont8.1'\n     'symbiont3.1'\n     'symbiont5.2'\n     'symbiont3.2'\n     'symbiont3.3'\n     'symbiont7.1'\n     'symbiont8.2'\n     ;\nEnd;\n\nBegin trees;\n   Translate\n     1 'symbiont5.1',\n      2 'symbiont8.1',\n      3 'symbiont3.1',\n      4 'symbiont5.2',\n      5 'symbiont3.2',\n      6 'symbiont3.3',\n      7 'symbiont7.1',\n      8 'symbiont8.2'\n       ;\ntree TREE1 = [&R] ((1[&cophylogeny.clock.rate=1.0,host.nodeRef=0]:0.2569508196399131,2[&cophylogeny.clock.rate=1.0,host.nodeRef=1]:0.2569508196399131)[&cophylogeny.clock.rate=1.0,host.nodeRef=8]:0.5257327028107965,((3[&cophylogeny.clock.rate=1.0,host.nodeRef=3]:0.07381854703766505,4[&cophylogeny.clock.rate=1.0,host.nodeRef=0]:0.07381854703766505)[&cophylogeny.clock.rate=1.0,host.nodeRef=0]:0.18313227260224807,(((5[&cophylogeny.clock.rate=1.0,host.nodeRef=3]:0.014063161669837392,6[&cophylogeny.clock.rate=1.0,host.nodeRef=3]:0.014063161669837392)[&cophylogeny.clock.rate=1.0,host.nodeRef=3]:0.06009630789231587,7[&cophylogeny.clock.rate=1.0,host.nodeRef=4]:0.07415946956215326)[&cophylogeny.clock.rate=1.0,host.nodeRef=9]:0.09861793517860144,8[&cophylogeny.clock.rate=1.0,host.nodeRef=1]:0.1727774047407547)[&cophylogeny.clock.rate=1.0,host.nodeRef=1]:0.08417341489915842)[&cophylogeny.clock.rate=1.0,host.nodeRef=8]:0.5257327028107965)[&cophylogeny.clock.rate=1.0,host.nodeRef=8];\nEnd;\n", testUltrametricity);
	    TestUtils.reduceDecimalPrecision((MutableTree) host); TestUtils.reduceDecimalPrecision((MutableTree) symbiont);
	    CophylogenyLikelihood likelihood = new CophylogenyLikelihood(host, (MutableTree) symbiont, new SimpleCophylogenyModel(new Parameter.Default(1), new Parameter.Default(1), new Parameter.Default(1), Units.Type.YEARS), new SimpleBranchRates(1.0), "host.nodeRef", "ActualCophylogenyLikelihood");
	    for (int i = 0; i < symbiont.getNodeCount(); ++i) {
	        final NodeRef n = symbiont.getNode(i);
	        likelihood.setStatesForNode(n, TestUtils.retrieveNodeByAttribute(host, "nodeRef", symbiont.getNodeAttribute(n, "host.nodeRef")));
	    }
	    likelihood.makeDirty();
	    final double actual = likelihood.getLogLikelihood();
	    System.err.println("Actual: " + actual);
	    host = TestUtils.treeFromNexus("#NEXUS\n\nBegin taxa;\n    Dimensions ntax=8;\n    Taxlabels\n     host1\n     host2\n     host3\n     host4\n     host5\n     host6\n     host7\n     host8\n     ;\nEnd;\n\nBegin trees;\n   Translate\n     1 host1,\n      2 host2,\n      3 host3,\n      4 host4,\n      5 host5,\n      6 host6,\n      7 host7,\n      8 host8\n       ;\ntree TREE1 = [&R] ((2[&height_95%_HPD={0.0,1.1102230246251565E-16},length_range={0.23046200386324736,0.7370906839020797},length_95%_HPD={0.3245025812394004,0.5659201918649017},nodeRef.prob=1.0,height_range={0.0,2.220446049250313E-16},height_median=0.0,rate=0.10000000000001512,height=3.2353238457857857E-17,length=0.44680239577727293,length_median=0.44853395709566735,nodeRef=5,nodeRef.set.prob={1.0},nodeRef.set={\"5\"}]:0.45223823177702316,((7[&height_95%_HPD={0.0,1.1102230246251565E-16},length_range={0.028677072527798487,0.27086623021191286},length_95%_HPD={0.04842216083255659,0.15747058775989778},nodeRef.prob=1.0,height_range={0.0,2.220446049250313E-16},height_median=0.0,rate=0.10000000000001512,height=3.2760275007270476E-17,length=0.0976330363664223,length_median=0.09297308739682483,nodeRef=6,nodeRef.set.prob={1.0},nodeRef.set={\"6\"}]:0.08292479074934718,3[&height_95%_HPD={0.0,1.1102230246251565E-16},length_range={0.028677072527798487,0.27086623021191286},length_95%_HPD={0.04842216083255659,0.15747058775989778},nodeRef.prob=1.0,height_range={0.0,2.220446049250313E-16},height_median=0.0,rate=0.10000000000001512,height=3.2760275007270476E-17,length=0.0976330363664223,length_median=0.09297308739682483,nodeRef=1,nodeRef.set.prob={1.0},nodeRef.set={\"1\"}]:0.08292479074934718)[&height_95%_HPD={0.04842216083255668,0.15747058775989786},length_range={0.028644503645089414,0.36511716473833045},length_95%_HPD={0.07434839816421698,0.24693063777263619},nodeRef.prob=1.0,height_range={0.028677072527798497,0.27086623021191286},height_median=0.09297308739682497,rate=0.10000000000001512,height=0.09763303636642233,posterior=1.0,length=0.16252964739790784,length_median=0.15957113027558315,nodeRef=17,nodeRef.set.prob={1.0},nodeRef.set={\"17\"}]:0.15134506665038544,((6[&height_95%_HPD={0.0,1.1102230246251565E-16},length_range={0.01584572219570872,0.21252168107847447},length_95%_HPD={0.030892854495931436,0.12078890755013967},nodeRef.prob=1.0,height_range={0.0,2.220446049250313E-16},height_median=0.0,rate=0.10000000000001512,height=3.255058951211852E-17,length=0.07512316652123423,length_median=0.07282731101171622,nodeRef=11,nodeRef.set.prob={1.0},nodeRef.set={\"11\"}]:0.07685449913525907,1[&height_95%_HPD={0.0,1.1102230246251565E-16},length_range={0.01584572219570872,0.21252168107847447},length_95%_HPD={0.030892854495931436,0.12078890755013967},nodeRef.prob=1.0,height_range={0.0,2.220446049250313E-16},height_median=0.0,rate=0.10000000000001512,height=3.255058951211852E-17,length=0.07512316652123423,length_median=0.07282731101171622,nodeRef=2,nodeRef.set.prob={1.0},nodeRef.set={\"2\"}]:0.07685449913525907)[&height_95%_HPD={0.030892854495931488,0.12078890755013971},length_range={0.005490802247484633,0.21613965892135223},length_95%_HPD={0.022901791812288685,0.1366821708056959},nodeRef.prob=1.0,height_range={0.015845722195708678,0.21252168107847447},height_median=0.07282731101171624,rate=0.10000000000001512,height=0.07512316652123427,posterior=1.0,length=0.0751611162906951,length_median=0.07194884709442785,nodeRef=15,nodeRef.set.prob={1.0},nodeRef.set={\"15\"}]:0.038813331591478484,4[&height_95%_HPD={0.0,1.1102230246251565E-16},length_range={0.052647993344648175,0.3194175996720713},length_95%_HPD={0.08665222181443676,0.21647701739422875},nodeRef.prob=1.0,height_range={0.0,2.220446049250313E-16},height_median=0.0,rate=0.10000000000001512,height=3.204487743557557E-17,length=0.15028673974749712,length_median=0.14794532896673684,nodeRef=0,nodeRef.set.prob={1.0},nodeRef.set={\"0\"}]:0.11566783072673756)[&height_95%_HPD={0.08665222181443677,0.21647701739422875},length_range={0.011105836059042196,0.3080057253439866},length_95%_HPD={0.03774966278502062,0.1880978395081588},nodeRef.prob=1.0,height_range={0.0526479933446482,0.3194175996720713},height_median=0.14794465093602532,rate=0.10000000000001512,height=0.15028042870315572,posterior=0.9998889012331963,length=0.10989060966361827,length_median=0.10596793734926917,nodeRef=16,nodeRef.set.prob={1.0},nodeRef.set={\"16\"}]:0.11860202667299505)[&height_95%_HPD={0.17864305554085924,0.34467862361301693},length_range={0.016395317535610215,0.4375434641839056},length_95%_HPD={0.075186577219643,0.2990028702391133},nodeRef.prob=1.0,height_range={0.1310076386668093,0.4588933523182276},height_median=0.25744334778398803,rate=0.10000000000001512,height=0.26016514069989816,posterior=1.0,length=0.18663725507737486,length_median=0.1853851227394095,nodeRef=18,nodeRef.set.prob={1.0},nodeRef.set={\"18\"}]:0.21796837437729055)[&height_95%_HPD={0.3245025812394005,0.5659201918649017},length_range={0.006877334608307306,0.7419751942087958},length_95%_HPD={0.15577271580854596,0.5179075611211502},nodeRef.prob=1.0,height_range={0.23046200386324733,0.7370906839020797},height_median=0.44853395709566735,rate=0.10000000000001512,height=0.446802395777273,posterior=1.0,length=0.3339984802826845,length_median=0.3290954598720502,nodeRef=19,nodeRef.set.prob={1.0},nodeRef.set={\"19\"}]:0.28795627541578195,(5[&height_95%_HPD={0.0,1.1102230246251565E-16},length_range={0.14473408501998142,0.6189374924242743},length_95%_HPD={0.18599701879387529,0.46853392361190904},nodeRef.prob=1.0,height_range={0.0,2.220446049250313E-16},height_median=0.0,rate=0.10000000000001512,height=3.2180556285379775E-17,length=0.3116409902175977,length_median=0.28981788686075277,nodeRef=10,nodeRef.set.prob={1.0},nodeRef.set={\"10\"}]:0.25857976483709205,8[&height_95%_HPD={0.0,1.1102230246251565E-16},length_range={0.14473408501998142,0.6189374924242743},length_95%_HPD={0.18599701879387529,0.46853392361190904},nodeRef.prob=1.0,height_range={0.0,2.220446049250313E-16},height_median=0.0,rate=0.10000000000001512,height=3.2180556285379775E-17,length=0.3116409902175977,length_median=0.28981788686075277,nodeRef=8,nodeRef.set.prob={1.0},nodeRef.set={\"8\"}]:0.25857976483709205)[&height_95%_HPD={0.1859970187938753,0.46853392361190904},length_range={0.1533789856598496,0.8420303140085965},length_95%_HPD={0.28167990818653676,0.6434244145803146},nodeRef.prob=1.0,height_range={0.14473408501998142,0.6189374924242744},height_median=0.2898178868607529,rate=0.10000000000001512,height=0.3116409902175977,posterior=1.0,length=0.4691598858423584,length_median=0.469975957441123,nodeRef=20,nodeRef.set.prob={1.0},nodeRef.set={\"20\"}]:0.48161474235571305)[&height_95%_HPD={0.6269360482675888,0.9412304082728955},height_median=0.778128753112989,height=0.7808008760599581,nodeRef.prob=1.0,posterior=1.0,height_range={0.5333184432079843,1.1498137263195063},length=0.0,nodeRef=14,nodeRef.set.prob={1.0},nodeRef.set={\"14\"}];\nEnd;\n", testUltrametricity);
	    symbiont = TestUtils.treeFromNexus("#NEXUS\n\nBegin taxa;\n    Dimensions ntax=8;\n    Taxlabels\n     'symbiont3.1'\n     'symbiont3.2'\n     'symbiont3.3'\n     'symbiont5.1'\n     'symbiont5.2'\n     'symbiont7.1'\n     'symbiont8.1'\n     'symbiont8.2'\n     ;\nEnd;\n\nBegin trees;\n   Translate\n     1 'symbiont3.1',\n      2 'symbiont3.2',\n      3 'symbiont3.3',\n      4 'symbiont5.1',\n      5 'symbiont5.2',\n      6 'symbiont7.1',\n      7 'symbiont8.1',\n      8 'symbiont8.2'\n       ;\ntree TREE1 = [&R] ((4[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.0,5.551115123125783E-17},length_range={0.07009040640376611,0.3716503013795565},length_95%_HPD={0.10655358812509273,0.25434129077706763},height_range={0.0,1.6653345369377348E-16},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=1.0,height_median=0.0,rate=0.10000000000001512,host.nodeRef=10,height=1.8082290346633478E-17,host.nodeRef.set.prob={1.0},length=0.18215182466676827,length_median=0.18127256705675637,host.nodeRef.set={\"10\"}]:0.1800683997312161,7[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.0,5.551115123125783E-17},length_range={0.07009040640376611,0.3716503013795565},length_95%_HPD={0.10655358812509273,0.25434129077706763},height_range={0.0,1.6653345369377348E-16},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=1.0,height_median=0.0,rate=0.10000000000001512,host.nodeRef=8,height=1.8082290346633478E-17,host.nodeRef.set.prob={1.0},length=0.18215182466676827,length_median=0.18127256705675637,host.nodeRef.set={\"8\"}]:0.1800683997312161)[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.10655358812509275,0.25434129077706763},length_range={0.05147873457595642,0.5537726979267441},length_95%_HPD={0.13812986310408165,0.33165616563399625},height_range={0.07009040640376618,0.3716503013795565},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=0.5383846239306743,height_median=0.18127256705675637,rate=0.10000000000001512,host.nodeRef=5,height=0.18215182466676832,posterior=1.0,host.nodeRef.set.prob={0.0012220864348405733,0.011665370514387291,0.0052216420397733585,1.1109876680368848E-4,0.0056660371069881124,0.5383846239306743,0.21897566937007,0.2180868792356405,6.665926008221309E-4},length=0.23622984008105674,length_median=0.23436865262793105,host.nodeRef.set={\"0\",\"17\",\"16\",\"19\",\"18\",\"5\",\"8\",\"10\",\"15\"}]:0.25026973308457523,((8[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.0,5.551115123125783E-17},length_range={0.05701531361339068,0.2851622867291413},length_95%_HPD={0.0894246842479037,0.20060687336287852},height_range={0.0,1.1102230246251565E-16},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=1.0,height_median=0.0,rate=0.10000000000001512,host.nodeRef=8,height=1.8076123126187834E-17,host.nodeRef.set.prob={1.0},length=0.1449442620329724,length_median=0.1435129534380185,host.nodeRef.set={\"8\"}]:0.11978843166872244,(6[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.0,5.551115123125783E-17},length_range={0.0073304247779785205,0.14394435311051096},length_95%_HPD={0.021224625346615025,0.08422317207986615},height_range={0.0,1.1102230246251565E-16},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=1.0,height_median=0.0,rate=0.10000000000001512,host.nodeRef=6,height=1.839065136891577E-17,host.nodeRef.set.prob={1.0},length=0.052247586945678624,length_median=0.051035501128845875,host.nodeRef.set={\"6\"}]:0.039780088811651106,(2[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.0,5.551115123125783E-17},length_range={1.343547058735914E-5,0.08961128010649873},length_95%_HPD={0.001762859609897134,0.05102525760617182},height_range={0.0,2.220446049250313E-16},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=1.0,height_median=0.0,rate=0.10000000000001512,host.nodeRef=1,height=1.8501661336937394E-17,host.nodeRef.set.prob={1.0},length=0.025116707577620377,length_median=0.02347737437976137,host.nodeRef.set={\"1\"}]:0.02648658874232559,3[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.0,5.551115123125783E-17},length_range={1.343547058735914E-5,0.08961128010649873},length_95%_HPD={0.0022393212204063905,0.0514966860316535},height_range={0.0,2.220446049250313E-16},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=1.0,height_median=0.0,rate=0.10000000000001512,host.nodeRef=1,height=1.8520162998274332E-17,host.nodeRef.set.prob={1.0},length=0.025104737379422936,length_median=0.02345882280108241,host.nodeRef.set={\"1\"}]:0.02648658874232559)[&cophylogeny.rate_median=0.0514409303850007,cophylogeny.rate=0.061328371187936416,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.0013831430166226788,0.050392295767577644},length_range={2.673242509979301E-4,0.11798811330756709},length_95%_HPD={0.0020297658993243227,0.0565488395883843},height_range={1.343547058740846E-5,0.08961128010649877},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14233663375031858},host.nodeRef.prob=0.25,height_median=0.02340983019625495,rate=0.10000000000001505,host.nodeRef=5,height=0.024987252388872942,posterior=0.9914453949561159,host.nodeRef.set.prob={0.14769161810847153,0.051434334379202154,0.07496638278798745,0.25,0.05569251456745854,0.13637382339757956,0.20887494397131331,0.07451815329448677,4.4822949350067237E-4},length=0.027412880901116037,length_median=0.025030063277988214,host.nodeRef.set={\"0\",\"1\",\"2\",\"5\",\"6\",\"8\",\"10\",\"11\",\"15\"}]:0.013293500069325517)[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.021024607624037572,0.08378630591072167},length_range={0.020731157240184372,0.21207463845952368},length_95%_HPD={0.04111396816739919,0.1484356920565057},height_range={0.007330424777978517,0.14394435311051096},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=0.295300522164204,height_median=0.05107408639623423,rate=0.10000000000001512,host.nodeRef=5,height=0.05233060707514856,posterior=1.0,host.nodeRef.set.prob={5.554938340184424E-4,0.12576380402177537,1.1109876680368848E-4,0.03499611154316187,0.019331185423841796,0.295300522164204,0.14509498944561716,0.11409843350738806,0.2438617931340962,0.01910898789023442,0.0017775802688590157},length=0.09258483352596092,length_median=0.09030021234818227,host.nodeRef.set={\"17\",\"0\",\"16\",\"1\",\"2\",\"5\",\"6\",\"8\",\"10\",\"11\",\"15\"}]:0.08000834285707134)[&cophylogeny.rate_median=0.05143882050455712,cophylogeny.rate=0.061246059292765546,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.08949330823199203,0.20045761842822718},length_range={0.0010031845197779887,0.18761854416671098},length_95%_HPD={0.007760628019016308,0.09578234703676698},height_range={0.057015313613390683,0.27564691723574447},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14225798069028062},host.nodeRef.prob=0.3988368191477463,height_median=0.14333887536786488,rate=0.10000000000001508,host.nodeRef=5,height=0.14473197429872023,posterior=0.9933340739917786,host.nodeRef.set.prob={0.04194161726876188,0.055586623420199086,0.004138239570517839,0.014092383402303992,4.473772508667934E-4,4.473772508667934E-4,0.3988368191477463,0.005032994072251426,0.285874063303881,0.18308913991723522,3.355329381500951E-4,0.01017783245721955},length=0.04817344654963988,length_median=0.04501582910494786,host.nodeRef.set={\"0\",\"17\",\"1\",\"16\",\"2\",\"18\",\"5\",\"6\",\"8\",\"10\",\"11\",\"15\"}]:0.0714425824056546,(1[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.0,5.551115123125783E-17},length_range={0.006305254893198639,0.1553358488759628},length_95%_HPD={0.01920283518465959,0.1009720144404646},height_range={0.0,1.1102230246251565E-16},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=1.0,height_median=0.0,rate=0.10000000000001512,host.nodeRef=1,height=1.8088457567079125E-17,host.nodeRef.set.prob={1.0},length=0.059866139884412625,length_median=0.05785937105881294,host.nodeRef.set={\"1\"}]:0.050189624987262184,5[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.0,5.551115123125783E-17},length_range={0.006305254893198639,0.1553358488759628},length_95%_HPD={0.01920283518465959,0.1009720144404646},height_range={0.0,1.1102230246251565E-16},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=1.0,height_median=0.0,rate=0.10000000000001512,host.nodeRef=10,height=1.8088457567079125E-17,host.nodeRef.set.prob={1.0},length=0.059866139884412625,length_median=0.05785937105881294,host.nodeRef.set={\"10\"}]:0.050189624987262184)[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.01920283518465965,0.10097201444046466},length_range={0.027478699745117394,0.2728910043506471},length_95%_HPD={0.06628590825741802,0.199723117193705},height_range={0.006305254893198664,0.1553358488759628},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=0.6735918231307633,height_median=0.057859371058813,rate=0.10000000000001512,host.nodeRef=10,height=0.05986613988441263,posterior=1.0,host.nodeRef.set.prob={0.003221864237306966,0.023108543495167205,0.06477058104655038,0.0019997778024663927,0.1330963226308188,0.00655482724141762,0.0896567048105766,0.6735918231307633,0.002777469170092212,0.0012220864348405733},length=0.13286089112434857,length_median=0.13122573351093186,host.nodeRef.set={\"17\",\"0\",\"1\",\"2\",\"5\",\"6\",\"8\",\"10\",\"11\",\"15\"}]:0.14104138908711486)[&cophylogeny.rate_median=0.05141583602754825,cophylogeny.rate=0.06124855753938849,cophylogeny.rate_range={1.2099144363748832E-5,0.6182447878878883},height_95%_HPD={0.13124223177467514,0.25416207774125},length_range={0.06627952819748323,0.6193851450870963},length_95%_HPD={0.1296717611755402,0.32712516677401005},height_range={0.07753861110433447,0.33016684660160206},cophylogeny.rate_95%_HPD={0.0071942748227629425,0.14200933981680006},host.nodeRef.prob=0.6701477613598489,height_median=0.191792820142629,rate=0.10000000000001512,host.nodeRef=5,height=0.19283732346407928,posterior=1.0,host.nodeRef.set.prob={0.0013331852016442618,0.014998333518497945,1.1109876680368848E-4,0.004221753138540162,0.004666148205754916,0.002110876569270081,0.6701477613598489,1.1109876680368848E-4,0.13542939673369625,0.16620375513831798,6.665926008221309E-4},length=0.2255443412837476,length_median=0.22423876650718924,host.nodeRef.set={\"0\",\"17\",\"1\",\"16\",\"18\",\"20\",\"5\",\"6\",\"8\",\"10\",\"15\"}]:0.23910711874141433)[&host.nodeRef.prob=0.3519608932340851,height_95%_HPD={0.3200994366901407,0.5162858841746638},height_median=0.41702726893519726,height=0.4183816647478268,host.nodeRef=18,host.nodeRef.set.prob={0.0026663704032885236,0.0036662593045217197,0.013665148316853682,0.3519608932340851,0.28485723808465724,0.13654038440173313,0.1060993222975225,0.10054438395733807},posterior=1.0,height_range={0.26235266320059347,0.7738115184924209},length=0.0,host.nodeRef.set={\"17\",\"16\",\"19\",\"18\",\"20\",\"5\",\"8\",\"10\"}];\nEnd;\n", testUltrametricity);
	    TestUtils.reduceDecimalPrecision((MutableTree) host); TestUtils.reduceDecimalPrecision((MutableTree) symbiont);
	    likelihood = new CophylogenyLikelihood(host, (MutableTree) symbiont, new SimpleCophylogenyModel(new Parameter.Default(1.8184E-9), new Parameter.Default(111.1071), new Parameter.Default(1.0), Units.Type.YEARS), new SimpleBranchRates(5.1416E-2), "host.nodeRef", "EstimatedCophylogenyLikelihood");
	    
        for (int i = 0; i < symbiont.getNodeCount(); ++i) {
            final NodeRef n = symbiont.getNode(i);
            likelihood.setStatesForNode(n, TestUtils.retrieveNodeByAttribute(host, "nodeRef", symbiont.getNodeAttribute(n, "host.nodeRef")));
        }
	    
	    final double estimated = likelihood.getLogLikelihood();
	    System.err.println("Estimated: " + estimated);
	    Assert.assertTrue(estimated > actual);
	    
	}
		
}
