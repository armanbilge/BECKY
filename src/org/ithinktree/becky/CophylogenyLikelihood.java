/**
 * CophylogenyLikelihood.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */

package org.ithinktree.becky;

import org.ithinktree.becky.xml.CophylogenyLikelihoodParser;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTrait.DefaultBehavior;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Units;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

/**
 * A likelihood function for cophylogenetic processes.
 * <p/>
 * Code form based on the SpeciationLikelihood class.
 * 
 * @author Arman D. Bilge
 *
 */
@SuppressWarnings("serial")
public class CophylogenyLikelihood extends AbstractModelLikelihood implements TreeTraitProvider, Units {

	public static final int NO_HOST = -1;
	
	final private Tree hostTree;
	final private MutableTree symbiontTree;
	final private CophylogenyModel cophylogenyModel;
	final private BranchRates branchRates;
	
	final private TreeTraitProvider.Helper treeTraits = new Helper();
		
	public CophylogenyLikelihood(final Tree hostTree, final MutableTree symbiontTree, final CophylogenyModel cophylogenyModel, final BranchRates branchRateModel, final String reconstructionTagName, final String id) {
		this(CophylogenyLikelihoodParser.COPHYLOGENY_LIKELIHOOD, hostTree, symbiontTree, cophylogenyModel, branchRateModel, reconstructionTagName);
		setId(id);
	}
	
	public CophylogenyLikelihood(final String name, final Tree hostTree, final MutableTree symbiontTree, final CophylogenyModel cophylogenyModel, final BranchRates branchRates, final String reconstructionTagName) {
		
		super(name);
		
		this.hostTree = hostTree;
		this.symbiontTree = symbiontTree;
		
		this.cophylogenyModel = cophylogenyModel;
		this.branchRates = branchRates;
		
		if (symbiontTree instanceof Model) {
			addModel((Model) symbiontTree);
		}
		if (hostTree instanceof Model) {
			addModel((Model) hostTree);
		}
		if (cophylogenyModel != null) {
			addModel(cophylogenyModel);
		}
		if (branchRates instanceof Model) {
			addModel((Model) branchRates);
		}
		
		reconstructedStates = new int[symbiontTree.getNodeCount()];
		storedReconstructedStates = new int[reconstructedStates.length];
		cophylogenyModel.initialize(symbiontTree);
		
		
		treeTraits.addTrait(reconstructionTagName, new NodeRefTrait() {
			@Override
			public String getTraitName() {
				return reconstructionTagName;
			}
			@Override
			public NodeRef getTrait(final Tree tree, final NodeRef node) {
				return getStatesForNode(node);
			}
			public String getTraitString(final Tree tree, final NodeRef node) {
				return formatTrait(getTrait(tree, node));
			}
		});
		
	}
	
	public CophylogenyLikelihood(String name) {
		super(name);
		symbiontTree = null;
		hostTree = null;
		cophylogenyModel = null;
		branchRates = null;
		final int[] empty = new int[0];
		reconstructedStates = empty;
		storedReconstructedStates = empty;
	}

	@Override
	public Model getModel() {
		return this;
	}

	public boolean evaluateEarly() {
		return true;
	}
	
	@Override
	public double getLogLikelihood() {
		if (!likelihoodKnown) {
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}
		return logLikelihood;
	}

	protected double calculateLogLikelihood() {
		
		double logL = 0.0;
		cophylogenyModel.updateVariables();
		NodeRef self, child1, child2;
		NodeRef selfHost, child1Host, child2Host;
		self = symbiontTree.getRoot();
		do {
			self = Tree.Utils.postorderSuccessor(symbiontTree, self);
			if (!symbiontTree.isExternal(self)) {
				
				child1 = symbiontTree.getChild(self, 0);
				child2 = symbiontTree.getChild(self, 1);
				
				selfHost = getStatesForNode(self);
				child1Host = getStatesForNode(child1);
				child2Host = getStatesForNode(child2);
				
 				logL += cophylogenyModel.calculateNodeLogLikelihood(symbiontTree, self, child1, child2, hostTree, selfHost, child1Host, child2Host, branchRates);
			}
		} while (!symbiontTree.isRoot(self) && logL != Double.NEGATIVE_INFINITY);

		return logL;
	}

	@Override
	public void makeDirty() {
		likelihoodKnown = false;
	}
	
	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		makeDirty();
	}

	@SuppressWarnings("rawtypes")
	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		// Nothing to do; no variables
	}

	@Override
	protected void storeState() {
				
		storedLikelihoodKnown = likelihoodKnown;
		storedLogLikelihood = logLikelihood;
		
		System.arraycopy(reconstructedStates, 0, storedReconstructedStates, 0, reconstructedStates.length);
		
	}

	@Override
	protected void restoreState() {
		
		likelihoodKnown = storedLikelihoodKnown;
		logLikelihood = storedLogLikelihood;
		
		System.arraycopy(storedReconstructedStates, 0, reconstructedStates, 0, storedReconstructedStates.length);
		
	}

	@Override
	protected void acceptState() {} // Nothing to do

	
	private final int[] reconstructedStates;
	private final int[] storedReconstructedStates;
	
	private double logLikelihood;
	private double storedLogLikelihood;
	private boolean likelihoodKnown = false;
	private boolean storedLikelihoodKnown = false;

	@Override
	public Type getUnits() {
		return cophylogenyModel.getUnits();
	}

	@Override
	public void setUnits(Type units) {
		cophylogenyModel.setUnits(units);
	}
	
	public NodeRef getStatesForNode(NodeRef node) {
		
		int hostIndex = reconstructedStates[node.getNumber()];
		if (hostIndex == NO_HOST) 
			return null;
		return hostTree.getNode(hostIndex);
		
	}
	
	public void setStatesForNode(NodeRef node, NodeRef state) {
		
		likelihoodKnown = false;
		reconstructedStates[node.getNumber()] = state == null ? NO_HOST : state.getNumber();
		fireModelChanged();

	}
	
	public abstract class NodeRefTrait extends DefaultBehavior implements TreeTrait<NodeRef> {
		
		@SuppressWarnings("rawtypes")
		public Class getTraitClass() {
			return NodeRef.class;
		}
		
		public String getTraitString(Tree tree, NodeRef node) {
			return formatTrait(getTrait(tree, node));
		}
		
		@Override
		public Intent getIntent() {
			return Intent.NODE;
		}
		
		public String formatTrait(NodeRef n) {
			if (n == null) {
				return Integer.toString(-1);
			}
			return Integer.toString(n.getNumber());
		}
		
	}
	
	@SuppressWarnings("rawtypes")
	@Override
	public TreeTrait[] getTreeTraits() {
		return treeTraits.getTreeTraits();
	}

	@SuppressWarnings("rawtypes")
	@Override
	public TreeTrait getTreeTrait(String key) {
		return treeTraits.getTreeTrait(key);
	}

	protected boolean isDirty() { return !likelihoodKnown; }
	
}
