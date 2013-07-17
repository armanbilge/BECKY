/**
 * CophylogenyLikelihood.java
 * 
 * BECKY - Bayesian Estimation of Coevolutionary KrYteria
 * 
 */

package org.ithinktree.becky;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTrait.DefaultBehavior;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
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

	public static final String STATES_KEY = "host.nodeRef";
	public static final int NO_HOST = -1;
	
	private Tree hostTree;
	private Tree symbiontTree;
	private CophylogenyModel cophylogenyModel;
	private BranchRateModel branchRateModel;
	private String hostAttributeName;
	private String reconstructionTagName;
	
	private TreeTraitProvider.Helper treeTraits = new Helper();
	
	/**
	 * 
	 */
	
	public CophylogenyLikelihood(Tree hostTree, Tree symbiontTree, CophylogenyModel cophylogenyModel, BranchRateModel branchRateModel, String reconstructionTagName, String hostAttributeName, String id) {
		this(CophylogenyLikelihoodParser.COPHYLOGENY_LIKELIHOOD, hostTree, symbiontTree, cophylogenyModel, branchRateModel, reconstructionTagName, hostAttributeName);
		setId(id);
	}
	
	public CophylogenyLikelihood(String name, Tree hostTree, Tree symbiontTree, CophylogenyModel cophylogenyModel, BranchRateModel branchRateModel, final String reconstructionTagName, String hostAttributeName) {
		
		super(name);
		
		this.hostTree = hostTree;
		this.symbiontTree = symbiontTree;
		this.cophylogenyModel = cophylogenyModel;
		this.branchRateModel = branchRateModel;
		if (hostTree instanceof TreeModel) {
			addModel((TreeModel) hostTree);
		}
		if (symbiontTree instanceof TreeModel) {
			addModel((TreeModel) symbiontTree);
		}
		if (cophylogenyModel != null) {
			addModel(cophylogenyModel);
		}
		if (branchRateModel != null) {
			addModel(branchRateModel);
		}
		this.hostAttributeName = hostAttributeName;
		this.reconstructionTagName = reconstructionTagName;
		
		reconstructedStates = new int[symbiontTree.getNodeCount()];
		storedReconstructedStates = new int[symbiontTree.getNodeCount()];
		
		treeTraits.addTrait(STATES_KEY, new NR() {

			@Override
			public String getTraitName() {
				return reconstructionTagName;
			}

			@Override
			public Intent getIntent() {
				return Intent.NODE;
			}

			@Override
			public NodeRef getTrait(Tree tree, NodeRef node) {
				return getStatesForNode(node);
			}
			
			public String getTraitString(Tree tree, NodeRef node) {
				return formatTrait(getStatesForNode(node));
			}
			
		});
		
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

	private double calculateLogLikelihood() {
						
		double logL = 0.0;
		
		NodeRef self, child1, child2, parent;
		NodeRef selfHost, child1Host, child2Host;
		double branchRate, selfNodeHeight, selfBranchTime, child1BranchTime, child2BranchTime;
		for (int i = 0; i < symbiontTree.getInternalNodeCount(); i++) {

			self = symbiontTree.getInternalNode(i);
						
			child1 = symbiontTree.getChild(self, 0);
			child2 = symbiontTree.getChild(self, 1);
			
			selfHost = getStatesForNode(self);
			child1Host = getStatesForNode(child1);
			child2Host = getStatesForNode(child2);
			
			selfNodeHeight = symbiontTree.getNodeHeight(self);
			if (symbiontTree.isRoot(self)) {
				selfBranchTime = -1.0;
			} else {
				parent = symbiontTree.getParent(self);
				branchRate = branchRateModel.getBranchRate(symbiontTree, self);
				selfBranchTime = branchRate * (symbiontTree.getNodeHeight(parent) - selfNodeHeight);
			}
			
			branchRate = branchRateModel.getBranchRate(symbiontTree, child1);
			child1BranchTime = branchRate * (selfNodeHeight - symbiontTree.getNodeHeight(child1));
				
			branchRate = branchRateModel.getBranchRate(symbiontTree, child2);
			child2BranchTime = branchRate * (selfNodeHeight - symbiontTree.getNodeHeight(child2));
			
			logL += cophylogenyModel.calculateNodeLogLikelihood(self, hostTree, selfHost, child1Host, child2Host, selfBranchTime, child1BranchTime, child2BranchTime);
				
			if (logL == Double.NEGATIVE_INFINITY) break;
			
		}
			
		return logL;
	}

	@Override
	public void makeDirty() {
		likelihoodKnown = false;
		
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		likelihoodKnown = false;
		
	}

	@SuppressWarnings("rawtypes")
	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {		
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
	protected void acceptState() {		
	} // Nothing to do

	
	private int[] reconstructedStates;
	private int[] storedReconstructedStates;
	
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
		if (state == null) {
			reconstructedStates[node.getNumber()] = NO_HOST;
		} else {
			reconstructedStates[node.getNumber()] = state.getNumber();
		}
	}
	
	public abstract class NR extends DefaultBehavior implements TreeTrait<NodeRef> {
		
		@SuppressWarnings("rawtypes")
		public Class getTraitClass() {
			return NodeRef.class;
		}
		
		public String getTraitString(Tree tree, NodeRef node) {
			return formatTrait(getTrait(tree, node));
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
	
}
