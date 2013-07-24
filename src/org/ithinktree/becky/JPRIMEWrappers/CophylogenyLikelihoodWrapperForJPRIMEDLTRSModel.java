package org.ithinktree.becky.JPRIMEWrappers;

import java.util.HashMap;

import org.ithinktree.becky.CophylogenyLikelihood;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import se.cbb.jprime.apps.dltrs.DLTRSModel;
import se.cbb.jprime.apps.dltrs.EpochDLTProbs;
import se.cbb.jprime.apps.dltrs.ReconciliationHelper;
import se.cbb.jprime.math.Continuous1DPDDependent;
import se.cbb.jprime.mcmc.ChangeInfo;
import se.cbb.jprime.mcmc.Dependent;
import se.cbb.jprime.mcmc.DoubleParameter;
import se.cbb.jprime.topology.DoubleMap;
import se.cbb.jprime.topology.LeafLeafMap;
import se.cbb.jprime.topology.RootedBifurcatingTreeParameter;

@SuppressWarnings("serial")
public class CophylogenyLikelihoodWrapperForJPRIMEDLTRSModel extends
		CophylogenyLikelihood {

	private final DLTRSModel dltrsModel;
	private final ReconciliationHelper reconciliationHelper;
	private final DoubleMap lengths;
	private final EpochDLTProbs dltProbs;
	private final DoubleParameter jprimeDuplicationRate;
	private final DoubleParameter jprimeLossRate;
	private final DoubleParameter jprimeTransferRate;
	
	// BEAST stuff
	private final Tree guest;
	private final Tree host;
	private final Parameter duplicationRate;
	private final Parameter lossRate;
	private final Parameter transferRate;
	
	@SuppressWarnings("unchecked")
	public CophylogenyLikelihoodWrapperForJPRIMEDLTRSModel(String name, RootedBifurcatingTreeParameter guest, RootedBifurcatingTreeParameter host,
			EpochDLTProbs dltProbs, Continuous1DPDDependent substPD, LeafLeafMap leafLeafMap) {
		
		super(name);
		
		if (guest instanceof WrappedBEASTObject) {
			Object o = ((WrappedBEASTObject<?>) guest).getWrappedBEASTObject();
			this.guest = (Tree) o;
			if (o instanceof TreeModel)
				addModel((TreeModel) o);
		}
		
		if (host instanceof WrappedBEASTObject) {
			Object o = ((WrappedBEASTObject<?>) host).getWrappedBEASTObject();
			this.host = (Tree) o;
			if (o instanceof TreeModel)
				addModel((TreeModel) o);

		}
		
		reconciliationHelper = new ReconciliationHelper(leafLeafMap);
		
		dltProbs = new EpochDLTProbs();
		
		lengths = new DoubleMap("lengths", this.guest.getNodeCount());
		
		dltrsModel = new DLTRSModel(guest, host, reconciliationHelper, lengths);
	}

	public double getLogLikelihood() {
		updateJPrIMEModels();
		dltrsModel.cacheAndUpdate(new HashMap<Dependent,ChangeInfo>(), true);
		return dltrsModel.getDataProbability().getLogValue();
	}

	private void updateJPrIMEModels() {
		for (int i = 0; i < guest.getNodeCount(); ++i)
			lengths.set(i, guest.getBranchLength(guest.getNode(i)));
		
	}
	
	@Override
	protected void storeState() {}

	protected void restoreState() {
		dltrsModel.restoreCache(true);
	}

	
}
