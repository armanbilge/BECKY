package org.ithinktree.becky.JrRIMEWrappers;

import java.util.HashMap;
import java.util.Iterator;

import org.ithinktree.becky.CophylogenyLikelihood;

import se.cbb.jprime.apps.dltrs.DLTRSModel;
import se.cbb.jprime.apps.dltrs.EpochDLTProbs;
import se.cbb.jprime.apps.dltrs.ReconciliationHelper;
import se.cbb.jprime.math.Continuous1DPDDependent;
import se.cbb.jprime.mcmc.ChangeInfo;
import se.cbb.jprime.mcmc.Dependent;
import se.cbb.jprime.topology.GuestHostMap;
import se.cbb.jprime.topology.LeafLeafMap;
import se.cbb.jprime.topology.NamesMap;
import se.cbb.jprime.topology.RBTree;
import se.cbb.jprime.topology.RBTreeEpochDiscretiser;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;

@SuppressWarnings("serial")
public class CophylogenyLikelihoodWrapperForJPRIMEDLTRSModel extends
		CophylogenyLikelihood {

	private final DLTRSModel dltrsModel;
	private final ReconciliationHelper reconciliationHelper;

	public CophylogenyLikelihoodWrapperForJPRIMEDLTRSModel(String name, Tree guest, Tree host,
			Parameter duplicationRate, Parameter lossRate, Parameter transferRate, boolean normalize, Continuous1DPDDependent substPD, Taxa guestTaxa, String hostAttributeName) {
		
		super(name);
		
		final RBTree jprimeHostRBTree = new JPrIMERBTreeWrapperForBEASTTree(host);
		final NamesMap jprimeHostNamesMap = new JPrIMENamesMapWrapperForBEASTTree(host);
		final RBTree jprimeGuestRBTree = new JPrIMERBTreeWrapperForBEASTTree(guest);
		final NamesMap jprimeGuestNamesMap = new JPrIMENamesMapWrapperForBEASTTree(guest);
		final RBTreeEpochDiscretiser rbTreeEpochDiscretiser = new RBTreeEpochDiscretiser(jprimeHostRBTree, jprimeHostNamesMap, new JPrIMETimesMapWrapperForBEASTTree(host));
		
		final GuestHostMap guestHostMap = new GuestHostMap();
		for (Iterator<Taxon> taxa = guestTaxa.iterator(); taxa.hasNext(); ) {
			Taxon t = taxa.next();
			guestHostMap.add(t.getId(), ((Taxon) t.getAttribute(hostAttributeName)).getId());
		}
		
		reconciliationHelper = new ReconciliationHelper(jprimeGuestRBTree, jprimeHostRBTree, rbTreeEpochDiscretiser, new LeafLeafMap(guestHostMap, jprimeGuestRBTree, jprimeGuestNamesMap, jprimeHostRBTree, jprimeHostNamesMap));
		
		EpochDLTProbs dltProbs = new EpochDLTProbs(rbTreeEpochDiscretiser,
				new JPrIMEDoubleParameterWrapperForBEASTParameter(duplicationRate),
				new JPrIMEDoubleParameterWrapperForBEASTParameter(lossRate),
				new JPrIMEDoubleParameterWrapperForBEASTParameter(transferRate),
				normalize
				);
				
		dltrsModel = new DLTRSModel(jprimeGuestRBTree, jprimeHostRBTree, reconciliationHelper, new JPrIMEDoubleMapWrapperForBEASTTree(guest), dltProbs, substPD);

	}

	public double getLogLikelihood() {
		dltrsModel.cacheAndUpdate(new HashMap<Dependent,ChangeInfo>(), true);
		return dltrsModel.getDataProbability().getLogValue();
	}
	
	@Override
	protected void storeState() {
		dltrsModel.cacheAndUpdate(new HashMap<Dependent,ChangeInfo>(), true);
	}

	protected void restoreState() {
		dltrsModel.restoreCache(true);
	}

	
}
