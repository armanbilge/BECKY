package org.ithinktree.becky.jprimewrappers;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.ithinktree.becky.CophylogenyLikelihood;

import se.cbb.jprime.apps.dltrs.DLTRSModel;
import se.cbb.jprime.apps.dltrs.EpochDLTProbs;
import se.cbb.jprime.apps.dltrs.ReconciliationHelper;
import se.cbb.jprime.math.Continuous1DPDDependent;
import se.cbb.jprime.math.PRNG;
import se.cbb.jprime.math.RealInterval;
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
			Parameter duplicationRate, Parameter lossRate, Parameter transferRate, boolean normalize, Taxa guestTaxa, String hostAttributeName) {
		
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
				
		dltrsModel = new DLTRSModel(jprimeGuestRBTree, jprimeHostRBTree, reconciliationHelper, new JPrIMEDoubleMapWrapperForBEASTTree(guest), dltProbs, 
				new Continuous1DPDDependent() {
					// Don't want this prior to interfere with the BEAST priors
					public double getCDF(double arg0) { return 1; } // Only function used by the DLTRS model
					public double getCV() { throw new UnsupportedOperationException(); }
					public RealInterval getDomainInterval() { throw new UnsupportedOperationException(); }
					public double getMean() { throw new UnsupportedOperationException(); }
					public double getMedian() { throw new UnsupportedOperationException(); }
					public double getMode() { throw new UnsupportedOperationException(); }
					public double getPDF(double arg0) { throw new UnsupportedOperationException(); }
					public double getProbability(double arg0, double arg1) { throw new UnsupportedOperationException(); }
					public double getQuantile(double arg0) { throw new UnsupportedOperationException(); }
					public double getStandardDeviation() { throw new UnsupportedOperationException(); }
					public double getVariance() { throw new UnsupportedOperationException(); }
					public double sampleValue(PRNG arg0) { throw new UnsupportedOperationException(); }
					public void setMean(double arg0) { throw new UnsupportedOperationException(); }
					public void setStandardDeviation(double arg0) { throw new UnsupportedOperationException(); }
					public void setVariance(double arg0) { throw new UnsupportedOperationException(); }
					public String getName() { throw new UnsupportedOperationException(); }
					public int getNoOfDimensions() { throw new UnsupportedOperationException(); }
					public int getNoOfParameters() { throw new UnsupportedOperationException(); }
					public void cacheAndUpdate(Map<Dependent, ChangeInfo> arg0, boolean arg1) { throw new UnsupportedOperationException(); }
					public void clearCache(boolean arg0) { throw new UnsupportedOperationException(); }
					public Dependent[] getParentDependents() { throw new UnsupportedOperationException(); }
					public void restoreCache(boolean arg0) { throw new UnsupportedOperationException(); }
			
		});

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
