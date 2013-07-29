package org.ithinktree.becky.JPrIMEWrappers;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

@SuppressWarnings("serial")
public class BEASTBranchRatesWrapperForJPrIMEContinuous1DPDDependent extends AbstractBranchRateModel {

	public BEASTBranchRatesWrapperForJPrIMEContinuous1DPDDependent(String arg0) {
		super(arg0);
		// TODO Auto-generated constructor stub
	}

	@Override
	public double getBranchRate(Tree arg0, NodeRef arg1) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected void acceptState() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void handleModelChangedEvent(Model arg0, Object arg1, int arg2) {
		// TODO Auto-generated method stub
		
	}

	@SuppressWarnings("rawtypes")
	@Override
	protected void handleVariableChangedEvent(Variable arg0, int arg1,
			ChangeType arg2) {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void restoreState() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void storeState() {
		// TODO Auto-generated method stub
		
	}


}
