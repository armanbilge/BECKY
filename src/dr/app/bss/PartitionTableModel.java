package dr.app.bss;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.table.AbstractTableModel;
import javax.swing.text.BadLocationException;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class PartitionTableModel extends AbstractTableModel {

	private PartitionDataList dataList;

	public final static int TREE_MODEL_INDEX = 0;
//	public final static int DATA_TYPE_INDEX = 1;
	public final static int FROM_INDEX = 1;
	public final static int TO_INDEX = 2;
	public final static int EVERY_INDEX = 3;
	public final static int DEMOGRAPHIC_MODEL_INDEX = 4;	
	public final static int BRANCH_SUBSTITUTION_MODEL_INDEX = 5;
	public final static int SITE_RATE_MODEL_INDEX = 6;
	public final static int CLOCK_RATE_MODEL_INDEX = 7;
	public final static int FREQUENCY_MODEL_INDEX = 8;
	public final static int ANCESTRAL_SEQUENCE_INDEX = 9;
	
	public static String[] COLUMN_NAMES = { "Tree Model", //
			// "Data Type",
			"From", //
			"To", //
			"Every", //
			"Demographic model", //
			"Branch Substitution Model", //
			"Site Rate Model", //
			"Clock Rate Model", //
			"Frequency Model", //
			"Ancestral Sequence" //
	};

	private static final Class<?>[] COLUMN_TYPES = new Class<?>[] {
			JComboBox.class, //
			// JComboBox.class, //
			Integer.class, //
			Integer.class, // 
			Integer.class, //
			JButton.class, //
			JButton.class, //
			JButton.class, //
			JButton.class, //
			JButton.class, //
			JButton.class //
	};

	private DemographicModelEditor demographicModelEditor;
	private BranchSubstitutionModelEditor branchSubstitutionModelEditor;
	private SiteRateModelEditor siteRateModelEditor;
	private ClockRateModelEditor clockRateModelEditor;
	private FrequencyModelEditor frequencyModelEditor;
	private AncestralSequenceEditor ancestralSequenceEditor;
	
	public PartitionTableModel(PartitionDataList dataList) {
		this.dataList = dataList;
	}// END: Constructor

	public void addRow(PartitionData row) {
		dataList.add(row);
		fireTableRowsInserted(dataList.size() - 1, dataList.size() - 1);
	}

	public void addDefaultRow() {
		dataList.add(new PartitionData());
		fireTableRowsInserted(dataList.size() - 1, dataList.size() - 1);
	}

	public void deleteRow(int row) {
		dataList.remove(row);
		fireTableDataChanged();
	}

	@Override
	public int getColumnCount() {
		return COLUMN_NAMES.length;
	}// END: getColumnCount

	@Override
	public int getRowCount() {
		return dataList.size();
	}// END: getRowCount

	@Override
	public Class<?> getColumnClass(int columnIndex) {
		return COLUMN_TYPES[columnIndex];
	}// END: getColumnClass

	public boolean isCellEditable(int row, int column) {
		switch (column) {
		case TREE_MODEL_INDEX:
			return true;
		case DEMOGRAPHIC_MODEL_INDEX:
			return false;
//		case DATA_TYPE_INDEX:
//			return true;
		case FROM_INDEX:
			return true;
		case TO_INDEX:
			return true;
		case EVERY_INDEX:
			return true;
		case BRANCH_SUBSTITUTION_MODEL_INDEX:
			return false;
		case SITE_RATE_MODEL_INDEX:
			return false;
		case CLOCK_RATE_MODEL_INDEX:
			return false;
		case FREQUENCY_MODEL_INDEX:
			return false;
		case ANCESTRAL_SEQUENCE_INDEX:
			return false;
		default:
			return false;
		}
	}// END: isCellEditable

	public String getColumnName(int column) {
		return COLUMN_NAMES[column];
	}// END: getColumnName

	@Override
	public Object getValueAt(final int row, final int column) {
		switch (column) {

		case TREE_MODEL_INDEX:
			
			return dataList.get(row).record == null ? new TreesTableRecord() : dataList
					.get(row).record;
			
		case DEMOGRAPHIC_MODEL_INDEX:

			final JButton topologyButton = new JButton(COLUMN_NAMES[column]);
			topologyButton
					.addActionListener(new ListenOpenDemographicModelEditor(row));
			return topologyButton;
			
//		case DATA_TYPE_INDEX:
//			return PartitionData.dataTypes[dataList.get(row).dataTypeIndex];
		
		case FROM_INDEX:
			return dataList.get(row).from;
		
		case TO_INDEX:
			return dataList.get(row).to;
		
		case EVERY_INDEX:
			return dataList.get(row).every;
		
		case BRANCH_SUBSTITUTION_MODEL_INDEX:

			final JButton branchSubstModelButton = new JButton(
					COLUMN_NAMES[column]);
			branchSubstModelButton
					.addActionListener(new ListenOpenBranchSubstitutionModelEditor(
							row));
			return branchSubstModelButton;

		case SITE_RATE_MODEL_INDEX:

			final JButton siteRateModelButton = new JButton(
					COLUMN_NAMES[column]);
			siteRateModelButton
					.addActionListener(new ListenOpenSiteRateModelEditor(row));

			return siteRateModelButton;

		case CLOCK_RATE_MODEL_INDEX:

			final JButton clockRateModelButton = new JButton(
					COLUMN_NAMES[column]);
			clockRateModelButton
					.addActionListener(new ListenOpenClockRateModelEditor(row));
			return clockRateModelButton;

		case FREQUENCY_MODEL_INDEX:

			final JButton frequencyModelButton = new JButton(
					COLUMN_NAMES[column]);
			frequencyModelButton
					.addActionListener(new ListenOpenFrequencyModelEditor(row));
			return frequencyModelButton;
			
		case ANCESTRAL_SEQUENCE_INDEX:
			
			final JButton ancestralSequenceButton = new JButton(COLUMN_NAMES[column]);
			ancestralSequenceButton.addActionListener(new ListenOpenAncestralSequenceEditor(row));
			return ancestralSequenceButton;
		default:
			return "Error";
		}

	}// END: getValueAt

	public void setValueAt(Object value, int row, int column) {

		switch (column) {

		case TREE_MODEL_INDEX:
			dataList.get(row).record = (TreesTableRecord) value;
			break;

//		case DATA_TYPE_INDEX:
//			dataList.get(row).dataTypeIndex = (Integer) Utils.arrayIndex(
//					PartitionData.dataTypes, (String) value);
//			break;

		case FROM_INDEX:
			dataList.get(row).from = (Integer) value;
			break;

		case TO_INDEX:
			dataList.get(row).to = (Integer) value;
			break;

		case EVERY_INDEX:
			dataList.get(row).every = (Integer) value;
			break;

		case BRANCH_SUBSTITUTION_MODEL_INDEX:
			dataList.get(row).substitutionModelIndex = (Integer) value;
			break;

		case SITE_RATE_MODEL_INDEX:
			dataList.get(row).siteRateModelIndex = (Integer) value;
			break;

		case CLOCK_RATE_MODEL_INDEX:
			dataList.get(row).clockModelIndex = (Integer) value;
			break;

		case FREQUENCY_MODEL_INDEX:
			dataList.get(row).frequencyModelIndex = (Integer) value;
			break;

		default:
			throw new RuntimeException("Invalid index.");
		}

		fireTableCellUpdated(row, column);

	}// END: setValueAt

	private class ListenOpenDemographicModelEditor implements ActionListener {

		private int row;

		public ListenOpenDemographicModelEditor(int row) {
			this.row = row;
		}// END: Constructor

		public void actionPerformed(ActionEvent ev) {

			demographicModelEditor = new DemographicModelEditor(dataList, row);
			demographicModelEditor.launch();

		}// END: actionPerformed
	}// END: ListenOpenDemographicModelEditor
	
	private class ListenOpenBranchSubstitutionModelEditor implements
			ActionListener {

		private int row;

		public ListenOpenBranchSubstitutionModelEditor(int row) {
			this.row = row;
		}// END: Constructor

		public void actionPerformed(ActionEvent ev) {

			branchSubstitutionModelEditor = new BranchSubstitutionModelEditor(
					dataList, row);
			branchSubstitutionModelEditor.launch();

		}// END: actionPerformed
	}// END: ListenOpenBranchSubstitutionModelEditor

	private class ListenOpenSiteRateModelEditor implements ActionListener {
		private int row;

		public ListenOpenSiteRateModelEditor(int row) {
			this.row = row;
		}

		public void actionPerformed(ActionEvent ev) {

			try {
				siteRateModelEditor = new SiteRateModelEditor(dataList, row);
				siteRateModelEditor.launch();
			} catch (NumberFormatException e) {
				Utils.handleException(e);
			} catch (BadLocationException e) {
				Utils.handleException(e);
			}

		}// END: actionPerformed
	}// END: ListenOpenSiteRateModelEditor

	private class ListenOpenClockRateModelEditor implements ActionListener {
		private int row;

		public ListenOpenClockRateModelEditor(int row) {
			this.row = row;
		}

		public void actionPerformed(ActionEvent ev) {

			clockRateModelEditor = new ClockRateModelEditor(dataList, row);
			clockRateModelEditor.launch();

		}// END: actionPerformed
	}// END: ListenOpenSiteRateModelEditor

	private class ListenOpenFrequencyModelEditor implements ActionListener {
		private int row;

		public ListenOpenFrequencyModelEditor(int row) {
			this.row = row;
		}

		public void actionPerformed(ActionEvent ev) {

			frequencyModelEditor = new FrequencyModelEditor(dataList, row);
			frequencyModelEditor.launch();

		}// END: actionPerformed
	}// END: ListenOpenSiteRateModelEditor

	
	private class ListenOpenAncestralSequenceEditor implements ActionListener {
		private int row;

		public ListenOpenAncestralSequenceEditor(int row) {
			this.row = row;
		}

		public void actionPerformed(ActionEvent ev) {

			ancestralSequenceEditor = new AncestralSequenceEditor(dataList, row);
			ancestralSequenceEditor.launch();

		}// END: actionPerformed
	}// END: ListenOpenAncestralSequenceEditor
	
	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}

}// END: class
