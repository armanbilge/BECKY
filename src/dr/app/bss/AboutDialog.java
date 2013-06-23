package dr.app.bss;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class AboutDialog extends JDialog {

	private static final int WIDTH = 700;
	private static final int HEIGHT = 650;
	private static final int FONT_SIZE = 15;

	private static final String FILIP_BIELEJEC = "Filip Bielejec";
	private static final String ANDREW_RAMBAUT = "Andrew Rambaut";
	private static final String MARC_SUCHARD = "Marc A. Suchard";
	private static final String PHILIPPE_LEMEY = "Philippe Lemey";
	private static final String LUIZ_MAX_CARVAHLO = "Luiz Max Carvahlo";
	private static final String GUY_BAELE = "Guy Baele";

	public AboutDialog() {
		initUI();
	}// END: Constructor

	public final void initUI() {

		JLabel label;
		JLabel contact;
		JLabel website;
		String addres;

		setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
		getContentPane().setBackground(Color.WHITE);
		setLocationRelativeTo(Utils.getActiveFrame());
		
		add(Box.createRigidArea(new Dimension(0, 10)));

		// Setup image
		label = new JLabel(Utils.createImageIcon(Utils.BSS_ICON));
		label.setAlignmentX(0.5f);
		add(label);

		add(Box.createRigidArea(new Dimension(0, 10)));

		// Setup name
		label = new JLabel(BeagleSequenceSimulatorApp.SHORT_NAME);
		label.setFont(new Font("Serif", Font.BOLD, FONT_SIZE));
		label.setAlignmentX(0.5f);
		add(label);

		add(Box.createRigidArea(new Dimension(0, 10)));

		// Setup long name
		label = new JLabel(BeagleSequenceSimulatorApp.LONG_NAME);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 2));
		label.setAlignmentX(0.5f);
		add(label);

		// Setup version
		label = new JLabel("Version v" + BeagleSequenceSimulatorApp.VERSION
				+ " Prerelease" + ", " + BeagleSequenceSimulatorApp.DATE);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 2));
		label.setAlignmentX(0.5f);
		add(label);

		add(Box.createRigidArea(new Dimension(0, 10)));

		// Setup authors
		label = new JLabel("by " + FILIP_BIELEJEC + ", " + ANDREW_RAMBAUT
				+ ", " + MARC_SUCHARD + ", " + GUY_BAELE + ", "  + LUIZ_MAX_CARVAHLO + " and "
				+ PHILIPPE_LEMEY);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 2));
		label.setAlignmentX(0.5f);
		add(label);

		add(Box.createRigidArea(new Dimension(0, 10)));

		// Setup about
		label = new JLabel("BEAST auxiliary software package");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		website = new JLabel();
		addres = "http://beast.bio.ed.ac.uk";
		website.setText("<html><p><a href=\"" + addres + "\">" + addres
				+ "</a></p></html>");
		website.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		website.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		// website.setVerticalAlignment(SwingConstants.CENTER);
		// website.setHorizontalAlignment(SwingConstants.CENTER);
		website.addMouseListener(new ListenBrowse(addres));
		add(website);

		add(Box.createRigidArea(new Dimension(0, 10)));

		label = new JLabel("Designed and developed by");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		label = new JLabel(FILIP_BIELEJEC + ", " + MARC_SUCHARD + " and " + ANDREW_RAMBAUT);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		add(Box.createRigidArea(new Dimension(0, 10)));

		label = new JLabel("Computational and Evolutionary Virology");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		contact = new JLabel();
		addres = "filip.bielejec@rega.kuleuven.be";
		contact.setText("<html><center><p><a href=\"mailto:" + addres + "\">"
				+ addres + "</a></p></center></html>");
		contact.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		contact.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		// contact.setAlignmentX(0.0f);
		contact.addMouseListener(new ListenSendMail(addres));
		add(contact);

		add(Box.createRigidArea(new Dimension(0, 10)));

		label = new JLabel(
				"Institute of Evolutionary Biology, University of Edinburgh");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		contact = new JLabel();
		addres = "a.rambaut@ed.ac.uk";
		contact.setText("<html><p><a href=\"mailto:" + addres + "\">" + addres
				+ "</a></p></html>");
		contact.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		contact.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		// contact.setAlignmentX(0.5f);
		contact.addMouseListener(new ListenSendMail(addres));
		add(contact);

		add(Box.createRigidArea(new Dimension(0, 10)));

		label = new JLabel("Source code distributed under the GNU LGPL");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		website = new JLabel();
		addres = "http://code.google.com/p/beast-mcmc";
		website.setText("<html><p><a href=\"" + addres + "\">" + addres
				+ "</a></p></html>");
		website.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		website.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		// website.setAlignmentX(0.5f);
		website.addMouseListener(new ListenBrowse(addres));
		add(website);

		add(Box.createRigidArea(new Dimension(0, 20)));

		label = new JLabel("In case of any problems please contact your local witch doctor or a shaman.");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);
		
		add(Box.createRigidArea(new Dimension(0, 20)));

		JButton close = new JButton("Close");
		close.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent event) {
				dispose();
			}
		});

		close.setAlignmentX(0.5f);
		add(close);

		setModalityType(ModalityType.APPLICATION_MODAL);

		setTitle("About " + BeagleSequenceSimulatorApp.SHORT_NAME);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setLocationRelativeTo(null);
		setSize(WIDTH, HEIGHT);
		setResizable(false);
	}// END: initUI

	private class ListenSendMail extends MouseAdapter {

		private String addres;

		public ListenSendMail(String addres) {
			this.addres = addres;
		}

		@Override
		public void mouseClicked(MouseEvent ev) {
			try {

				Desktop.getDesktop().mail(new URI("mailto:" + addres));

			} catch (IOException e) {
				Utils.handleException(
						e,
						"Problem occurred while trying to open this address in your system's standard email client.");
			} catch (URISyntaxException e) {
				Utils.handleException(
						e,
						"Problem occurred while trying to open this address in your system's standard email client.");
			}// END: try-catch block

		}// END: mouseClicked

	}// END: ListenSendMail

	private class ListenBrowse extends MouseAdapter {

		private String website;

		public ListenBrowse(String website) {
			this.website = website;
		}

		@Override
		public void mouseClicked(MouseEvent ev) {
			try {

				Desktop.getDesktop().browse(new URI(website));

			} catch (IOException e) {
				Utils.handleException(
						e,
						"Problem occurred while trying to open this link in your system's standard browser.");
			} catch (URISyntaxException e) {
				Utils.handleException(
						e,
						"Problem occurred while trying to open this link in your system's standard browser.");
			}// END: try-catch block

		}// END: mouseClicked

	}// END: ListenSendMail

}// END: class