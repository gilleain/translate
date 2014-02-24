package translation.hbondviewer;

import java.awt.Graphics;
import javax.swing.JPanel;
import translation.model.Chain;

public class HBondViewPanel extends JPanel {
    private int width;
    private int height;
    private HBondDiagramDrawer hBondDiagramDrawer;
    
    public HBondViewPanel() {
        this.width = 0;
        this.height = 0;
        this.hBondDiagramDrawer = null;
    }

    public HBondViewPanel(int width, int height) {
        this.width = width;
        this.height = height;
        this.hBondDiagramDrawer = new HBondDiagramDrawer(this.width, this.height);
    }

    public void setChain(Chain chain) {
        this.hBondDiagramDrawer.setChain(chain);
        this.repaint();
    }

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        if ((this.width == 0) || (this.height == 0)) {
            this.width = this.getWidth();
            this.height = this.getHeight();
        } else {
            if (this.hBondDiagramDrawer == null) {
                this.hBondDiagramDrawer = new HBondDiagramDrawer(this.width, this.height);
            }   
            this.hBondDiagramDrawer.paint(g);
        }
    }

}
