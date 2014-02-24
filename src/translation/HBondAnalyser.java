package translation;

import java.io.FileInputStream;
import java.io.InputStream;
import java.io.IOException;
import java.io.OutputStream;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

import javax.vecmath.Point3d;

import translation.model.Chain;
import translation.model.HBond;
import translation.model.Protein;
import translation.model.Residue;

public class HBondAnalyser {

    private Properties properties;

    private int threeTenHelixStart;
    private int threeTenHelixEnd;

    private int alphaHelixStart;
    private int alphaHelixEnd;

    private int piHelixStart;
    private int piHelixEnd;

    private int strandStart;
    private int strandEnd;

    public HBondAnalyser() {
        this.properties = new Properties();
    }

    public HBondAnalyser(Properties properties) {
        this.properties = properties;
    }

    public void setDefaultProperties() {
        this.properties.setProperty("CALCULATE_BACKBONE_AMIDE_HYDROGENS", "TRUE");
        this.properties.setProperty("MAX_HO_DISTANCE", "3.0");
        this.properties.setProperty("MIN_NHO_ANGLE", "120");
        this.properties.setProperty("MIN_HOC_ANGLE", "90");
    }

    public void resetEndpoints() {
        this.threeTenHelixStart = -1;
        this.threeTenHelixEnd   = -1;

        this.alphaHelixStart    = -1;
        this.alphaHelixEnd      = -1;

        this.piHelixStart       = -1;
        this.piHelixEnd         = -1;

        this.strandStart        = -1;
        this.strandEnd          = -1;
    }

    public String endpoints() {
        //return String.format("%d %d %d %d %d %d %d %d", threeTenHelixStart, threeTenHelixEnd, alphaHelixStart, alphaHelixEnd, piHelixStart, piHelixEnd, strandStart, strandEnd);
        return "";
    }

    public void setProperty(String key, String value) {
        this.properties.setProperty(key, value);
    }

    public void analyse(Protein protein) throws PropertyException {
        Iterator<Chain> chains = protein.chainIterator();
        boolean calculateBackboneHydrogens = this.properties.getProperty("CALCULATE_BACKBONE_AMIDE_HYDROGENS").equals("TRUE");

        while (chains.hasNext()) {
            Chain chain = (Chain) chains.next();

            if (calculateBackboneHydrogens) {
                chain.addBackboneAmideHydrogens();
            }
            this.analyse(chain);
        }
    }

    public void analyse(Chain chain) throws PropertyException {
        double MAX_HO_DISTANCE = 0.0;
        double MIN_NHO_ANGLE   = 0.0;
        double MIN_HOC_ANGLE   = 0.0;

        try {
            MAX_HO_DISTANCE = Double.parseDouble(this.properties.getProperty("MAX_HO_DISTANCE"));
            MIN_NHO_ANGLE   = Double.parseDouble(this.properties.getProperty("MIN_NHO_ANGLE"));
            MIN_HOC_ANGLE   = Double.parseDouble(this.properties.getProperty("MIN_HOC_ANGLE"));
            //System.out.println("HO " + MIN_HO_DISTANCE + " NHO " + MIN_NHO_ANGLE + " HOC " + MIN_HOC_ANGLE);
        } catch (NumberFormatException nfe) {
            throw new PropertyException("Error in properties!");
        }

        int index = -1;
        Iterator<Residue> residues = chain.residueIterator();

        while (residues.hasNext()) {
            Residue first = (Residue) residues.next();
            index++;

            // we ignore gamma turns!
            int position = first.getAbsoluteNumber();
            int nextPosition = position + 3;

            Point3d firstN = first.getCoordinates("N");
            Point3d firstH = first.getCoordinates("H");
            Point3d firstO = first.getCoordinates("O");
            Point3d firstC = first.getCoordinates("C");

            // FIXME : unfortunately, this misses out on PRO residues (also below)
            if (firstN == null || firstH == null || firstO == null || firstC == null) {
                continue;
            }

            // allow for chain breaks, or return if we have reached the end
            if (!chain.hasResidueByAbsoluteNumbering(nextPosition)) {
                Residue second = chain.getNextResidue(position);
                if (second == null) {   //probably reached the end of the chain!
                    break;
                }
                nextPosition = second.getAbsoluteNumber();
            }

            // now, compare the first residue to the residues further on in the chain
            Iterator<Residue> itr = chain.residueIterator(nextPosition);
            while (itr.hasNext()) {
                int secondPosition = itr.next().getAbsoluteNumber();

                // we must still check this, since a chain break might move us to i + 2
                if (secondPosition < (position + 3)) {
                    continue;
                }

                Residue second;
                try {
                    second = chain.getResidueByAbsoluteNumbering(secondPosition);
                    if (!second.isStandardAminoAcid()) {
                        continue;
                    }
                } catch (IndexOutOfBoundsException i) {
                    break;
                }

                Point3d secondN = second.getCoordinates("N");
                Point3d secondH = second.getCoordinates("H");
                Point3d secondO = second.getCoordinates("O");
                Point3d secondC = second.getCoordinates("C");

                // FIXME : PRO residues...
                if (secondN == null || secondH == null || secondO == null || secondC == null) {
                    continue;
                }

                // bonds from first N-H to second C=O
                double firstHODistance = firstH.distance(secondO);
                double firstNHOAngle = Geometer.angle(firstN, firstH, secondO);
                double firstHOCAngle = Geometer.angle(firstH, secondO, secondC);
            
                HBond firstSecondBond = null;
                if (firstHODistance < MAX_HO_DISTANCE && firstNHOAngle > MIN_NHO_ANGLE && firstHOCAngle > MIN_HOC_ANGLE) {
                    firstSecondBond = new HBond(first, second, firstHODistance, firstNHOAngle, firstHOCAngle);
                }

                if (firstSecondBond != null) {
                    chain.addHBond(firstSecondBond);
                    first.addHBond(firstSecondBond);
                    second.addHBond(firstSecondBond);
                }

                // bonds from second N-H to first C=O
                double secondHODistance = secondH.distance(firstO);
                double secondNHOAngle = Geometer.angle(secondN, secondH, firstO);
                double secondHOCAngle = Geometer.angle(secondH, firstO, firstC);

                HBond secondFirstBond = null;
                if (secondHODistance < MAX_HO_DISTANCE && secondNHOAngle > MIN_NHO_ANGLE && secondHOCAngle > MIN_HOC_ANGLE) {
                    secondFirstBond = new HBond(second, first, secondHODistance, secondNHOAngle, secondHOCAngle);
                }

                if (secondFirstBond != null) {
                    chain.addHBond(secondFirstBond);
                    first.addHBond(secondFirstBond);
                    second.addHBond(secondFirstBond);
                }
            }
            
            // now, use these hbond assignments to determine the residue's environment
            List<String> tags = this.convertBondsToTags(first);
            //System.out.println(first.toFullString() + " " + tags);

            this.updateSSEEndpoints(index, tags, chain);
            //System.out.println(index + " " + first + " " + this.endpoints() + " " + tags);
        }

        // finally, finish off the SSEs
        this.finishSSES(chain);
        chain.sortBackboneSegments();
        chain.mergeHelices();
        chain.addTerminii();
    }

    public List<String> convertBondsToTags(Residue residue) {
    	List<HBond> nTerminalHBonds = residue.getNTerminalHBonds();
    	List<HBond> cTerminalHBonds = residue.getCTerminalHBonds();

        List<String> tags = new ArrayList<String>();

        for (int i = 0; i < nTerminalHBonds.size(); i++) {
            int n = nTerminalHBonds.get(i).getResidueSeparation();
            String nTag = this.convertBondsToTag(n, 0);
            if (nTag != null) {
                tags.add(nTag);
            }
        }

        for (int j = 0; j < cTerminalHBonds.size(); j++) {
            int c = cTerminalHBonds.get(j).getResidueSeparation();
            String cTag = this.convertBondsToTag(0, c);
            if (cTag != null) {
                tags.add(cTag);
            }
        }

        for (int i = 0; i < nTerminalHBonds.size(); i++) {
            int n = nTerminalHBonds.get(i).getResidueSeparation();

            for (int j = 0; j < cTerminalHBonds.size(); j++) {
                int c = cTerminalHBonds.get(j).getResidueSeparation();

                String ncTag = this.convertBondsToTag(n, c);
                if (ncTag != null) {
                    tags.add(ncTag);
                }
            }
        }
        
        return tags;
    }

    public String convertBondsToTag(int n, int c) {

        // neither has a bond
        if (n == 0 && c == 0) {
            return "Loop";
        }

        // one or the other has a bond
        else if (n != 0 && c == 0) {
            if (n == 3) {
                return "End of a 310 helix";
            } else if (n == 4) {
                return "End of an alpha helix";
            } else if (n == 5) {
                return "End of a pi helix";
            } else if (n < -5 || n > 5) {
                return "Single Strand Bond";
            }
        }

        else if (n == 0 && c != 0) {
            if (c == 3) {
                return "Start of a 310 helix";
            } else if (c == 4) {
                return "Start of an alpha helix";
            } else if (c == 5) {
                return "Start of a pi helix";
            } else if (c < -5 || c > 5) {
                return "Single Strand Bond";
            }
        }

        // both have one bond
        else if (n != 0 && c != 0) {

            // a standard : middle of a helix
            if (n == 4 && c == 4) {
                return "Middle of an alpha helix";
            // we counter-intuitively take the SUM here, because one will always be negative, and the other positive 
            } else if ((Math.abs(n - c) == 2) && (n > 5 || n < -5)) {
                return "Parallel Strand";
            } else if ((n == c)  && (n > 5 || n < -5)) {
                return "Antiparallel Strand";
            }
        }

        return null;
    }

    public void updateSSEEndpoints(int index, List<String> tags, Chain chain) {
        if (tags.contains("Start of a 310 helix")) {

            // not seen any three ten bond before
            if (this.threeTenHelixStart == -1) {
                this.threeTenHelixStart = index;

            // this three ten bond does not overlap with the previous one
            } else if (index + 3 > this.threeTenHelixEnd + 1) {

                // only create three-ten helices with more than one bond
                if (this.threeTenHelixEnd - this.threeTenHelixStart > 3) {
                    chain.createHelix(this.threeTenHelixStart + 1, this.threeTenHelixEnd - 1); 
                    //System.err.println("310 : " + this.threeTenHelixStart + ":" +  this.threeTenHelixEnd); 
                }

                // start a potential new helix with this bond
                this.threeTenHelixStart = index;
            }

            // in any case, the helix will start three places further on
            this.threeTenHelixEnd = index + 3;
        }

        if (tags.contains("Start of an alpha helix")) {

            // not seen any alpha bonds before
            if (this.alphaHelixStart == -1) {
                this.alphaHelixStart = index;

            // this three ten bond does not overlap with the previous one
            } else if (index + 4 > this.alphaHelixEnd + 1) {

                // only create three-ten helices with more than one bond
                if (this.alphaHelixEnd - this.alphaHelixStart > 4) {
                    chain.createHelix(this.alphaHelixStart + 1, this.alphaHelixEnd - 1); 
                    //System.err.println("Alpha: " + this.alphaHelixStart + ":" +  this.alphaHelixEnd); 
                }

                // start a potential new helix with this bond
                this.alphaHelixStart = index;
            }

            // in any case, the helix will start three places further on
            this.alphaHelixEnd = index + 4;
        }

        if (tags.contains("Start of a pi helix")) {

            // not seen any three ten bond before
            if (this.piHelixStart == -1) {
                this.piHelixStart = index;

            // this three ten bond does not overlap with the previous one
            } else if (index + 5 > this.piHelixEnd + 1) {

                // only create three-ten helices with more than one bond
                if (this.piHelixEnd - this.piHelixStart > 5) {
                    chain.createHelix(this.piHelixStart + 1, this.piHelixEnd - 1); 
                    //System.err.println("Pi: " + this.piHelixStart + ":" +  this.piHelixEnd); 
                }

                // start a potential new helix with this bond
                this.piHelixStart = index;
            }

            // in any case, the helix will start three places further on
            this.piHelixEnd = index + 5;
        }
        
        if (tags.contains("Antiparallel Strand") || tags.contains("Parallel Strand") || tags.contains("Single Strand Bond")) {
            // System.err.println("Strand like bond at residue " + chain.getResidueByAbsoluteNumbering(index) + " indices = " + this.endpoints() + " tags = " + tags);
            if (this.strandStart == -1) {
                this.strandStart = index;
                this.strandEnd = index;
            } else {
                boolean partnerStrandContinuous = this.residuePartnerNeighboursPrevious(index, chain);
                if (this.strandEnd < index - 2 || !partnerStrandContinuous) {
                    if (this.strandEnd - this.strandStart > 0) {
                        chain.createStrand(this.strandStart, this.strandEnd);
                    }
                    this.strandStart = index;
                    this.strandEnd = index;
                } else {
                    if (this.strandStart != -1) {
                        this.strandEnd = index;
                    }
                }
            }
        }
    }

    // there must be a better way. avoiding tags altogether, storing the n and c terminal hbond separately in Residue,
    // aaaaaaannnnnd all sorts of other optimisations
    public boolean residuePartnerNeighboursPrevious(int index, Chain chain) {
        // 2 or 1?
        if (index < 2) {
            return true;
        }

        // could, of course, just pass the residue along to this method...
        Residue residue = chain.getResidueByAbsoluteNumbering(index);
        Residue twoBehind = chain.getResidueByAbsoluteNumbering(index - 2);

        // get the hbond partners for these two residues as arrays of indices
        int[] residueHBondPartners = residue.getHBondPartners();
        int[] twoBehindHBondPartners = twoBehind.getHBondPartners();

        // this allows bond in the other register (odd/even) to start without starting a new strand
        if (twoBehindHBondPartners.length == 0) {
            return true;
        }

        // check to see if the indices match up
        for (int i = 0; i < residueHBondPartners.length; i++) {
            int iIndex = residueHBondPartners[i];
            for (int j = 0; j < twoBehindHBondPartners.length; j++) {
                int jIndex = twoBehindHBondPartners[j];
                if (this.indicesArePartners(iIndex, jIndex)) {
                    return true;
                }
            }
        }
        return false;
    }

    public boolean indicesArePartners(int i, int j) {
        return (Math.abs(i - j) <= 2);
    }

    public void finishSSES(Chain chain) {
        if (this.threeTenHelixEnd != -1 && this.threeTenHelixEnd - this.threeTenHelixStart > 3) {
            //System.err.println("Last 310 : " + this.threeTenHelixStart + ":" +  this.threeTenHelixEnd); 
            chain.createHelix(this.threeTenHelixStart + 1, this.threeTenHelixEnd - 1);
        }

        if (this.alphaHelixEnd != -1 && this.alphaHelixEnd - this.alphaHelixStart > 4) {
            //System.err.println("Last Alpha: " + this.alphaHelixStart + ":" +  this.alphaHelixEnd); 
            chain.createHelix(this.alphaHelixStart + 1, this.alphaHelixEnd - 1);
        }

        if (this.piHelixEnd != -1 && this.piHelixEnd - this.piHelixStart > 5) {
            //System.err.println("Last Pi: " + this.piHelixStart + ":" +  this.piHelixEnd); 
            chain.createHelix(this.piHelixStart + 1, this.piHelixEnd - 1);
        }

        if (this.strandStart != -1 && this.strandEnd - this.strandStart > 0) {
            chain.createStrand(this.strandStart, this.strandEnd);
        }

        this.resetEndpoints();
    }

    public void storeProperties(OutputStream out) throws IOException {
        String header = "HBond Analysis Parameters";

        this.properties.store(out, header);
    }

    public void loadProperties(InputStream in) throws IOException {
        this.properties.load(in); 
    }

    public String toString() {
        return this.properties.toString();
    }

    public static void main(String[] args) {
        String pdbFilename        = args[0];
        String propertiesFilename = args[1];

        try {
            HBondAnalyser hBondAnalyser = new HBondAnalyser();
            hBondAnalyser.loadProperties(new FileInputStream(propertiesFilename));
            hBondAnalyser.storeProperties(System.err);

            Protein protein = PDBReader.read(pdbFilename);

            hBondAnalyser.analyse(protein);
            System.out.println(protein);
        } catch (IOException ioe) {
            System.err.println(ioe);
        } catch (PropertyException pe) {
            System.err.println(pe);
        }
    }
}
