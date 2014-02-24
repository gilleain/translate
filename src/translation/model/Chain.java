package translation.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import translation.Axis;
import translation.ChainDomainMap;
import translation.Geometer;

public class Chain {
    private String label;
    private String type;
    private Point3d center;
    private List<Residue> residues;
    private List<HBond> hbonds;
    private List<Sheet> sheets;
    private List<BackboneSegment> backboneSegments;
    private List<Edge> chiralities;
    private List<Domain> domains;

    public Chain() {
        this.residues = new ArrayList<Residue>();
        this.hbonds = new ArrayList<HBond>();
        this.sheets = new ArrayList<Sheet>();
        this.backboneSegments = new ArrayList<BackboneSegment>();
        this.chiralities = new ArrayList<Edge>();
        this.domains = new ArrayList<Domain>();
        this.center = null;
    }

    public Chain(String label) {
        this();
        this.label = label;
    }

    public int length() {
        return this.residues.size();
    }
    
    public Point3d getCenter() {
    	return this.center;
    }
    
    public String getType() {
    	return this.type;
    }
    
    public List<Domain> getDomains() {
    	return this.domains;
    }

    public boolean isDNA() {
        return ((Residue) this.residues.get(0)).isDNA();
    }

    // since this doesn't check, in the event that a structure actually HAS
    // backbone amide hydrogens, it will overwrite them...
    public void addBackboneAmideHydrogens() {
        for (int i = 1; i < this.residues.size(); i++) {
            Residue lastResidue = (Residue) this.residues.get(i - 1);
            Residue residue     = (Residue) this.residues.get(i);

            if (residue.isPro()) {
                continue;
            }

            Point3d lastC       = lastResidue.getCoordinates("C");
            Point3d lastO       = lastResidue.getCoordinates("O");
            Point3d N           =     residue.getCoordinates("N");

            Vector3d OC = new Vector3d();
            OC.sub(lastC, lastO);
            OC.normalize();

            Point3d H = new Point3d();
            H.add(N, OC);

            residue.setAtom("H", H);
        }
    } 

    public Residue createResidue(int pdbNumber, String residueType) {
        Residue r = new Residue(this.residues.size(), pdbNumber, residueType);
        this.residues.add(r);
        return r;
    }

    public void createHelix(int helixStartIndex, int helixEndIndex) {
        BackboneSegment helix = new Helix(this.getResidueByAbsoluteNumbering(helixStartIndex));
        for (int i = helixStartIndex + 1; i < helixEndIndex + 1; i++) {
            helix.expandBy(this.getResidueByAbsoluteNumbering(i));
        }
        //System.err.println("Created " + helix + " from " + helixStartIndex + " " + helixEndIndex);
        this.backboneSegments.add(helix);
    }

    public void createStrand(int strandStartIndex, int strandEndIndex) {
        BackboneSegment strand = new Strand(this.getResidueByAbsoluteNumbering(strandStartIndex));
        for (int i = strandStartIndex + 1; i < strandEndIndex + 1; i++) {
            strand.expandBy(this.getResidueByAbsoluteNumbering(i));
        }
        //System.err.println("Created " + strand + " from " + strandStartIndex + " " + strandEndIndex);
        this.backboneSegments.add(strand);
    }

    public void createLoop(int startIndex, int endIndex) {
        BackboneSegment unstructured = new UnstructuredSegment(this.getResidueByAbsoluteNumbering(startIndex));
        for (int i = startIndex + 1; i < endIndex + 1; i++) {
            unstructured.expandBy(this.getResidueByAbsoluteNumbering(i));
        }
        this.backboneSegments.add(unstructured);
    }

    // firstly, assume that the sses are sorted so that for i < j; i.end < j.end
    // secondly assume that, if two helices overlap, they should be merged
    public void mergeHelices() {
        for (int i = 0; i < this.backboneSegments.size() - 1; i++) {
            BackboneSegment sseA = (BackboneSegment) this.backboneSegments.get(i);
            if (sseA instanceof Helix) {
                BackboneSegment sseB = (BackboneSegment) this.backboneSegments.get(i + 1);
                if (sseB instanceof Helix) {
                    if (sseA.overlaps(sseB)) {
                        //System.err.println("Merging : " + sseA + " and " + sseB); 
                        sseB.mergeWith(sseA);

                        //System.err.println("Removing: " + sseA); 
                        this.backboneSegments.remove(i);
                        i--;
                    } else {
                        //System.err.println("No overlap : " + sseA + " and " + sseB); 
                    }
                } 
            }
        }
    }

    public void addTerminii() {
        if (this.backboneSegments.size() == 0) {
            return;
        }

        if (!((BackboneSegment) this.backboneSegments.get(0) instanceof Terminus)) {
            this.backboneSegments.add(0, new Terminus("N Terminus", 'N'));
        }

        if (!((BackboneSegment) this.backboneSegments.get(this.backboneSegments.size() - 1) instanceof Terminus)) {
            this.backboneSegments.add(new Terminus("C Terminus", 'C'));
        }
    }

    public void sortBackboneSegments() {
        Collections.sort(this.backboneSegments);
    }

    public Iterator<BackboneSegment> backboneSegmentIterator() {
        return this.backboneSegments.iterator();
    }

    public boolean hasResidueByPDBNumbering(int pdbResidueNumber) {
        for (Residue residue : this.residues) {
            if (residue.getPDBNumber() == pdbResidueNumber) {
                return true;
            }
        }
        return false;
    }

    public Residue getResidueByPDBNumbering(int pdbResidueNumber) {
    	for (Residue residue : this.residues) {
            if (residue.getPDBNumber() == pdbResidueNumber) {
                return residue;
            }
        }
        return null;
    }

    public Residue getResidueByAbsoluteNumbering(int i) {
        return this.residues.get(i);
    }

    public boolean hasResidueByAbsoluteNumbering(int i) {
        return i < this.residues.size();
    }

    public Residue firstResidue() {
        return this.residues.get(0);
    }

    public Residue lastResidue() {
        return this.residues.get(this.residues.size() - 1);
    }

    public Residue getNextResidue(int i) {
        if (i + 1 >= this.residues.size()) {
            return null;
        } else {
            return this.residues.get(i + 1);
        }
    }

    public Iterator<Residue> residueIterator() {
        return this.residues.iterator();
    }

    public Iterator<Residue> residueIterator(int i) {
        return this.residues.subList(i, this.residues.size() - 1).iterator();
    }

    public Iterator<HBond> hbondIterator() {
        Collections.sort(this.hbonds);
        return this.hbonds.iterator();
    }

    public void addSheet(Sheet sheet) {
        this.sheets.add(sheet);
    }

    public void createSheet(BackboneSegment strand, BackboneSegment otherStrand) {
        if (this.sheets.size() == 0) {
            this.sheets.add(new Sheet(1, strand, otherStrand));
        } else {
            int lastSheetNumber = ((Sheet) this.sheets.get(this.sheets.size() - 1)).getNumber();
            this.sheets.add(new Sheet(lastSheetNumber + 1, strand, otherStrand));
        }
    }

    public void removeSheet(Sheet sheet) {
        this.sheets.remove(sheet);
    }

    public Sheet getSheetContaining(BackboneSegment strand) {
    	for (Sheet sheet : this.sheets) {
            if (sheet.contains(strand)) {
                return sheet;
            }
        }
        return null;
    }

    public Iterator<Sheet> sheetIterator() {
        return this.sheets.iterator();
    }

    public int numberOfSheets() {
        return this.sheets.size();
    }

    public void addChirality(BackboneSegment first, BackboneSegment second, char chirality) {
        this.chiralities.add(new Edge(first, second, chirality));
    }

    public Axis getAxis() {
        List<Point3d> centroids = new ArrayList<Point3d>();
        for (BackboneSegment nextBackboneSegment : this.backboneSegments) {
            Axis a = nextBackboneSegment.getAxis();
            //System.out.println("Centroid for BS " + nextBackboneSegment + " is " + a.getCentroid());
            centroids.add(a.getCentroid());
        }
        //System.out.println("averagePoints = " + Geometer.averagePoints(centroids));
        Axis axis = Geometer.leastSquareAxis(centroids);

        //this is a slight hack to get start and end coordinates into the chain axis
        //might be better to use the N terminus as the start?
        Point3d start = new Point3d(axis.getCentroid());
        start.scale(5, axis.getAxisVector());
        axis.setStart(start);
        Point3d end = new Point3d(axis.getCentroid());
        end.scale(-5, axis.getAxisVector());
        axis.setEnd(end);

        return axis;
    }

    public void addHBond(HBond hbond) {
        this.hbonds.add(hbond);
    }

    public void setLabel(String label) {
        this.label = label;
    }

    public String getCathCompatibleLabel() {
        if (this.label.equals("  ")) {
            return "0";
        } else {
            return this.label.trim();
        }
    }

    public String getLabel() {
        return this.label;
    }

    public void addBackboneSegment(BackboneSegment backboneSegment) {
        this.backboneSegments.add(backboneSegment);
    }

    public ListIterator<BackboneSegment> backboneSegmentListIterator() {
        return this.backboneSegments.listIterator();
    }

    public ListIterator<BackboneSegment> backboneSegmentListIterator(int i) {
        return this.backboneSegments.listIterator(i);
    }

    public ListIterator<BackboneSegment> backboneSegmentListIterator(BackboneSegment backboneSegment) {
        return this.backboneSegments.listIterator(this.backboneSegments.indexOf(backboneSegment));
    }

    public ListIterator<BackboneSegment> backboneSegmentListIterator(BackboneSegment startSegment, BackboneSegment endSegment) {
        return this.backboneSegments.subList(this.backboneSegments.indexOf(startSegment) + 1, this.backboneSegments.indexOf(endSegment)).listIterator();
    }

    public int numberOfSegmentsInBetween(BackboneSegment startSegment, BackboneSegment endSegment) {
        return this.backboneSegments.indexOf(endSegment) - (this.backboneSegments.indexOf(startSegment) + 1);
    }
    
    public List<BackboneSegment> filterBackboneSegmentsByDomain(Domain domain) {
        return domain.filter(this.backboneSegments);
    }

    public void calculateTorsions() {
        Iterator<Residue> residueIterator = this.residueIterator();
        Residue previousResidue = null;
        Residue thisResidue = null;

        if (residueIterator.hasNext()) {
            previousResidue = (Residue) residueIterator.next();
        }

        if (residueIterator.hasNext()) {
            thisResidue = (Residue) residueIterator.next();
        }

        while (residueIterator.hasNext()) {
            Point3d prevC = previousResidue.getCoordinates("C");
            Point3d thisN = thisResidue.getCoordinates("N");
            Point3d thisCA = thisResidue.getCoordinates("CA");
            Point3d thisC = thisResidue.getCoordinates("C");

            double phi = Geometer.torsion(prevC, thisN, thisCA, thisC);
            thisResidue.setPhi(phi);

            previousResidue = thisResidue;
            if (residueIterator.hasNext()) {
                Residue nextResidue = (Residue) residueIterator.next();

                Point3d nextN = nextResidue.getCoordinates("N");

                double psi = Geometer.torsion(thisN, thisCA, thisC, nextN);
                thisResidue.setPsi(psi);

                thisResidue = nextResidue;
            } else {
                break;
            }
        }
    }

    public void findOrientations() {
        Axis chainAxis = this.getAxis();

        for (BackboneSegment nextSegment : this.backboneSegments) {
            nextSegment.determineOrientation(chainAxis);
        }
    }

    public String toPymolScript() {
        StringBuffer script = new StringBuffer();

        int segmentNumber = 1;
        for (BackboneSegment bs : this.backboneSegments) {
            String name = bs.getTopsSymbol() + "" + segmentNumber + "(" + bs.firstPDB() + "-" + bs.lastPDB() + ")";
            String selection = "not hetatm and name ca+c+n+o+h and resi " + bs.firstPDB() + "-" + bs.lastPDB();
            script.append("cmd.select(\"" + name + "\", \"" + selection + "\")\n");
            String color = "green";
            if (bs instanceof Strand) {
                script.append("line(");
//                Axis axis = bs.getAxis();
//                Point3d start = axis.getStart();
//                Point3d end = axis.getEnd();
//                script.append(String.format("%6.2f, %6.2f, %6.2f, ", start.x, start.y, start.z));
//                script.append(String.format("%6.2f, %6.2f, %6.2f, ", end.x, end.y, end.z));
                script.append("\"").append(name).append("axis\")\n");
                color = "yellow";
            } else if (bs instanceof Helix) {
                color = "red";
            }
            script.append("cmd.color(\"" + color + "\", \"" + name + "\")\n");
            segmentNumber++;
        }
        return script.toString(); 
    }

    public Map<String, String> toTopsDomainStrings(ChainDomainMap chainDomainMap) {
        if (chainDomainMap != null && chainDomainMap.containsKey(this.getCathCompatibleLabel())) {
            List<Domain> domains = chainDomainMap.get(this.getCathCompatibleLabel());
            return this.toTopsDomainStrings(domains);
        } else {
            Map<String, String> h = new HashMap<String, String>();
            h.put("0", this.toTopsString(new Domain(0)));
            return h;
        }
    }

    public Map<String, String> toTopsDomainStrings(List<Domain> domains) {
        Map<String, String> domainStrings = new HashMap<String, String>();

        for (Domain d : domains) {
            //System.err.println("Getting tops string for domain " + d);
            domainStrings.put(d.getID(), this.toTopsString(d));
        }

        return domainStrings;
    }

    public String toTopsString(Domain domain) {
        //name
        StringBuffer s = new StringBuffer();
        s.append(this.getCathCompatibleLabel() + domain.getNumber() + " ");

        //vertexstring
        Iterator<BackboneSegment> backboneSegmentIterator;
        if (domain.isEmpty()) {
            backboneSegmentIterator = this.backboneSegments.iterator();
        } else {
            List<BackboneSegment> segmentsFilteredByDomain = domain.filter(this.backboneSegments);
            segmentsFilteredByDomain.add(0, new Terminus("N Terminus", 'N'));
            segmentsFilteredByDomain.add(new Terminus("C Terminus", 'C'));
            backboneSegmentIterator = segmentsFilteredByDomain.iterator();
        }

        int vertexNumber = 0;
        while (backboneSegmentIterator.hasNext()) {
            BackboneSegment nextBackboneSegment = (BackboneSegment) backboneSegmentIterator.next();
            if (nextBackboneSegment instanceof UnstructuredSegment) {
                continue;
            } else {
                nextBackboneSegment.setNumber(vertexNumber);
                s.append(nextBackboneSegment.getTopsSymbol());
                vertexNumber++;
            }
        }

        s.append(" ");
        
        //edgestring
        List<Edge> edges = new ArrayList<Edge>();
        for (Sheet sheet : this.sheets) {
            edges.addAll(sheet.toTopsEdges(domain));
        }

        Collections.sort(edges);
        
        //merge the chirals with the hbonds
        Iterator<Edge> chiralIterator = this.chiralities.iterator();
        Iterator<Edge> edgeIterator;
        while (chiralIterator.hasNext()) {
            boolean merged = false;
            Edge chiral = chiralIterator.next();
            edgeIterator = edges.iterator();
            while (edgeIterator.hasNext()) {
                Edge hbond = edgeIterator.next();
                if (hbond.equals(chiral)) {
                    hbond.mergeWith(chiral);
                    merged = true;
                    break;
                }
            }
            // those chirals that are not merged are simply added
            if (!merged && chiral.containedIn(domain)) {
                edges.add(chiral);
            } 
        }

        Collections.sort(edges);

        edgeIterator = edges.iterator();
        while (edgeIterator.hasNext()) {
            Edge edge = (Edge) edgeIterator.next();
            s.append(edge.toString());
        }

        return s.toString();
    }

    public String toString() {
        StringBuffer s = new StringBuffer();
        s.append("Chain : " + this.getLabel() + " residue " + this.firstResidue().getPDBNumber() + " to " + this.lastResidue().getPDBNumber() + "\n");

        for (Residue r : this.residues) {
            s.append(r.toFullString());
            s.append("\n");
        }

        for (HBond nextBond : this.hbonds) {
            s.append(nextBond + "\n");
        }

        for (BackboneSegment nextBackboneSegment : this.backboneSegments) {
            s.append(nextBackboneSegment.toFullString() + "\n");
        }

        for (Sheet sheet : this.sheets) {
            System.out.println("Sheet paths " + sheet.getSheetPaths());
            s.append(sheet + "\n");
        }

        return s.toString();
    }
}
