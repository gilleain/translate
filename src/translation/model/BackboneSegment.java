package translation.model;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Logger;

import javax.vecmath.Point3d;

import translation.Axis;

public abstract class BackboneSegment implements Comparable<BackboneSegment> {
	
    protected int number;
    
    protected SortedSet<Residue> residues;
    
    protected Axis axis;
    
    protected String orientation;

    public BackboneSegment() {
        this.residues = new TreeSet<Residue>();
        this.axis = null;
        this.orientation = "None";
    }

    public BackboneSegment(Residue first) {
        this();
        this.residues.add(first);
    }

    public abstract String toFullString();

    public abstract char getTypeChar();

    public int compareTo(BackboneSegment other) {
        try {
            return this.firstResidue().compareTo(other.firstResidue());
        } catch (NoSuchElementException n) {
            return 0;
        }
    }

    public void setNumber(int number) {
        this.number = number;
    }

    public int getNumber() {
        return this.number;
    }

    public int length() {
        return this.residues.size();
    }

    public int firstPDB() {
        try {
            return this.firstResidue().getPDBNumber();
        } catch (NoSuchElementException n) {
            return -1;
        }
    }

    public int lastPDB() {
        try {
            return this.lastResidue().getPDBNumber();
        } catch (NoSuchElementException n) {
            return -1;
        }
    }

    public Residue firstResidue() throws NoSuchElementException {
        return (Residue) this.residues.first();
    }

    public Residue lastResidue() throws NoSuchElementException {
        return (Residue) this.residues.last();
    }

    public Iterator<Residue> residueIterator() {
        return this.residues.iterator();
    }

    public boolean bondedTo(Residue otherResidue) {
        Iterator<Residue> residueIterator = this.residueIterator();
        while (residueIterator.hasNext()) {
            Residue residue = (Residue) residueIterator.next();
            if (residue.bondedTo(otherResidue)) {
                return true;
            }
        } 
        return false;
    }

    public boolean continuousWith(BackboneSegment other) {
        if (this.getClass() == other.getClass()) {
            if (this.getAxis().approximatelyLinearTo(other.getAxis())) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    public void mergeWith(BackboneSegment other) {
        this.residues.addAll(other.residues);
        this.axis = null;
    }

    public boolean overlaps(BackboneSegment other) {
        try {
            int oS = other.firstResidue().getAbsoluteNumber();
            int oE = other.lastResidue().getAbsoluteNumber();
            return this.containsAbsoluteNumber(oS) || this.containsAbsoluteNumber(oE);
        } catch (NoSuchElementException n) {
            return false;
        }
    }

    public List<Point3d> getCAlphaCoordinates() {
        List<Point3d> cAlphas = new ArrayList<Point3d>();
        for (Residue nextResidue : this.residues) {
            cAlphas.add(nextResidue.getCoordinates("CA"));
        }
        return cAlphas;
    }

    public char getTopsSymbol() {
        if (this.getOrientation().equals("UP")) {
            return Character.toUpperCase(this.getTypeChar());
        } else {
            return Character.toLowerCase(this.getTypeChar());
        }
    }

    public String getOrientation() {
        return this.orientation;
    }


    public void expandBy(Residue r) {
        this.residues.add(r);
    }

    public boolean containedInPDBNumberRange(int pdbResidueNumberStart, int pdbResidueNumberEnd) {
        try {
            return this.firstPDB() >= pdbResidueNumberStart && this.lastPDB() <= pdbResidueNumberEnd;
        } catch (NoSuchElementException n) {
            return false;
        }
    }
    
    public boolean overlapsPDBNumberRange(int pdbResidueNumberStart, int pdbResidueNumberEnd) {
        try {
            int first = this.firstPDB();
            int last = this.lastPDB();
            return (first >= pdbResidueNumberStart && first <= pdbResidueNumberEnd) ||
                    (last >= pdbResidueNumberStart && last <= pdbResidueNumberEnd);
        } catch (NoSuchElementException n) {
            return false;
        }
    }

    public boolean containsPDBNumber(int pdbResidueNumber) {
       for (Residue r : this.residues) {
            if (r.getPDBNumber() == pdbResidueNumber) {
                return true;
            }
        }
        return false;
    }

    public boolean containsAbsoluteNumber(int absoluteResidueNumber) {
        try {
            Residue first = this.firstResidue();
            Residue last = this.lastResidue();
            return first.getAbsoluteNumber() <= absoluteResidueNumber && last.getAbsoluteNumber() >= absoluteResidueNumber;
        } catch (NoSuchElementException n) {
            return false;
        }
    }

    public boolean contains(Residue r) {
        return this.containsAbsoluteNumber(r.getAbsoluteNumber());
    }

    public Axis getAxis() {
        if (this.axis == null) {
            this.calculateAxis();
        }
        return this.axis;
    }

    public void calculateAxis() {
        //this.axis = Geometer.leastSquareAxis(this.getCAlphaCoordinates());
        // if this segment is only one or no residues, make a zero axis
        if (this.residues.size() < 2) {
            this.axis = new Axis();
            return;
        }

        // otherwise, diff the centers of the first and last residues
        try {
            Point3d start = this.firstResidue().getCenter();
            Point3d end = this.lastResidue().getCenter();
            this.axis = new Axis(start, end);
            this.axis.setStart(start);
            this.axis.setEnd(end);
        } catch (NoSuchElementException n) {
            System.out.println("NoSuchElementException for " + this);
        }
        //System.out.println("setting axis of " + this + " to : " + this.axis);
    }

    // the axis we pass into the function is considered to be "UP"
    public void determineOrientation(Axis axis) {
        double angle = this.getAxis().angle(axis);
        Logger.getLogger("translation.FoldAnalyser").info("Angle of " + this + " with axis " + axis + " is " + angle);
        if (angle > 90) {
            this.orientation = "DOWN";
        } else {
            this.orientation = "UP";
        }
    }

    public void setOrientation(String orientation) {
        this.orientation = orientation;
    }

    public char getRelativeOrientation(BackboneSegment other) {
        return this.getRelativeOrientation(other.getAxis());
    }

    public char getRelativeOrientation(Axis other) {
        double angle = this.getAxis().angle(other);
        //System.out.println("angle between " + this + " and " + other + " =  " + angle);
        if (angle > 90) {
            return 'A';
        } else {
            return 'P';
        }
    }
}
