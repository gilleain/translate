package translation.model;


public class UnstructuredSegment extends BackboneSegment {

    public UnstructuredSegment() {
        super();
    }

    public UnstructuredSegment(Residue r) {
        super(r);
    }

    public char getTypeChar() {
        return 'U';
    }

    public String toString() {
        return "Unstructured from " + this.firstResidue().getPDBNumber() + " to " + this.lastResidue().getPDBNumber();
    }

    public String toFullString() {
        return this.toString();
    }
}
