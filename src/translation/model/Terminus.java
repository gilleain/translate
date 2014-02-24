package translation.model;


public class Terminus extends BackboneSegment {
    private String label;
    private char typeChar;

    public Terminus(String label, char typeChar) {
        this.label = label;
        this.typeChar = typeChar;
    }

    public char getTypeChar() {
        return this.typeChar;
    }

    public String getOrientation() {
        return "UP";
    }

    public String toString() {
        return this.label;
    }    

    public String toFullString() {
        return this.toString();
    }

}
