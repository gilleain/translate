package translation;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

public class Axis {
    private Point3d centroid;
    private Vector3d axisVector;
    private Point3d start;
    private Point3d end;

    public static double ANGLE_DELTA = 5.0;

    public Axis() {
        this.centroid = new Point3d();
        this.axisVector = new Vector3d();
        this.start = this.centroid;
        this.end = this.centroid;
    }

    public Axis(Point3d start, Point3d end) {
        this.centroid = new Point3d();
        this.centroid.interpolate(start, end, 0.5);
        this.axisVector = new Vector3d();
        this.axisVector.sub(end, start);
        this.axisVector.normalize();
    }

    public Axis(Point3d centroid, Vector3d axisVector) {
        this.centroid = centroid;
        this.axisVector = axisVector;
    }

    public Point3d getCentroid() {
        return this.centroid;
    }

    public Vector3d getAxisVector() {
        return this.axisVector;
    }

    public void setStart(Point3d start) {
        this.start = Geometer.scalePoint(start, this.axisVector, this.centroid);
    }

    public Point3d getStart() {
        return this.start;
    }

    public void setEnd(Point3d end) {
        this.end = Geometer.scalePoint(end, this.axisVector, this.centroid);
    }

    public Point3d getEnd() {
        return this.end;
    }

    public double angle(Vector3d otherVector) {
        return Math.toDegrees(this.axisVector.angle(otherVector));
    }

    public double angle(Axis other) {
        return this.angle(other.axisVector);
    }

    public boolean approximatelyLinearTo(Axis other) {
        //return (this.angle(other) - Axis.ANGLE_DELTA) < 0;
        double angleDiff = this.angle(other) - Axis.ANGLE_DELTA;
        //System.err.println("angle diff = " + angleDiff);
        return angleDiff < 0;
    }

    public String toString() {
        if (this.start == null) {
            return "[], []";
        } else {
            //String startS = String.format("(%.2f, %.2f, %.2f)", this.start.x, this.start.y, this.start.z);
            //String endS = String.format("(%.2f, %.2f, %.2f)", this.end.x, this.end.y, this.end.z);

            //String centroidS = String.format("(%.2f, %.2f, %.2f)", this.centroid.x, this.centroid.y, this.centroid.z);
            //String axisS = String.format("(%.2f, %.2f, %.2f)", this.axisVector.x, this.axisVector.y, this.axisVector.z);

            //return String.format("[%s, %s]", centroidS, axisS);
            //return String.format("[%s, %s], [%s, %s]", startS, endS, centroidS, axisS);
            
            return "";
        }
    }

}
