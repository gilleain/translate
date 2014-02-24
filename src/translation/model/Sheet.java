package translation.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.TreeMap;

import java.util.logging.Logger;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import translation.Axis;
import translation.Geometer;

public class Sheet {

    private int number;
    //private ArrayList strands;

    // The strand map has as the keys all the strands in the map;
    // the values are the other strands in the sheet that the key is attached to.
    // In theory, a nice pure sheet would only have one value per key. In theory...
    private TreeMap strandMap;
    private Axis axis;

    public Sheet(int number) {
        this.number = number;
        //this.strands = new ArrayList();
        this.strandMap = new TreeMap();
        this.axis = null;
    }

    public Sheet(int number, BackboneSegment first, BackboneSegment second) {
        this(number);
        this.addPair(first, second);
    }

    public void addPair(BackboneSegment first, BackboneSegment second) {
        if (first.compareTo(second) < 0) {
            this.map(first, second);
        } else {
            this.map(second, first);
        }
    }

    public void map(BackboneSegment keyStrand, BackboneSegment partner) {
        ArrayList values;
        if (this.strandMap.containsKey(keyStrand)) {
            values = (ArrayList) this.strandMap.get(keyStrand);
        } else {
            values = new ArrayList();
            this.strandMap.put(keyStrand, values);
        }
        values.add(partner);
    }

    public ArrayList getPartners(BackboneSegment key) {
        return (ArrayList) this.strandMap.get(key);
    }

    public Iterator getPartnerIterator(BackboneSegment key) {
        return this.getPartners(key).iterator();
    }

    public int getNumber() {
        return this.number;
    }

    public void setAxis(Axis axis) {
        this.axis = axis;
    }

    public void setAxis(Vector3d axisVector) {
        this.setAxis(new Axis(this.calculateCentroid(), axisVector));
    }

    public Axis getAxis() {
        return this.axis;
    } 

    public Point3d calculateCentroid() {
        ArrayList centers = new ArrayList();
        Iterator iterator = this.strandMap.keySet().iterator();
        while (iterator.hasNext()) {
            BackboneSegment strand = (BackboneSegment) iterator.next();
            centers.add(strand.getAxis().getCentroid());
        }
        return Geometer.averagePoints((Collection)centers);
    }

    public int size() {
        int size = 0;
        Iterator keyIterator = this.iterator();
        while (keyIterator.hasNext()) {
            size++;
            size += ((ArrayList) this.strandMap.get(keyIterator.next())).size();
        }
        return size;
    }

    public void extend(Sheet other) {
        Iterator keyIterator = other.iterator();
        while (keyIterator.hasNext()) {
            BackboneSegment key = (BackboneSegment) keyIterator.next();
            ArrayList otherValues = other.getPartners(key);

            if (this.strandMap.containsKey(key)) {
                ArrayList thisValues = this.getPartners(key);
                thisValues.addAll(otherValues);
            } else {
                this.strandMap.put(key, otherValues);
            }
        }
    }

    public Iterator iterator() {
        return this.strandMap.keySet().iterator();
    }

    public Iterator chainOrderIterator() {
        ArrayList chainOrder = new ArrayList();
        chainOrder.addAll(this.strandMap.keySet());
        Iterator iterator = this.iterator();
        while (iterator.hasNext()) {
            ArrayList partners = (ArrayList) this.strandMap.get(iterator.next());
            for (int i = 0; i < partners.size(); i++) {
                BackboneSegment partner = (BackboneSegment) partners.get(i);
                if (!chainOrder.contains(partner)) {
                    chainOrder.add(partner);
                }
            }
        }
        Collections.sort(chainOrder);
        return chainOrder.iterator();
    }

    public boolean contains(BackboneSegment strand) {
        if (this.strandMap.containsKey(strand)) {
            return true;
        } else {
            Iterator keyIterator = this.iterator();
            while (keyIterator.hasNext()) {
                if (((ArrayList) this.strandMap.get(keyIterator.next())).contains(strand)) {
                    return true;
                }
            }
        }
        return false;
    }

    public ArrayList getSheetPaths() {
        ArrayList paths = new ArrayList();
        Iterator iterator = this.iterator();
        while (iterator.hasNext()) {
            BackboneSegment key = (BackboneSegment) iterator.next();
            paths.add(this.traverseSheetPath(key, new ArrayList()));
        }
        return paths;
    }

    public ArrayList traverseSheetPath(BackboneSegment currentStrand, ArrayList path) {
        ArrayList partners = (ArrayList) this.strandMap.get(currentStrand);
        for (int i = 0; i < partners.size(); i++) {
            BackboneSegment partner = (BackboneSegment) partners.get(i);
            path.add(partner);
            if (this.strandMap.containsKey(partner)) {
                this.traverseSheetPath(partner, path);            
            } else {
                return path;
            }
        } 
        return path;
    }

    public void assignOrientationsToStrands() {

        // while we're at it, we might as well calculate the sheet axis
        Vector3d sheetVector = new Vector3d();

        // reset the iterator
        Iterator iterator = this.iterator();
        while (iterator.hasNext()) {
            BackboneSegment strand = (BackboneSegment) iterator.next();
            if (strand.getOrientation().equals("None")) {
                strand.setOrientation("UP");
            }

            Iterator partnerIterator = this.getPartnerIterator(strand);
            while (partnerIterator.hasNext()) {
                BackboneSegment partner = (BackboneSegment) partnerIterator.next();
                this.assignOrientations(strand, partner);

                // add to the average vector, or subtract if DOWN
                if (partner.getOrientation().equals("UP")) {
                    sheetVector.add(partner.getAxis().getAxisVector());
                } else {
                    sheetVector.sub(partner.getAxis().getAxisVector());
                }
                Logger.getLogger("translation.FoldAnalyser").info("orientation " + strand + " -> " + partner);
            }
        }

        sheetVector.normalize();
        this.setAxis(sheetVector);
    }

    private void assignOrientations(BackboneSegment strand, BackboneSegment partner) {
        char relativeOrientation = strand.getRelativeOrientation(partner);
        if (relativeOrientation == 'P') {
            String strandOrientation = strand.getOrientation();
            if (strandOrientation.equals("None")) {
                String partnerOrientation = partner.getOrientation();
                if (partnerOrientation.equals("None")) {
                    Logger.getLogger("translation.FoldAnalyser").info("No orientation known for " + strand + " and " + partner);
                } else {
                    Logger.getLogger("translation.FoldAnalyser").info("Assigning orientation : " + partnerOrientation + " to " + strand);
                    strand.setOrientation(partnerOrientation);
                }
            } else {
                Logger.getLogger("translation.FoldAnalyser").info("Assigning orientation : " + strandOrientation + " to " + partner);
                partner.setOrientation(strandOrientation);
            }
        } else {
            String strandOrientation = strand.getOrientation();
            if (strandOrientation.equals("None")) {
                String partnerOrientation = partner.getOrientation();
                if (partnerOrientation.equals("None")) {
                    Logger.getLogger("translation.FoldAnalyser").info("No orientation known for " + strand + " and " + partner);
                } else {
                    Logger.getLogger("translation.FoldAnalyser").info("Assigning orientation : " + partnerOrientation + " to " + strand);
                    strand.setOrientation(partnerOrientation);
                }
            } else if (strandOrientation.equals("UP")) {
                Logger.getLogger("translation.FoldAnalyser").info("Assigning orientation : UP to " + partner);
                partner.setOrientation("DOWN");
            } else if (strandOrientation.equals("DOWN")) {
                Logger.getLogger("translation.FoldAnalyser").info("Assigning orientation : DOWN to " + partner);
                partner.setOrientation("UP");
            }
        }
    }


    public ArrayList toTopsEdges(Domain domain) {
        ArrayList edges = new ArrayList();
        Iterator keyIterator = this.strandMap.keySet().iterator();

        while (keyIterator.hasNext()) {
            BackboneSegment strand = (BackboneSegment) keyIterator.next();
            if (!domain.contains(strand)) {
                continue;
            }

            Iterator partnerIterator = this.getPartnerIterator(strand);
            while (partnerIterator.hasNext()) {
                BackboneSegment partner = (BackboneSegment) partnerIterator.next();
                if (!domain.contains(partner)) {
                    continue;
                }

                char relativeOrientation = strand.getRelativeOrientation(partner); 
                Edge edge;
                if (strand.compareTo(partner) < 0) {
                    edge = new Edge(strand, partner, relativeOrientation);
                } else {
                    edge = new Edge(partner, strand, relativeOrientation);
                }
                //System.err.println("Made edge : " + edge);

                // since the mapping is symmetric, we have to discard half the edges we make!
                if (!edges.contains(edge)) {
                    edges.add(edge);
                }
            }
        } 

        return edges;
    }

    public String toString() {
        StringBuffer returnValue = new StringBuffer();
        returnValue.append("Sheet (" + this.number + ") [");

        //Iterator iterator = this.strands.iterator();
        Iterator iterator = this.strandMap.keySet().iterator();

        while (iterator.hasNext()) {
            BackboneSegment strand = (BackboneSegment) iterator.next();
            returnValue.append(strand);

            Iterator partnerIterator = this.getPartnerIterator(strand);
            while (partnerIterator.hasNext()) {
                BackboneSegment partner = (BackboneSegment) partnerIterator.next();
                returnValue.append(" -> ").append(partner); 
            }
            returnValue.append("\n");
        }
        returnValue.append("]");
        return returnValue.toString();
    }

}
