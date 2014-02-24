package translation;

import java.io.IOException;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.HashMap;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import translation.model.BackboneSegment;
import translation.model.Chain;
import translation.model.HBond;
import translation.model.Helix;
import translation.model.Protein;
import translation.model.Residue;
import translation.model.Sheet;
import translation.model.Strand;
import translation.model.Terminus;
import translation.model.UnstructuredSegment;

public class StructureFinder {
    private Protein protein;

    public StructureFinder(String filename) throws IOException {
        this.protein = PDBReader.read(filename);
    }

    public Protein getProtein() {
        this.protein.findStructure(this);
        return this.protein;
    }

    // The 'interface' to the function of this class
    // I suppose that rigorously, this should be the only public method...
    public void findStructure(Chain chain) {
        chain.calculateTorsions();
        this.calculateHBondPartners(chain);
        //this.assignTorsionsAndHBondsToTypes(chain);
        //this.convertTorsionsToRepetitiveStructure(chain);
        //this.cleanStructure(chain);
        //this.cleanStructure(chain);
        this.buildSSES(chain);
        this.findSheets(chain);
        this.assignOrientationsDependingOnArchitecture(chain);
        this.determineChiralities(chain);
    }

    // go through the fixed structures, finding helices connecting parallel strands not more than 5 sses away
    public void determineChiralities(Chain chain) {
        if (chain.numberOfSheets() > 0) {
            Iterator sheets = chain.sheetIterator();
            while (sheets.hasNext()) {
                Sheet sheet = (Sheet) sheets.next();
                Iterator sheetIterator = sheet.iterator();
                BackboneSegment strand = (BackboneSegment) sheetIterator.next();

                while (sheetIterator.hasNext()) {
                    BackboneSegment partner = (BackboneSegment) sheetIterator.next();
                    boolean isCorrectOrder = strand.compareTo(partner) < 0;
                    Vector3d upVector = strand.getAxis().getAxisVector();

                    // if the strand and its partner are less than x segments away, and parallel
                    int distance;
                    if (isCorrectOrder) {
                        distance = chain.numberOfSegmentsInBetween(strand, partner);
                    } else {
                        distance = chain.numberOfSegmentsInBetween(partner, strand);
                    }

                    if ((distance < 11) && (strand.getOrientation().equals(partner.getOrientation()))) {
                        Point3d strandCentroid = strand.getAxis().getCentroid();
                        Point3d partnerCentroid = partner.getAxis().getCentroid();

                        // get the sses we will use for the chirality calculation
                        ListIterator inBetweeners;
                        if (isCorrectOrder) {
                            inBetweeners = chain.backboneSegmentListIterator(strand, partner);
                        } else {
                            inBetweeners = chain.backboneSegmentListIterator(partner, strand);
                        }

                        // find the average center of these sses
                        ArrayList centroids = new ArrayList();
                        while (inBetweeners.hasNext()) {
                            BackboneSegment segment = (BackboneSegment) inBetweeners.next();
                            centroids.add(segment.getAxis().getCentroid()); 
                        }
                        Point3d averageCentroid = Geometer.averagePoints((Collection) centroids);
                        //System.out.println("InBetweener centroid = " + averageCentroid);

                        //finally, do the calculation
                        char chirality;
                        if (isCorrectOrder) {
                            chirality = Geometer.chirality(strandCentroid, averageCentroid, partnerCentroid, upVector);
                        } else {
                            chirality = Geometer.chirality(partnerCentroid, averageCentroid, strandCentroid, upVector);
                        }
                    
                        if (isCorrectOrder) {
                            chain.addChirality(strand, partner, chirality);
                        } else {
                            chain.addChirality(partner, strand, chirality);
                        }
                    }
                    strand = partner;
                }
            }
        }
    }

    //basically, using the least-square plane to find the best chain axis is a rubbish way to do things!
    //unfortunately, the only alternative is to have some rather ad-hoc rules based on architecture
    public void assignOrientationsDependingOnArchitecture(Chain chain) {
        //for a single sheet, get the orientation of the sheet, and use that
        int sheetCount = chain.numberOfSheets();
        if (sheetCount == 1) {
            // get the sheet, and assign the orientations internally, based on relative orientations
            Iterator sheets = chain.sheetIterator();
            Sheet sheet = (Sheet) sheets.next();

            //System.out.println("Assigning using a single sheet");
            sheet.assignOrientationsToStrands();

            Axis sheetAxis = sheet.getAxis();
            //System.out.println("Sheet axis = " + sheetAxis.getCentroid() + ", " + sheetAxis.getAxisVector());

            // assign the helices to orientations dependant on the sheet 
            ListIterator segments = chain.backboneSegmentListIterator();
            while (segments.hasNext()) {
                BackboneSegment segment = (BackboneSegment) segments.next();
                if ((segment instanceof Strand) || (segment instanceof Terminus)) {
                    continue;
                } else {
                    segment.determineOrientation(sheetAxis); 
                }
            }
        } else if (sheetCount > 1) {
            Iterator sheets = chain.sheetIterator();
            while (sheets.hasNext()) {
                Sheet sheet = (Sheet) sheets.next();
                sheet.assignOrientationsToStrands();
            }
        }
    }

    public void assignTorsionsAndHBondsToTypes(Chain chain) {
        //run through once, using hbonds and torsions to attempt individual residue assignments
        Iterator residueIterator = chain.residueIterator();
        HashMap residueAssignments = new HashMap();

        //temporary buffers to store char assignents
        StringBuffer torsionBuffer = new StringBuffer();
        StringBuffer hbondBuffer = new StringBuffer();
        StringBuffer assignBuffer = new StringBuffer();

        StringBuffer tmpTorsion = new StringBuffer();
        StringBuffer tmpHBond = new StringBuffer();

        int spaceCounter = 1;
        while (residueIterator.hasNext()) {
            Residue residue = (Residue) residueIterator.next();

            char typeFromTorsion = this.determineTypeFromTorsion(residue);
            torsionBuffer.append(typeFromTorsion);

            tmpTorsion.append(typeFromTorsion);
            if (spaceCounter % 10 == 0) { tmpTorsion.append(' '); }

            char typeFromHBonds  = this.determineTypeFromHBonds(residue);
            hbondBuffer.append(typeFromHBonds);

            tmpHBond.append(typeFromHBonds);
            if (spaceCounter % 10 == 0) { tmpHBond.append(' '); }
            
            char assignment = this.assign(typeFromTorsion, typeFromHBonds);
            assignBuffer.append(assignment);
            if (spaceCounter % 10 == 0) { assignBuffer.append(' '); }

            residueAssignments.put(residue, new Character(assignment));

            spaceCounter++;
        }
        //System.out.println("torsions : " + tmpTorsion);
        //System.out.println("hbonds   : " + tmpHBond);
        //System.out.println("final    : " + assignBuffer);

        //now, run through the assignments, merging into sses
        int sseStart = 1;
        int sseEnd = -1;
        char currentSSEType = 'U';
        ArrayList backboneSegments = new ArrayList();
        for (int i = 0; i < torsionBuffer.length(); i++) {
            char sseChar = torsionBuffer.charAt(i);
            // same type : extend the end of current
            if (sseChar == currentSSEType) {
                sseEnd = i + 1;
            // change of type : finish previous and start new
            } else {
                //System.out.println("SSE : " + currentSSEType + " from " + sseStart + " to " + sseEnd);
                backboneSegments.add(this.createSegment(sseStart, sseEnd, currentSSEType, chain));
                sseStart = i + 1;
                sseEnd = i + 1;
            } 
            currentSSEType = sseChar;
        }
        //System.out.println("SSE : " + currentSSEType + " from " + sseStart + " to " + sseEnd);
//        Iterator itr = backboneSegments.iterator();
//        while (itr.hasNext()) {
//            BackboneSegment segment = (BackboneSegment) itr.next();
//            System.out.println(segment);
//        }
        this.trimByHBonds(backboneSegments, hbondBuffer); 
//        itr = backboneSegments.iterator();
//        while (itr.hasNext()) {
//            BackboneSegment segment = (BackboneSegment) itr.next();
//            System.out.println(segment);
//        }
    }
    
    public void trimByHBonds(ArrayList backboneSegments, StringBuffer hbondBuffer) {
        ListIterator itr = backboneSegments.listIterator();
        BackboneSegment previousSegment = null;

        while (itr.hasNext()) {
            // get the current and the next segments in the list
            BackboneSegment currentSegment = (BackboneSegment) itr.next();
            BackboneSegment nextSegment = null;
            if (itr.hasNext()) {
                nextSegment = (BackboneSegment) itr.next();
                itr.previous();
            }

            // helices
            if ((currentSegment instanceof Helix) && (currentSegment.length() < 4)) {
                if (previousSegment != null && previousSegment instanceof UnstructuredSegment) {
                    if (nextSegment != null && nextSegment instanceof UnstructuredSegment) {
                        this.mergeThreeSegments(itr, previousSegment, currentSegment, nextSegment);            
                    } else {
                        previousSegment.mergeWith(currentSegment);
                        itr.remove();
                        currentSegment = previousSegment;
                    }
                } else {
                    if (nextSegment != null && nextSegment instanceof UnstructuredSegment) {
                        nextSegment.mergeWith(currentSegment);
                        itr.previous();
                        itr.remove();
                        currentSegment = nextSegment;
                    } else {
                        // a short helix surrounded by strand? => merge?
                        if (nextSegment != null) {
                            this.mergeThreeSegments(itr, previousSegment, currentSegment, nextSegment);            
                        }
                    }
                }
            }

            // strands
    
            previousSegment = currentSegment;
            //System.out.println("Setting previous segment to : " + previousSegment);
        }
    }

    public BackboneSegment createSegment(int start, int end, int sseType, Chain chain) {
        // make a new segment of the appropriate type
        BackboneSegment newSegment = null;
        if (sseType == 'U') {
            newSegment = new UnstructuredSegment();
        } else if (sseType == 'E') {
            newSegment = new Strand();
        } else if (sseType == 'H') {
            newSegment = new Helix();
        }

        // fill the new segment with residues
        for (int i = start; i <= end; i++) {
            Residue residue = chain.getResidueByAbsoluteNumbering(i);
            newSegment.expandBy(residue);
        }

        // return this new, filled segment
        return newSegment;
    }

    public char assign(char typeFromTorsion, char typeFromHBonds) {
        //firstly, if we agree, then all is fine and dandy
        if (typeFromTorsion == typeFromHBonds) {
            return typeFromTorsion;
        } else {
            //torsions may fall outside the narrow bounds defined for them
            if (typeFromTorsion == 'U') {
                return typeFromHBonds;
            } else {
                //it is possible for a residue in an SSE not to have any bonds
                if (typeFromHBonds == 'U') {
                    return typeFromTorsion;
                } else {
                    //torsion == 'H' and hbonds == 'E'
                    if (typeFromTorsion == 'H') {
                        return 'E';
                    //torsion == 'E' and hbonds == 'H' - can this happen?
                    } else {
                        return 'H';
                    }
                }
            }
        }
    }

    public char determineTypeFromTorsion(Residue residue) {
        if (Strand.torsionsMatch(residue)) {
            return 'E';
        } else if (Helix.torsionsMatch(residue)) {
            return 'H';
        } else {
            return 'U';
        }
    }

    public char determineTypeFromHBonds(Residue residue) {
        if (Helix.hbondsMatch(residue)) {
            return 'H';
        } else if (Strand.hbondsMatch(residue)) {
            return 'E';
        } else {
            return 'U';
        }
    }

    public void findSheets(Chain chain) {
        //for each strand, examine all strands in front (so we don't do the same comparison twice)
        //each examination is a simple centroid-centroid distance
        ListIterator firstSegments = chain.backboneSegmentListIterator();
    
        while (firstSegments.hasNext()) {
            // get the first segment, reject if not a strand
            BackboneSegment firstSegment = (BackboneSegment) firstSegments.next();
            if (!(firstSegment instanceof Strand)) {
                continue;
            }
           
            // get the segments after the current one
            ListIterator secondSegments = chain.backboneSegmentListIterator(firstSegment);
            while (secondSegments.hasNext()) {
                BackboneSegment secondSegment = (BackboneSegment) secondSegments.next();
                if ((secondSegment != firstSegment) && (secondSegment instanceof Strand)) {
                    //make a crude distance check
                    if (this.closeApproach(firstSegment, secondSegment)) {
                        //if this passes, make a finer bonding check
                        if (this.bonded(firstSegment, secondSegment)) {
                            //this.addStrandPair(firstSegment, secondSegment, chain);
                        }
                    }
                }
            }
        }
    }

    /*
    public void addStrandPair(BackboneSegment first, BackboneSegment second, Chain chain) {
        Sheet firstSheet = chain.getSheetContaining(first);
        Sheet secondSheet = chain.getSheetContaining(second);
        
        if (firstSheet == null && secondSheet == null) {
            chain.createSheet(first, second);
            //System.out.println("Adding " + first + " and " + second + " to new sheet");
        } else {
            if (firstSheet == null) {
                //System.out.println("Adding " + first + " and " + second + " to " + secondSheet);
                secondSheet.insert(second, first);
            } else if (secondSheet == null) {
                //System.out.println("Adding " + first + " and " + second + " to " + firstSheet);
                firstSheet.insert(first, second);
            } else {
                // one possibility is that both sheets are the same, and this edge closes a barrel
                if (firstSheet == secondSheet) {
                    firstSheet.closeBarrel(first, second);
                } else {
                    // otherwise, we have to join the sheets
                    this.joinSheets(first, second, firstSheet, secondSheet, chain);
                }
            }
        }
    }

    public void joinSheets(BackboneSegment first, BackboneSegment second, Sheet firstSheet, Sheet secondSheet, Chain chain) {
        // check that the strands are not somehow in the middle of the sheet (bifurcated sheets)
        if (firstSheet.strandInMiddle(first) || secondSheet.strandInMiddle(second)) {
            System.err.println("bifurcation in : Sheet (" + firstSheet + ") and Strand (" + first + ") or Sheet (" + secondSheet + ") and Strand (" + second + " )");
            return;
        }

        int indexOfFirst = firstSheet.indexOf(first);
        int indexOfSecond = secondSheet.indexOf(second);

        if (indexOfSecond == 0) {
            if (indexOfFirst == 0) {
                //System.out.println("Reversing : " + firstSheet + " and adding : " + secondSheet + " onto the end");
                firstSheet.reverse();
                firstSheet.extend(secondSheet);
            } else {
                //System.out.println("Adding : " + secondSheet + " onto the end of " + firstSheet);
                firstSheet.extend(secondSheet);
            }
            chain.removeSheet(secondSheet);
        } else {
            if (indexOfFirst == 0) {
                //System.out.println("Adding : " + firstSheet + " onto the end of " + secondSheet);
                secondSheet.extend(firstSheet);
                chain.removeSheet(firstSheet);
            } else {
                //System.out.println("Reversing : " + secondSheet + " and adding it onto the end of " + firstSheet);
                secondSheet.reverse();
                firstSheet.extend(secondSheet);
                chain.removeSheet(secondSheet);
            }
        }
    }
    */

    public boolean closeApproach(BackboneSegment a, BackboneSegment b) {
        Vector3d distanceVector = new Vector3d();
        distanceVector.sub(a.getAxis().getCentroid(), b.getAxis().getCentroid());
        double length = distanceVector.length();
        //System.out.println("Distance between " + a + " and " + b + " = " + Math.rint(length));
        return length < 10.0;
    }

    public boolean bonded(BackboneSegment strand, BackboneSegment otherStrand) {
        //basically, run through the residues, checking the list of hbonds to find residues that might be in the other strand
        int numberOfHBonds = 0;
        Iterator strandResidues = strand.residueIterator();
        while (strandResidues.hasNext()) {
            Residue nextResidue = (Residue) strandResidues.next();
            if (otherStrand.bondedTo(nextResidue)) {
                numberOfHBonds++;
            }
        }

        //all we want to know is, are there enough hbonds for these strands to qualify as bonded
        //WARNING! this is number of RESIDUES not number of HBONDS
        if (numberOfHBonds > 1) { //?is two enough?
            return true;
        } else {
            //System.out.println(strand + " not bonded to " + otherStrand);
            return false;
        }
    }


    //merge RepetitiveStructure separated by only a single unstructured residue
    //also, delete single-residue RepetitiveStructure surrounded by UnstructuredRegions
    public void cleanStructure(Chain chain) {
        ListIterator backboneSegmentIterator = chain.backboneSegmentListIterator();
        BackboneSegment currentSegment = null;

        while (backboneSegmentIterator.hasNext()) {
            BackboneSegment previousSegment = currentSegment;
            currentSegment = (BackboneSegment) backboneSegmentIterator.next();

            //Don't bother about the terminii
            if (currentSegment instanceof Terminus) {
                continue;
            }

            //Consider the unstructured parts in between for deletion
            if (currentSegment instanceof UnstructuredSegment) {
                //only merge if the unstructured segment is short
                if (currentSegment.length() > 2) {
                    continue;
                }
                //check the segments 'fore and 'aft - are they continuous?
                if (backboneSegmentIterator.hasNext()) {
                    BackboneSegment nextSegment = (BackboneSegment) backboneSegmentIterator.next();
                    //System.err.println("Checking : " + previousSegment + " and " + currentSegment + " and " + nextSegment);
                    if (previousSegment.continuousWith(nextSegment)) {
                        this.mergeThreeSegments(backboneSegmentIterator, previousSegment, currentSegment, nextSegment);
                        continue;
                    } else {
                        backboneSegmentIterator.previous();
                    }
                }
            //Or, if it is repetitive, and surrounded by unstructured, and short, delete it.
            } else {
                //only delete if the Repetitive segment is short
                if (currentSegment.length() > 2) {
                    continue;
                }

                //check the segments 'fore and 'aft - are they unstructured?
                if (backboneSegmentIterator.hasNext()) {
                    BackboneSegment nextSegment = (BackboneSegment) backboneSegmentIterator.next();
                    //System.err.println("Checking : " + previousSegment + " and " + currentSegment + " and " + nextSegment);
                    if (previousSegment instanceof UnstructuredSegment && nextSegment instanceof UnstructuredSegment) {
                        this.mergeThreeSegments(backboneSegmentIterator, previousSegment, currentSegment, nextSegment);
                        continue;
                    //} else if (previousSegment.getClass() == ){
                    } else {
                        backboneSegmentIterator.previous();
                    }
                }
                
            }
        }
    }

    public void mergeThreeSegments(ListIterator backboneSegmentIterator, BackboneSegment previousSegment, BackboneSegment currentSegment, BackboneSegment nextSegment) {
        //System.err.println("Merging : " + previousSegment + " and " + currentSegment + " and " + nextSegment);
        //merge the previous segment with the current segment and the next segment
        previousSegment.mergeWith(currentSegment);
        previousSegment.mergeWith(nextSegment);
        //now, delete the current and next segments
        backboneSegmentIterator.previous();
        backboneSegmentIterator.remove();
        backboneSegmentIterator.previous();
        backboneSegmentIterator.remove();
    }

    public BackboneSegment fitNextResidue(Residue residue, BackboneSegment currentBackboneSegment) {
        if (currentBackboneSegment instanceof Strand) {
            if (Strand.torsionsMatch(residue)) {
                currentBackboneSegment.expandBy(residue);
                return currentBackboneSegment;
            } else if (Helix.torsionsMatch(residue)) {
                return new Helix(residue);
            } else {
                return new UnstructuredSegment(residue);
            }
        } else if (currentBackboneSegment instanceof Helix) {
            if (Helix.torsionsMatch(residue)) {
                currentBackboneSegment.expandBy(residue);
                return currentBackboneSegment;
            } else if (Strand.torsionsMatch(residue)) {
                return new Strand(residue);
            } else {
                return new UnstructuredSegment(residue);
            }
        } else {
            if (Strand.torsionsMatch(residue)) {
                return new Strand(residue);
            } else if (Helix.torsionsMatch(residue)) {
                return new Helix(residue);
            } else {
                currentBackboneSegment.expandBy(residue);
                return currentBackboneSegment;
            }
        }
    }

    public void convertTorsionsToRepetitiveStructure(Chain chain) {
        Iterator residueIterator = chain.residueIterator();
        Terminus nterminus = new Terminus("N Terminus", 'N');
        chain.addBackboneSegment(nterminus);
        BackboneSegment currentBackboneSegment = new UnstructuredSegment();

        while (residueIterator.hasNext()) {
            Residue residue = (Residue) residueIterator.next();
            BackboneSegment nextBackboneSegment = this.fitNextResidue(residue, currentBackboneSegment);
            if (nextBackboneSegment != currentBackboneSegment) {
                chain.addBackboneSegment(currentBackboneSegment);   //store the previous segment
            }
            currentBackboneSegment = nextBackboneSegment;
        }
        chain.addBackboneSegment(currentBackboneSegment);       //store the final segment
        chain.addBackboneSegment(new Terminus("C Terminus", 'C'));    //and add a C Terminus for luck
    }

    public void calculateHBondPartners(Chain c) {
        Iterator itr = c.residueIterator();
        while (itr.hasNext()) {
            Residue first = (Residue) itr.next();
            this.searchForwards(first, c);
        }
    }

    public void searchForwards(Residue first, Chain c) {

        // check that it makes sense to be doing this
        if (!first.isStandardAminoAcid()) {
            return;
        }

        int position = first.getAbsoluteNumber();
        int nextPosition = position + 1;
        Point3d firstN = first.getCoordinates("N");
        Point3d firstH = first.getCoordinates("H");
        Point3d firstO = first.getCoordinates("O");
        Point3d firstC = first.getCoordinates("C");

        // FIXME : unfortunately, this misses out on PRO residues (also below)
        if (firstN == null || firstH == null || firstO == null || firstC == null) {
            return;
        }

        // allow for chain breaks, or return if we have reached the end
        if (!c.hasResidueByAbsoluteNumbering(nextPosition)) {
            Residue second = c.getNextResidue(position);
            if (second == null) {   //probably reached the end of the chain!
                return ;
            }
            nextPosition = second.getAbsoluteNumber();
        }

        // now, compare the first residue to the residues further on in the chain
        Iterator itr = c.residueIterator(nextPosition);
        while (itr.hasNext()) {
            int secondPosition = ((Residue) itr.next()).getAbsoluteNumber();
            if (secondPosition < (position + 3)) {
                continue;
            }

            Residue second;
            try {
                second = c.getResidueByAbsoluteNumbering(secondPosition);
                if (!second.isStandardAminoAcid()) {
                    continue;
                }
            } catch (IndexOutOfBoundsException i) {
                break;
            }

            try {
                Point3d secondN = second.getCoordinates("N");
                Point3d secondH = second.getCoordinates("H");
                Point3d secondO = second.getCoordinates("O");
                Point3d secondC = second.getCoordinates("C");

                // FIXME : PRO residues...
                if (secondN == null || secondH == null || secondO == null || secondC == null) {
                    continue;
                }

                // bonds from first N-H to second C=O
                HBond firstSecondBond = this.createBond(first, second, firstN, firstH, secondO, secondC);
                if (firstSecondBond != null) {
                    c.addHBond(firstSecondBond);
                    first.addHBond(firstSecondBond);
                    second.addHBond(firstSecondBond);
                }

                // bonds from second N-H to first C=O
                HBond secondFirstBond = this.createBond(second, first, secondN, secondH, firstO, firstC);
                if (secondFirstBond != null) {
                    c.addHBond(secondFirstBond);
                    first.addHBond(secondFirstBond);
                    second.addHBond(secondFirstBond);
                }
            } catch (NullPointerException npe) {
                System.err.println(npe + " for " + first + " and " + second);
                continue;
            }
        }
        // finally, as an optimisation measure, use this loop over the residues to analyse the bonding environment
        this.analyzeHBonds(first);
    }

    public HBond createBond(Residue donor, Residue acceptor, Point3d nitrogen, Point3d hydrogen, Point3d oxygen, Point3d carbon) {
        double distance = hydrogen.distance(oxygen);
        double nhoAngle = Geometer.angle(nitrogen, hydrogen, oxygen);
        double hocAngle = Geometer.angle(hydrogen, oxygen, carbon);
        
        if (distance < 3.5 && nhoAngle > 120.0 && hocAngle > 90.0) {
            return new HBond(donor, acceptor, distance, nhoAngle, hocAngle);
        } else {
            return null;
        }
    }

    public void analyzeHBonds(Residue residue) {
        ArrayList nTerminalHBonds = residue.getNTerminalHBonds();
        ArrayList cTerminalHBonds = residue.getCTerminalHBonds();
        int numberOfNTerminalHBonds = nTerminalHBonds.size();
        int numberOfCTerminalHBonds = cTerminalHBonds.size();

        //System.err.println(numberOfNTerminalHBonds + ", " + numberOfCTerminalHBonds + " for " + residue);

        if (numberOfNTerminalHBonds == 0 && numberOfCTerminalHBonds == 0) {
            residue.setEnvironment("Loop");
        } else if (numberOfNTerminalHBonds == 1 && numberOfCTerminalHBonds == 0) {
            if (((HBond) nTerminalHBonds.get(0)).hasHelixResidueSeparation()) {
                residue.setEnvironment("End of a Helix");
            }
        } else if (numberOfNTerminalHBonds == 1 && numberOfCTerminalHBonds == 0) {
            if (((HBond) cTerminalHBonds.get(0)).hasHelixResidueSeparation()) {
                residue.setEnvironment("Start of a Helix");
            }
        } else if (numberOfNTerminalHBonds == 1 && numberOfCTerminalHBonds == 1) {
            int n = ((HBond) nTerminalHBonds.get(0)).getResidueSeparation();
            int c = ((HBond) cTerminalHBonds.get(0)).getResidueSeparation();
            //System.err.println(n + ", " + c + " for " + residue);

            if (n == 4 && c == 4) {
                residue.setEnvironment("Middle of a Helix");
            } else if (Math.abs(n - c) == 2) {
                residue.setEnvironment("Parallel Strand");
            } else if (n - c == 0) {
                residue.setEnvironment("AntiParallel Strand");
            } else if (n == 4) {
                residue.setEnvironment("End of a Helix");
            } else if (c == 4) {
                residue.setEnvironment("Start of a Helix");
            }
        }
    }

    public void buildSSES(Chain chain) {
        chain.addBackboneSegment(new Terminus("N Terminus", 'N'));
        this.buildHelices(chain);
        this.buildStrands(chain);
        //this.buildLoops(chain);
        chain.sortBackboneSegments();
        chain.addBackboneSegment(new Terminus("C Terminus", 'C'));    //and add a C Terminus for luck
    }

    public void buildHelices(Chain chain) {
        Iterator residueIterator = chain.residueIterator();

        int helixStartIndex = -1;
        int helixEndIndex = -1;
        while (residueIterator.hasNext()) {
            Residue residue = (Residue) residueIterator.next();
            int residueAbsoluteNumber = residue.getAbsoluteNumber();

            if (residue.getEnvironment().equals("Middle of a Helix")) {

                // this may be the first residue in the helix core
                if (helixStartIndex == -1) {
                    helixStartIndex = residueAbsoluteNumber - 4;
                    helixEndIndex = residueAbsoluteNumber + 4;

                // otherwise, see if we have to extend the helix
                } else {
                    if (residueAbsoluteNumber > helixEndIndex) {
                        helixEndIndex++;
                    }
                }
            } else {
                if (helixEndIndex != -1 && residueAbsoluteNumber > helixEndIndex) {
                    chain.createHelix(helixStartIndex, helixEndIndex);
                    helixStartIndex = helixEndIndex = -1;
                }
            }
        }
    }

    public void buildStrands(Chain chain) {
        Iterator residueIterator = chain.residueIterator();

        boolean lastResidueWasStrand = false;
        BackboneSegment currentStrand = new Strand();
        while (residueIterator.hasNext()) {
            Residue residue = (Residue) residueIterator.next();
            String residueEnvironment = residue.getEnvironment();
            if (residueEnvironment.equals("Parallel Strand") || residueEnvironment.equals("AntiParallel Strand")) {
                lastResidueWasStrand = true;
                currentStrand.expandBy(residue);
            } else if (residueEnvironment.equals("Loop")) {
                if (lastResidueWasStrand) {
                    currentStrand.expandBy(residue);
                } else {
                    if (currentStrand.length() != 0) {
                        chain.addBackboneSegment(currentStrand);
                        currentStrand = new Strand();
                    }
                }
                lastResidueWasStrand = false;
            } else {
                if (currentStrand.length() != 0) {
                    chain.addBackboneSegment(currentStrand);
                    currentStrand = new Strand();
                }
                lastResidueWasStrand = false;
            }
        }
    }

    public void buildLoops(Chain chain) {
        ListIterator itr = chain.backboneSegmentListIterator();

        while (itr.hasNext()) {
            BackboneSegment backboneSegment = (BackboneSegment) itr.next();
            if (backboneSegment.length() < 2) {
                continue;
            }

            if (itr.hasNext()) {
                BackboneSegment nextSegment = (BackboneSegment) itr.next();
                int startOfLoop = backboneSegment.lastResidue().getAbsoluteNumber() + 1;
                int endOfLoop = nextSegment.firstResidue().getAbsoluteNumber() - 1;
                chain.createLoop(startOfLoop, endOfLoop);
            } 
        }
    }

    public static void main(String[] args) {
        try {

            StructureFinder structureFinder = new StructureFinder(args[0]);
            Protein protein = structureFinder.getProtein();
            System.out.println(protein.toString());

            HashMap cathChainDomainMap = CATHDomainFileParser.parseUpToParticularID(args[1], protein.getID());
            HashMap chainDomainStringMap = protein.toTopsDomainStrings(cathChainDomainMap);

            Iterator itr = chainDomainStringMap.keySet().iterator();
            while (itr.hasNext()) {
                String chainID = (String) itr.next();
                HashMap domainStrings = (HashMap) chainDomainStringMap.get(chainID);
                Iterator itr2 = domainStrings.keySet().iterator();
                while (itr2.hasNext()) {
                    String domainString = (String) domainStrings.get((String) itr2.next());
                    System.out.println(protein.getID() + domainString);
                }
            }

        } catch (IOException ioe) {
            System.out.println(ioe);
        }
    }
}
