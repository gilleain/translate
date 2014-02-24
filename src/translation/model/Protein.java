package translation.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import translation.ChainDomainMap;
import translation.StructureFinder;

public class Protein {

    private String id;
    private ArrayList chains;

    public Protein() {
        this.id = "";
        this.chains = new ArrayList();
    }

    public Protein(String id) {
        this();
        this.id = id.toLowerCase();
    }

    public String getID() {
        return this.id;
    }

    public void addChain(Chain chain) {
        this.chains.add(chain);
    }

    public Iterator chainIterator() {
        return this.chains.iterator();
    }

    public void findStructure(StructureFinder structureFinder) {
        for (int i = 0; i < this.chains.size(); i++) {
            structureFinder.findStructure((Chain) this.chains.get(i));
        }
    }
    
    public HashMap getBackboneSegmentsByDomain(HashMap chainDomainMap) {
        HashMap chainBackboneSegmentMap = new HashMap();
        for (int i = 0; i < this.chains.size(); i++) {
            Chain chain = (Chain) this.chains.get(i);
            String label = chain.getCathCompatibleLabel();
            
            ArrayList domains;
            if (chainDomainMap.containsKey(label)) {
                domains = (ArrayList) chainDomainMap.get(label);
            } else {
                domains = (ArrayList) chainDomainMap.get("A");
            }
            HashMap domainSegmentMap = new HashMap();
            for (int j = 0; j < domains.size(); j++) {
                Domain domain = (Domain) domains.get(j);
                ArrayList segments = chain.filterBackboneSegmentsByDomain(domain);
                //System.err.println(segments.size() + " segments for domain " + domain);
                domainSegmentMap.put(domain.getID(), segments);
            }
            chainBackboneSegmentMap.put(label, domainSegmentMap);
        }
        return chainBackboneSegmentMap;
    }

    public Map toTopsDomainStrings(ChainDomainMap chainDomainMap) {
        Map chainDomainStringMap = new HashMap();
        for (int i = 0; i < this.chains.size(); i++) {
            Chain chain = (Chain) this.chains.get(i);
            chainDomainStringMap.put(chain.getCathCompatibleLabel(), 
            						 chain.toTopsDomainStrings(chainDomainMap));
        }
        return chainDomainStringMap;
    }

    public String[] toTopsChainStringArray() {
        String[] chainStrings = new String[this.chains.size()];
        for (int i = 0; i < this.chains.size(); i++) {
            Chain chain = (Chain) this.chains.get(i);
            chainStrings[i] = chain.getCathCompatibleLabel() 
            				+ chain.toTopsDomainStrings(new ChainDomainMap());
        }
        return chainStrings;
    }

    public String toString() {
        StringBuffer stringBuffer = new StringBuffer();
        for (int i = 0; i < this.chains.size(); i++) {
            Chain chain = (Chain) this.chains.get(i);
            stringBuffer.append(chain.toString());
        }
        return stringBuffer.toString();
    }

}
