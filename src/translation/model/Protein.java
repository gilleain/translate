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
    private List<Chain> chains;

    public Protein() {
        this.id = "";
        this.chains = new ArrayList<Chain>();
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

    public Iterator<Chain> chainIterator() {
        return this.chains.iterator();
    }

    public void findStructure(StructureFinder structureFinder) {
        for (int i = 0; i < this.chains.size(); i++) {
            structureFinder.findStructure((Chain) this.chains.get(i));
        }
    }
    
    public Map<String, Map<String, List<BackboneSegment>>> getBackboneSegmentsByDomain(ChainDomainMap chainDomainMap) {
        Map<String, Map<String, List<BackboneSegment>>> chainBackboneSegmentMap = 
        		new HashMap<String, Map<String, List<BackboneSegment>>>();
        for (int i = 0; i < this.chains.size(); i++) {
            Chain chain = (Chain) this.chains.get(i);
            String label = chain.getCathCompatibleLabel();
            
            List<Domain> domains;
            if (chainDomainMap.containsKey(label)) {
                domains = chainDomainMap.get(label);
            } else {
                domains = chainDomainMap.get("A");
            }
            
            Map<String, List<BackboneSegment>> domainSegmentMap = new HashMap<String, List<BackboneSegment>>();
            for (int j = 0; j < domains.size(); j++) {
                Domain domain = (Domain) domains.get(j);
                List<BackboneSegment> segments = chain.filterBackboneSegmentsByDomain(domain);
                //System.err.println(segments.size() + " segments for domain " + domain);
                domainSegmentMap.put(domain.getID(), segments);
            }
            chainBackboneSegmentMap.put(label, domainSegmentMap);
        }
        return chainBackboneSegmentMap;
    }

    public Map<String, Map<String, String>> toTopsDomainStrings(ChainDomainMap chainDomainMap) {
        Map<String, Map<String, String>> chainDomainStringMap = new HashMap<String, Map<String, String>>();
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
