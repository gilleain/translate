package translation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import tops.engine.Result;
import tops.engine.drg.Explorer;
import translation.model.Chain;
import translation.model.Protein;

/*
Intended solely for testing the translation code by running on a directory
*/

public class TestRunner {

    public static void main(String[] args) {
        String directoryName = args[0];
        String strfileName   = args[1];
        String cathfileName  = args[2];
        String loggingLevel  = args[3];

        Logger.getLogger("translation.FoldAnalyser").setLevel(Level.parse(loggingLevel));

        // get the filenames from the pdb file directory
        File structureDirectory = new File(directoryName);
        String[] fileList = structureDirectory.list();

        // get the strings from the string file
        HashMap idchaindomainMap = new HashMap();
        BufferedReader bufferer = null;
        Map<String, ChainDomainMap> cathPDBIDChainDomainMap = null;

        try {
            bufferer = new BufferedReader(new FileReader(strfileName));
            String line;
            while ((line = bufferer.readLine()) != null) {
                String head = line.substring(0,6);
                String body = line.substring(6);
                String id   = head.substring(0,4);
                String ch   = head.substring(4,5);
                String dom  = head.substring(5,6);

                if (idchaindomainMap.containsKey(id)) {
                    HashMap chainMap = (HashMap) idchaindomainMap.get(id);
                    if (chainMap.containsKey(ch)) {
                        HashMap domainBodyMap = (HashMap) chainMap.get(ch);
                        if (domainBodyMap.containsKey(dom)) {
                            System.err.println("Domain->Body map already contains " + dom + " for id " + id + " chain " + ch);
                        } else {
                            domainBodyMap.put(dom, body);
                        }
                    } else {
                        HashMap newDomainMap = new HashMap();
                        newDomainMap.put(dom, body);
                        chainMap.put(ch, newDomainMap);
                    }
                } else {
                    HashMap newChainMap = new HashMap();
                    HashMap newDomainMap = new HashMap();
                    newDomainMap.put(dom, body);
                    newChainMap.put(ch, newDomainMap);
                    idchaindomainMap.put(id, newChainMap);
                }
            }

            // get the cath domain definitions
            cathPDBIDChainDomainMap = CATHDomainFileParser.parseWholeFile(cathfileName);

        } catch (IOException ioe) {
            System.err.println(ioe.toString());
        } finally {
            if (bufferer != null) {
                try {
                    bufferer.close();
                } catch (IOException ioe) {
                    System.err.println(ioe.toString());
                }
            }
        }


        // translate and match
        FoldAnalyser foldAnalyser = new FoldAnalyser();
        Explorer explorer = new Explorer();

        for (int i = 0; i < fileList.length; i++) {

            try {
                // from the name of the file, get the id 
                String filename = directoryName + "/" + fileList[i];
                String pdbid = fileList[i].substring(0, 4);

                // use this id to get the chainmaps for this id
                HashMap chainMap = (HashMap) idchaindomainMap.get(pdbid);
                ChainDomainMap cathChainMap = cathPDBIDChainDomainMap.get(pdbid);

                // translate the pdbfile
                Protein protein = foldAnalyser.analyse(PDBReader.read(filename));
                Map chainDomainStringMap = protein.toTopsDomainStrings(cathChainMap);

                Iterator itr = protein.chainIterator();

                while (itr.hasNext()) {
                    Chain chain = (Chain) itr.next();

                    // ignore DNA chains
                    if (chain.isDNA()) {
                        System.err.println(pdbid + chain.getLabel() + " is DNA");
                        continue;
                    } 

                    // otherwise, get tops domain strings, and compare
                    String chainID = chain.getCathCompatibleLabel();
                    HashMap domainStringMap = (HashMap) chainDomainStringMap.get(chainID);

                    // find the chain that has been translated in the map of dssptops strings
                    if (chainMap != null && chainMap.containsKey(chainID)) {
                        HashMap domainBodyMap = (HashMap) chainMap.get(chainID);
                        Iterator itr2 = domainBodyMap.keySet().iterator();

                        // go through the domains in the dssptops strings, getting corresponding translated versions
                        while (itr2.hasNext()) {
                            String domainID = (String) itr2.next();
                            String dsspTopsDomain = pdbid + chainID + domainID + (String) domainBodyMap.get(domainID);
                            if (domainStringMap.containsKey(domainID)) {
                                String myTopsDomain = pdbid + chainID + (String) domainStringMap.get(domainID);
                                Result result = explorer.comparePair(dsspTopsDomain, myTopsDomain);
                                System.out.println(result + " for " + dsspTopsDomain + " <=> " + myTopsDomain);
                            } else {
                                System.err.println("Domain " + domainID + " not found for pdbid " + pdbid);
                            }
                        }
                    // if this chain isn't in the map, complain
                    } else {
                        System.err.println("Chain " + chain.getCathCompatibleLabel() + " not found for pdbid " + pdbid);
                    }
                }
            } catch (Exception ex) {
                ex.printStackTrace();
                System.err.println(ex.toString());
            }
        }
    }

}
