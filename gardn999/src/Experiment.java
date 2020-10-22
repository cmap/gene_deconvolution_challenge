
import java.io.*;
import java.util.*;


public class Experiment{
  String id;
  HashMap<Integer, Pair> pairMap;
  int n = 0, nHighRight = 0;
  double mean = 0, median = 0, highPeak = 0, lowPeak = 0;
  
  public Experiment(String id, String filePath){
    this.id = id;
    pairMap = new HashMap<>();
    
    // add pairs to pairMap from data file
    try (BufferedReader br = new BufferedReader(new FileReader(filePath))){
      for (String line; (line = br.readLine()) != null; ){
        String[] s = line.split("\t");
        if (!s[0].equals("barcode_id")){
          int barcode = Integer.parseInt(s[0]),
              FI = Integer.parseInt(s[1]);
          if (Main.barcodeSet.contains(barcode)){
            if (!pairMap.containsKey(barcode))
              pairMap.put(barcode, new Pair(barcode, id));
            pairMap.get(barcode).add(FI);
          }
        }
      }
    }catch (Exception e){ System.err.println(e); }
    
    // calculate quantities used in Pair.getFeatures
    for (Integer i : pairMap.keySet()){
      Pair p = pairMap.get(i);
      p.calcFeatures();
      if (p.highPeak > p.lowPeak) nHighRight++;
      mean+= p.mean; median+= p.median;
      highPeak+= p.highPeak; lowPeak+= p.lowPeak;
    }
    n = pairMap.size(); 
    mean/= n; median/= n; highPeak/= n; lowPeak/= n;
  }
  
  // use random forest to make peak predictions
  void calcPeaks(RFRegressor[] rf, Set set){
    for (Integer i : pairMap.keySet()){
      Pair p = pairMap.get(i);
      float[] features = p.getFeatures(set);
      p.predLow = rf[0].predict(features);
      p.predHigh = rf[1].predict(features);
    }
  }
  
  @Override public String toString(){
    String s = id + ": " + pairMap.size() + " pairs\n";
    for (int barcode : pairMap.keySet()) 
      s+= barcode + ": " + pairMap.get(barcode).size() + "\n";
    return s;
  }
}
