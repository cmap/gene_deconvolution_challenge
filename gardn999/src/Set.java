
import java.io.*;
import java.util.*;

public class Set{
  // maximum number threads for utilizing multiple cpu cores
  static final int nThreads = 12;
  
  String dirPath;
  TreeMap<String, Experiment> experimentMap;
  TreeMap<Integer, Barcode> barcodeMap;
  
  // for feature creation
  int n = 0, nHighRight = 0;
  double mean = 0, median = 0, highPeak = 0, lowPeak = 0;
  
  public void init(String dirPath){
    this.dirPath = dirPath;
    
    File[] files = new File(dirPath).listFiles();
    if (files == null) return;
    
    // speed up creation of new Experiment objects by using multiple threads
    Thread[] threads = new Thread[nThreads];
    
    experimentMap = new TreeMap<>();
    for (int iThread = 0; iThread < nThreads; iThread++){
      final int iStart = iThread;
      threads[iThread] = new Thread(){
        @Override public void run(){
          for (int i = iStart; i < files.length; i+= nThreads){
            String path = files[i].getPath(),
                   id = path.substring(path.length()-7, path.length()-4);
            Experiment exp = new Experiment(id, path);
            synchronized(Set.this){ experimentMap.put(id, exp); }
          }
        }
      };
      threads[iThread].start();
    }
      
    // wait for all threads to complete
    try{
      for (int iThread = 0; iThread < nThreads; iThread++){ 
        threads[iThread].join(); threads[iThread] = null;
      }
    }catch (InterruptedException e){
      System.err.println(e);
    }
    
    // fill map<barcode id, Barcode>
    barcodeMap = new TreeMap<>();
    for (String e : experimentMap.keySet()){
      for (int barcode : Main.barcodeSet){
        Pair pair = experimentMap.get(e).pairMap.get(barcode);
        if (!barcodeMap.containsKey(barcode))
          barcodeMap.put(barcode, new Barcode());
        Barcode b = barcodeMap.get(barcode);
        b.n++;
        if (pair.highPeak > pair.lowPeak) b.nHighRight++;
        b.highPeak+= pair.highPeak; b.lowPeak+= pair.lowPeak;
        b.mean+= pair.mean; b.median+= pair.median;
      }
    }
    
    // calculate quantities used in Pair.getFeatures
    n = experimentMap.size()*barcodeMap.size();
    for (Integer i : barcodeMap.keySet()){
      Barcode b = barcodeMap.get(i);
      b.calc();
      nHighRight+= b.nHighRight;
      mean+= b.mean*b.n; median+= b.median*b.n;
      highPeak+= b.highPeak*b.n; lowPeak+= b.lowPeak*b.n;
    }
    mean/= n; median/= n; highPeak/= n; lowPeak/= n;
  }

  public void calcPeaks(RFRegressor[] rf){
    Set set = this;
    
    // speed up the random forest peak predictions by using multiple threads
    Thread[] threads = new Thread[nThreads];
    
    for (int iThread = 0; iThread < nThreads; iThread++){
      final int iStart = iThread;
      final List<String> experiments = new ArrayList<>(experimentMap.keySet());
      threads[iThread] = new Thread(){
        @Override public void run(){
          for (int i = iStart; i < experiments.size(); i+= nThreads){
            Experiment exp = experimentMap.get(experiments.get(i));
            exp.calcPeaks(rf, set);
          }
        }
      };
      threads[iThread].start();
    }
    
    // wait for all threads to complete
    try{
      for (int iThread = 0; iThread < nThreads; iThread++){ 
        threads[iThread].join(); threads[iThread] = null;
      }
    }catch (InterruptedException e){
      System.err.println(e);
    }
  }
  
  @Override public String toString(){
    String s = dirPath + "\n";
    for (Experiment e : experimentMap.values())
      s+= e.id + ": " + e.pairMap.size() + " pairs\n";
    return s;
  }
  
}
