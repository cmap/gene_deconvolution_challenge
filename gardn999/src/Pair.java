
import java.util.*;

public class Pair{
  String expID;
  int barcode, highGene, lowGene;
  ArrayList<Integer> fiList = new ArrayList<>();
  double mean, median, highPeak, lowPeak, predHigh, predLow;
  
  static final int nx = 50, // total number of FI values to use
                   nFeatures = nx + 10;
  
  public Pair(int barcode, String expID){ 
    this.barcode = barcode;
    this.expID = expID;
    highGene = Main.barcodeHighMap.get(barcode);
    lowGene = Main.barcodeLowMap.get(barcode);
  }
  
  // return features used in random forest
  float[] getFeatures(Set set){
    Barcode bar = set.barcodeMap.get(barcode);
    Experiment exp = set.experimentMap.get(expID);
    int n = fiList.size();
    float[] ret = new float[nFeatures];
    
    // use nx number of FI values
    for (int i = 0; i < nx; i++) ret[i] = fiList.get(i*(n-1)/(nx-1));
    
    ret[nx] = n;
    ret[nx+1] = (float)n*set.nHighRight/set.n/bar.nHighRight;
    ret[nx+2] = (float)bar.nHighRight;
    ret[nx+3] = (float)(bar.highPeak/bar.lowPeak);
    ret[nx+4] = (float)(bar.mean - bar.median);
    ret[nx+5] = (float)(bar.mean - mean);
    ret[nx+6] = (float)(bar.median - median);
    ret[nx+7] = (float)(exp.mean - exp.median);
    ret[nx+8] = (float)(exp.highPeak- exp.lowPeak);
    ret[nx+9] = (float)(set.highPeak/set.lowPeak);
    
    return ret;
  }
  
  // sort FI values and calculate mean, median, highPeak and lowPeak
  void calcFeatures(){
    int n = fiList.size();
    Collections.sort(fiList);
    median = fiList.get(n/2);
    mean = 0;
    for (int i = 0; i < n; i++) mean+= fiList.get(i);
    mean/= n;
    if (mean > median){
      highPeak = fiList.get(n/3);
      lowPeak = fiList.get(n*5/6);
    }else{
      highPeak = fiList.get(n*2/3);
      lowPeak = fiList.get(n/6);
    }
  }
  
  public void add(int fi){ fiList.add(fi); }
  
  public int size(){ return fiList.size(); }
  
  @Override public String toString(){
    String s = "barcode: " + barcode + ", ";
    s+= "high: " + highGene + ", low: " + lowGene + "\n";
    for (int fi : fiList) s+= fi + "\n";
    return s;
  }
}
