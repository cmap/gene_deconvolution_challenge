
public class Barcode{
  int n = 0,  nHighRight = 0;
  double mean = 0, median = 0, highPeak = 0, lowPeak = 0;
  
  // calculate quantities used in Pair.getFeatures
  void calc(){ mean/= n; median/= n; highPeak/= n; lowPeak/= n; }
  
  @Override public String toString(){
    return String.format(
            "high right: %3d/%3d, high: %8.2f, low: %8.2f",
            nHighRight, n, highPeak, lowPeak);
  }
}
