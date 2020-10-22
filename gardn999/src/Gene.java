
import java.util.*;

public class Gene{
  TreeMap<String, Double> peakMap = new TreeMap<>();
  int id;
  
  public Gene(int id){ this.id = id; }
  
  void addPeak(String expID, double peak){ peakMap.put(expID, peak); }
  
  @Override public String toString(){
    String s = "" + id;
    for (String expID : peakMap.keySet()){
      s+= String.format("\t%5.4f", peakMap.get(expID));
    }
    s+= "\n";
    return s;
  }
  
}
