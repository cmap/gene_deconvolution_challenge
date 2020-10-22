
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import javax.imageio.*;

public class Main {
  // random forest parameters
  static final int nTrees = 50, minRows = 15;
  
  // default paths
  static String 
          inputDPK = "input/DPK.CP001_A549_24H_X1_B42/",
          inputLITMUS = "input/LITMUS.KD017_A549_96H_X1_B42/",
          subDPK = "DPK.CP001_A549_24H_X1_B42.gct",
          subLITMUS = "LITMUS.KD017_A549_96H_X1_B42.gct",
          truthDPK = "ground-truth/DPK.CP001_A549_24H_X1_B42_DECONV_UNI.gct",
          truthLITMUS = 
            "ground-truth/LITMUS.KD017_A549_96H_X1_B42_DECONV_UNI.gct";
  
  // for barcode_to_gene_map.txt
  static final String pathBarcodeToGene = "barcode_to_gene_map.txt";
  
  static final HashMap<Integer, Integer> barcodeHighMap = new HashMap<>(),
                                         barcodeLowMap = new HashMap<>();
  static final HashSet<Integer> barcodeSet = new HashSet<>();
  
  // initialize barcode to gene maps
  static{
    try (BufferedReader br = 
            new BufferedReader(new FileReader(pathBarcodeToGene))){
      for (String line; (line = br.readLine()) != null; ){
        String[] s = line.split("\t");
        if (!s[0].equals("barcode_id")){
          int barcode_id = Integer.parseInt(s[0]),
              gene_id = Integer.parseInt(s[1]),
              high_prop = Integer.parseInt(s[2]);
          if (high_prop == 0) barcodeLowMap.put(barcode_id, gene_id);
          else barcodeHighMap.put(barcode_id, gene_id);
        }
        
        for (int barcode_id : barcodeLowMap.keySet())
          if (barcodeHighMap.containsKey(barcode_id))
            barcodeSet.add(barcode_id);
      }
    }catch (Exception e){ System.err.println(e); }
  }
  
  // create map<gene id, Gene> from Set object
  static TreeMap<Integer, Gene> getPredMap(Set set){
    TreeMap<Integer, Gene> geneMap = new TreeMap<>();
    for (String id : set.experimentMap.keySet()){
      Experiment e = set.experimentMap.get(id);
      for (Integer pairId : e.pairMap.keySet()){
        Pair pair = e.pairMap.get(pairId);
        if (!geneMap.containsKey(pair.lowGene))
          geneMap.put(pair.lowGene, new Gene(pair.lowGene));
        geneMap.get(pair.lowGene).addPeak(id, pair.predLow);
        
        if (!geneMap.containsKey(pair.highGene))
          geneMap.put(pair.highGene, new Gene(pair.highGene));
        geneMap.get(pair.highGene).addPeak(id, pair.predHigh);
      }
    }
    return geneMap;
  }
  
  // read truth file into map<gene id, Gene>
  static TreeMap<Integer, Gene> getTruthMap(String truthPath){
    TreeMap<Integer, Gene> truthMap = new TreeMap<>();
    try (BufferedReader br = new BufferedReader(new FileReader(truthPath))){
      String[] expIDs = null;
      for (String line; (line = br.readLine()) != null; ){
        String[] s = line.trim().split("\t");
        if (s[0].equals("id")){
          expIDs = s;
        }else if (s.length > 4){
          int id = Integer.parseInt(s[0]);
          truthMap.put(id, new Gene(id));
          Gene gene = truthMap.get(id);
          for (int i = 1; i < s.length; i++)
            gene.addPeak(expIDs[i], Double.parseDouble(s[i]));
        }
      }
    }catch (Exception e){ System.err.println(e); }
    return truthMap;
  }
  
  // plot the FI distribution for a pair
  // - truth peaks are black lines and predictions are blue
  // - low peaks are half height and high peaks are full height
  static void plotPair(int barcode, String exp, String inputFile, 
                       String truthFile, String rfFile){
    RFRegressor[] rf = new RFRegressor[2];
    RFRegressor.loadFromGZip(rfFile, rf);
    Set set = new Set(); 
    set.init(inputFile); 
    set.calcPeaks(rf);
    TreeMap<Integer, Gene> truthMap = getTruthMap(truthFile);
      
    Experiment e = set.experimentMap.get(exp);
    Pair p = e.pairMap.get(barcode);
    float[] f = p.getFeatures(set);
    double predLow = rf[0].predict(f), predHigh = rf[1].predict(f);
    
    System.out.println("low gene: " + p.lowGene + ", high gene: " + p.highGene);
    double truthLow = truthMap.get(p.lowGene).peakMap.get(exp),
           truthHigh = truthMap.get(p.highGene).peakMap.get(exp);
    System.out.println(String.format("low peak  - pred:%3.2f, truth:%3.2f",
            predLow, truthLow));
    System.out.println(String.format("high peak - pred:%3.2f, truth:%3.2f", 
            predHigh, truthHigh));
    
    // plot parameters
    int W = 800, H = 800, nBins = 50;
    double maxFI = 5000;
    BufferedImage image = 
            new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);
    Graphics g = image.getGraphics();
    
    // set background white
    g.setColor(Color.white);
    g.fillRect(0, 0, W, H);
    
    // draw FI counts
    g.setColor(Color.orange);
    int[] n = new int[nBins];
    for (int i : p.fiList) if (i < maxFI) n[(int)(nBins*i/maxFI)]++;
    for (int i = 0; i < nBins; i++)
      g.fillRect(i*W/nBins, H - H*n[i]/20 + 1, W/nBins - 2, H*n[i]/20);
    
    // draw peaks
    g.setColor(new Color(0, 200, 255));
    g.fillRect((int)(predLow*W/maxFI - 2), H/2, 5, H/2);
    g.fillRect((int)(predHigh*W/maxFI - 2), 0, 5, H);
    
    g.setColor(Color.black);
    g.fillRect((int)(truthLow*W/maxFI - 1), H/2, 3, H/2);
    g.fillRect((int)(truthHigh*W/maxFI - 1), 0, 3, H);
    
    // draw border
    g.setColor(Color.black);
    g.drawRect(0, 0, W-1, H-1);
    
    g.dispose();
    
    try{ 
      ImageIO.write(image, "png", new File("pair.png")); 
    }catch(IOException ex){ 
      System.err.println(ex); 
    }
  }
  
  // output the .gct file
  static void doSubmission(TreeMap<Integer, Gene> geneMap, String subFile){
    try (PrintWriter pw = new PrintWriter(subFile)){
      pw.write("#1.3\n");
      Gene first = geneMap.firstEntry().getValue();
      pw.write(geneMap.size() + "\t" + first.peakMap.size() + "\t0\t0\n");
      pw.write("id");
      for (String id : first.peakMap.keySet()) pw.write("\t" + id);
      pw.write("\n");
      for (Integer id : geneMap.keySet()) pw.write(geneMap.get(id).toString());
    }catch (Exception exception){ System.err.println(exception); }
  }

  // Makes a rough estimate of prediction quality by displaying how close
  // on average the peak predictions are to the truth and the percent of the
  // time the peaks are on the correct side relative to each other.
  static void doScore(TreeMap<Integer, Gene> predMap, 
                      TreeMap<Integer, Gene> truthMap, String label){
    int nCorrectSide = 0, 
        n = predMap.size()*predMap.firstEntry().getValue().peakMap.size();
    double[] diffSum = new double[2];
    
    for (int barcode : barcodeSet){
      int idHigh = barcodeHighMap.get(barcode),
          idLow = barcodeLowMap.get(barcode);
      Gene genePredHigh = predMap.get(idHigh), 
           genePredLow = predMap.get(idLow),
           geneTruthHigh = truthMap.get(idHigh),
           geneTruthLow = truthMap.get(idLow);
      for (String exp : genePredHigh.peakMap.keySet()){
        double predHigh = genePredHigh.peakMap.get(exp), 
               predLow = genePredLow.peakMap.get(exp),
               truthHigh = geneTruthHigh.peakMap.get(exp),
               truthLow = geneTruthLow.peakMap.get(exp);
        if ((truthHigh-truthLow)*(predHigh-predLow) >= 0) nCorrectSide++;
      
        diffSum[0]+= (truthLow - predLow)*(truthLow - predLow);
        diffSum[1]+= (truthHigh - predHigh)*(truthHigh - predHigh);
      }
    }
    
    System.out.println(String.format(label + 
            " low:%7.2f, high:%7.2f, avg:%7.2f, side:%6.2f%%", 
            Math.sqrt(diffSum[0]*2/n), Math.sqrt(diffSum[1]*2/n), 
            Math.sqrt((diffSum[0] + diffSum[1])/n), 
            100.0*nCorrectSide*2/n ));
  }

  // add data to features and truth arrays for training
  static void doTrainArrays(TreeMap<Integer, Gene> truthMap, Set set, 
                                float[][] features, float[][] truth,
                                int start){
    int row = start;
    for (String id : set.experimentMap.keySet()){
      Experiment e = set.experimentMap.get(id);
      for (Integer pairId : e.pairMap.keySet()){
        Pair pair = e.pairMap.get(pairId);
        truth[0][row] = 
                (float)(double)truthMap.get(pair.lowGene).peakMap.get(e.id);
        truth[1][row] = 
                (float)(double)truthMap.get(pair.highGene).peakMap.get(e.id);
        float[] f = pair.getFeatures(set);
        for (int i = 0; i < f.length; i++) features[i][row] = f[i];
        row++;
      }
    }
  }
  
  // for creation and testing of random forest model using both sets
  static void runRF(boolean train){
    TreeMap<Integer, Gene> truthDPKMap = getTruthMap(truthDPK),
                           truthLITMUSMap = getTruthMap(truthLITMUS);
    RFRegressor[] rfAll = new RFRegressor[2];
    Set setDPK = new Set(), setLITMUS = new Set();
    setDPK.init(inputDPK); setLITMUS.init(inputLITMUS);
    
    // train the model and write out to "rfAll.gz"
    if (train){
      int nDPK = truthDPKMap.size()*setDPK.experimentMap.size()/2,
          nLITMUS = truthLITMUSMap.size()*setLITMUS.experimentMap.size()/2,
          nFeatures = Pair.nFeatures;
      float[][] features = new float[nFeatures][nDPK+nLITMUS], 
                truth = new float[2][nDPK+nLITMUS];
      doTrainArrays(truthDPKMap, setDPK, features, truth, 0);
      doTrainArrays(truthLITMUSMap, setLITMUS, features, truth, nDPK);
      
      for (int i = 0; i < 2; i++){
        rfAll[i] = new RFRegressor();
        rfAll[i].train(features, truth[i], nDPK+nLITMUS, nTrees, minRows, 4);
      }
      RFRegressor.saveToGZip("rfAll.gz", rfAll);
      
    // load the model from "rfAll.gz"
    }else{
      RFRegressor.loadFromGZip("rfAll.gz", rfAll);
    }
    
    // make peak predictions
    setDPK.calcPeaks(rfAll); setLITMUS.calcPeaks(rfAll);
    
    // write out submission files print out score estimates
    TreeMap<Integer, Gene> predMapLITMUS = getPredMap(setLITMUS);
    doSubmission(predMapLITMUS, subLITMUS);
    doScore(predMapLITMUS, truthLITMUSMap, "train:all/test:LIT)");
    
    TreeMap<Integer, Gene> predMapDPK = getPredMap(setDPK);
    doSubmission(predMapDPK, subDPK);
    doScore(predMapDPK, truthDPKMap, "train:all/test:DPK)");
  }
  
  // for testing random forest feature selection 
  // train -> test with other set: DPK -> LITMUS and LITMUS -> DPK
  static void runRFSplit(boolean train){
    TreeMap<Integer, Gene> truthDPKMap = getTruthMap(truthDPK),
                           truthLITMUSMap = getTruthMap(truthLITMUS);
    RFRegressor[] rfDPK = new RFRegressor[2];
    RFRegressor[] rfLITMUS = new RFRegressor[2];
    Set setDPK = new Set(), setLITMUS = new Set();
    setDPK.init(inputDPK); setLITMUS.init(inputLITMUS);
    
    // train the models and write out to "rfDPK.gz" and "rfLITMUS.gz"
    if (train){
      int nDPK = truthDPKMap.size()*setDPK.experimentMap.size()/2,
          nLITMUS = truthLITMUSMap.size()*setLITMUS.experimentMap.size()/2,
          nFeatures = Pair.nFeatures;
      float[][] dpkFeatures = new float[nFeatures][nDPK], 
                dpkTruth = new float[2][nDPK],
                litmusFeatures = new float[nFeatures][nLITMUS], 
                litmusTruth = new float[2][nLITMUS];
      doTrainArrays(truthDPKMap, setDPK, dpkFeatures, dpkTruth, 0);
      doTrainArrays(truthLITMUSMap, setLITMUS, litmusFeatures, litmusTruth, 0);
      
      for (int i = 0; i < 2; i++){
        rfDPK[i] = new RFRegressor();
        rfDPK[i].train(dpkFeatures, dpkTruth[i], nDPK, nTrees, minRows, 4);
        rfLITMUS[i] = new RFRegressor();
        rfLITMUS[i].train(
                litmusFeatures, litmusTruth[i], nLITMUS, nTrees, minRows, 4);
      }
      RFRegressor.saveToGZip("rfDPK.gz", rfDPK);
      RFRegressor.saveToGZip("rfLITMUS.gz", rfLITMUS);
      
    // load the models from "rfDPK.gz" and "rfLITMUS.gz"
    }else{
      RFRegressor.loadFromGZip("rfDPK.gz", rfDPK);
      RFRegressor.loadFromGZip("rfLITMUS.gz", rfLITMUS);
    }
    
    // make peak predictions using other model file
    setDPK.calcPeaks(rfLITMUS); setLITMUS.calcPeaks(rfDPK);
    
    // write out submission files print out score estimates
    TreeMap<Integer, Gene> predMapLITMUS = getPredMap(setLITMUS);
    doSubmission(predMapLITMUS, subLITMUS);
    doScore(predMapLITMUS, truthLITMUSMap, "train:DPK/test:LIT)");

    TreeMap<Integer, Gene> predMapDPK = getPredMap(setDPK);
    doSubmission(predMapDPK, subDPK);
    doScore(predMapDPK, truthDPKMap, "train:LIT/test:DPK)");
  }
  
  // make sure directory has a slash at the end
  static String checkSlash(String path){
    char endChar = path.charAt(path.length()-1);
    return endChar == '\\' || endChar == '/' ? path : path + "/";
  }
  
  // train or make predictions using command-line arguments
  static void runCommandLine(String[] args){
    long t0 = System.nanoTime();
    
    // make predictions
    // arguments: <input directory> <output directory> <plate name>
    if (args.length == 3){
      System.out.println("Starting tester");
      
      String inputDir = args[0], outputDir = checkSlash(args[1]), 
      fileName = args[2], outPath = outputDir + fileName;
      System.out.println("Input directory: " + inputDir + 
                       ", Creating " + outPath);
      
      System.out.println("Reading data");
      RFRegressor[] rfAll = new RFRegressor[2];
      RFRegressor.loadFromGZip("rfAll.gz", rfAll);
      Set set = new Set(); set.init(inputDir); set.calcPeaks(rfAll);
      
      System.out.println("Making predictions");
      TreeMap<Integer, Gene> predMap = getPredMap(set);
      doSubmission(predMap, outPath);
      
      System.out.println(String.format("%s created (%3.2fs)", 
              outPath, (System.nanoTime()-t0)/1.0e9));
      
    // create "rfAll.gz" using DPK and LITMUS training sets and truth files
    // arguments: <DPK input directory> <DPK truth file>
    //            <LITMUS input directory> <LITMUS truth file>
    }else if (args.length == 4){
      System.out.println("Starting trainer");
      
      inputDPK = args[0]; truthDPK = args[1];
      inputLITMUS = args[2]; truthLITMUS = args[3];
      
      System.out.println("Reading data");
      
      TreeMap<Integer, Gene> truthDPKMap = getTruthMap(truthDPK),
                             truthLITMUSMap = getTruthMap(truthLITMUS);
      Set setDPK = new Set(), setLITMUS = new Set();
      setDPK.init(inputDPK); setLITMUS.init(inputLITMUS);
      
      // create training arrays
      int nDPK = truthDPKMap.size()*setDPK.experimentMap.size()/2,
          nLITMUS = truthLITMUSMap.size()*setLITMUS.experimentMap.size()/2,
          nFeatures = Pair.nFeatures;
      float[][] features = new float[nFeatures][nDPK+nLITMUS], 
                truth = new float[2][nDPK+nLITMUS];
      doTrainArrays(truthDPKMap, setDPK, features, truth, 0);
      doTrainArrays(truthLITMUSMap, setLITMUS, features, truth, nDPK);
      
      // train the model and write out to "rfAll.gz"
      System.out.println("Training model");
      
      RFRegressor[] rfAll = new RFRegressor[2];
      for (int i = 0; i < 2; i++){
        rfAll[i] = new RFRegressor();
        rfAll[i].train(features, truth[i], nDPK+nLITMUS, nTrees, minRows, 4);
      }
      RFRegressor.saveToGZip("rfAll.gz", rfAll);
      
      System.out.println(String.format("rfAll.gz created (%3.2fs)", 
              (System.nanoTime()-t0)/1.0e9));
    
    // display instructions if incorrect number of arguments
    }else{
      System.out.println("Incorrect number of arguments");
      System.out.println("Trainer: java -Xmx2G -jar gardn999_CMap_DPeak.jar " +
              "<DPK input dir> <DPK truth file> <LITMUS input dir> " +
              "<LITMUS truth file>");
      System.out.println("Tester: java -Xmx2G -jar gardn999_CMap_DPeak.jar " +
              "<input dir> <output dir> <plate name>");
    }
  }
  
  public static void main(String[] args){
    //// train or make predictions using command-line arguments ////////////////
    runCommandLine(args);
    
    //// for testing ///////////////////////////////////////////////////////////
    //runRF(true); 
    //runRFSplit(true);
    //plotPair(62, "C15", inputDPK, truthDPK, "rfLITMUS.gz");
  }
  
}
