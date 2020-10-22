
import java.io.*;
import java.util.*;
import java.util.zip.*;

public class RFRegressor{
  public static boolean displayProgress = false;
  public static final int nSplitVals = 10;
  
  private int nTrees, nNodes;
  private int[] roots, nodeLeft;
  private short[] splitFeature;
  private float[] splitVal;
  
  public void train(final float[][] features, final float[] values, 
       final int nRows, final int maxTrees, final int minRowsPerNode, 
       final int nThreads){
    nTrees = nNodes = 0;
    final int maxNodesPerTree = nRows/minRowsPerNode*2;
    roots = new int[maxTrees];
    int maxNodes = maxTrees*maxNodesPerTree;
    nodeLeft = new int[maxNodes];
    splitFeature = new short[maxNodes];
    splitVal = new float[maxNodes];
    
    Thread[] threads = new Thread[nThreads];
    for (int i = 0; i < nThreads; i++){
      final int iStart = i;
      threads[i] = new Thread(){
        @Override public void run(){
          for (int i = iStart; i < roots.length; i+= nThreads){
            RegressionNode root = new RegressionTree(features, values,
                 nRows, i, minRowsPerNode, maxNodes).root;
            synchronized (RFRegressor.this){ add(root); }
          }
        }
      };
      threads[i].start();
    }
    
    try{
      for (int i = 0; i < nThreads; i++){ 
        threads[i].join(); threads[i] = null;
      }
    }catch (InterruptedException e){
      System.err.println(e);
    }
  }
  
  public double predict(float[] features){
    double ret = 0;
    for (int root : roots){
      ret+= classify(root, features);
    }
    return ret/roots.length;
  }
  
  private double classify(int pos, float[] features){
    while (true){
      int sf = splitFeature[pos];
      if (sf < 0) return splitVal[pos];
      pos = features[sf] < splitVal[pos] ? nodeLeft[pos] : nodeLeft[pos] + 1;
    }
  }
  
  public static void saveToGZip(String rfFile, RFRegressor...rfList){
    try (DataOutputStream out = new DataOutputStream(new GZIPOutputStream(
         new FileOutputStream(rfFile), 1_000_000))){
      for (RFRegressor rf : rfList) rf.write(out);
    }catch(Exception e){
      System.err.println(e);
    }
  }
  
  public static void loadFromGZip(String rfFile, RFRegressor...rfList){
    try (DataInputStream in = new DataInputStream(new GZIPInputStream(
         new FileInputStream(rfFile), 1_000_000))){
      for (int i = 0; i < rfList.length; i++){ 
        rfList[i] = new RFRegressor();
        rfList[i].read(in);
      }
    }catch(Exception e){
      System.err.println(e);
    }
  }
  
  public void read(DataInputStream in) throws Exception{
    nTrees = in.readInt();
    roots = new int[nTrees];
    for (int i = 0; i < nTrees; i++) roots[i] = in.readInt();
    nNodes = in.readInt();
    nodeLeft = new int[nNodes];
    splitFeature = new short[nNodes];
    splitVal = new float[nNodes];
    for (int i = 0; i < nNodes; i++){
      splitFeature[i] = in.readShort();
      nodeLeft[i] = in.readInt();
      splitVal[i] = in.readFloat();
    }
  }
  
  public void write(DataOutputStream out) throws Exception{
    out.writeInt(nTrees);
    for (int i = 0; i < nTrees; i++) out.writeInt(roots[i]);
    out.writeInt(nNodes);
    for (int i = 0; i < nNodes; i++){
      out.writeShort(splitFeature[i]);
      out.writeInt(nodeLeft[i]);
      out.writeFloat(splitVal[i]);
    }
  }
  
  public synchronized void add(RegressionNode root){
    int rt = roots[nTrees++] = nNodes;
    nNodes++;
    expand(root, rt);
  }
  
  private void expand(RegressionNode node, int pos){
    if (node.left == null){
      nodeLeft[pos] = -1;
      splitFeature[pos] = -1;
      splitVal[pos] = node.average;
    }else{
      int l = nodeLeft[pos] = nNodes;
      nNodes+= 2;
      splitFeature[pos] = (short) node.splitFeature;
      splitVal[pos] = node.splitVal;
      expand(node.left, l);
      expand(node.right, l + 1);
    }
  }
  
  public class RegressionTree{
    private RegressionNode root;
    private final SplittableRandom rnd;
    
    RegressionTree(float[][] features, float[] values, int nRows, 
         int index, int minRowsPerNode, int maxNodes){
      RegressionNode[] nodes = new RegressionNode[maxNodes + 2];
      int nUsedRows = Math.min(nRows, Math.max(5, nRows/3)),
          nUsedFeatures = 
           Math.min(features.length, Math.max(5, features.length/3));
      
      rnd = new SplittableRandom(197209091220L + index);
      int[] weight = new int[nRows];
      for (int i = 0; i < nUsedRows; i++) weight[rnd.nextInt(nRows)] = 1;
      int nSelected = 0;
      for (int i = 0; i < nRows; i++) if (weight[i] > 0) nSelected++;
      int[] selectedRows = new int[nSelected];
      nSelected = 0;
      float sum = 0, sum2 = 0;
      int wSum  = 0;
      for (int i = 0; i < nRows; i++){
        if (weight[i] > 0){
          selectedRows[nSelected++] = i;
          float v = values[i];
          sum+= v; sum2+= v*v; wSum++;
        }
      }
      float average = sum/wSum, error = wSum==0 ? 0 : sum2 - sum*sum/wSum;
      root = new RegressionNode(nRows, average, error, 0, nSelected-1);
      int nodeCnt = 0;
      nodes[nodeCnt++] = root;
      float[] splitVals = new float[nSplitVals];
      for (int i = 0; i < nodeCnt && nodeCnt < maxNodes; i++){
        RegressionNode node = nodes[i];
        if (node.nRows < minRowsPerNode*2) continue;
        float maxSplitGain = 0, bestSplitVal = 0;
        int bestSplitFeature = -1;
        for (int j = 0; j < nUsedFeatures; j++){
          int splitFeature = rnd.nextInt(features.length);
          float[] featuresSplitFeature = features[splitFeature];
          for (int k = 0; k < nSplitVals; k++) 
            splitVals[k] = 
                 featuresSplitFeature[randomNodeRow(selectedRows, node, rnd)];
          Arrays.sort(splitVals);
          float[] sumLeft = new float[nSplitVals+1],
                  sum2Left = new float[nSplitVals+1];
          int[] wLeft = new int[nSplitVals+1];
          for (int r = node.startRow; r <= node.endRow; r++){
            int row = selectedRows[r];
            float rowVal = featuresSplitFeature[row],
                  v = values[row];
            for (int k = 0; k < nSplitVals+1; k++){
              if (k==nSplitVals || rowVal < splitVals[k]){
                sumLeft[k]+= v;
                sum2Left[k]+= v*v;
                wLeft[k]++;
                break;
              }
            }
          }
          for (int k = 0; k < nSplitVals; k++){
            int wSumLeft = 0, wSumRight = 0;
            for (int kk = 0; kk <= k; kk++) wSumLeft+= wLeft[kk];
            if (wSumLeft < minRowsPerNode) continue;
            for (int kk = k+1; kk < nSplitVals+1; kk++) wSumRight+= wLeft[kk];
            if (wSumRight < minRowsPerNode) continue;
            
            float sumLeftSum = 0, sum2LeftSum = 0;
            for (int kk = 0; kk <= k; kk++){
              sumLeftSum+= sumLeft[kk]; sum2LeftSum+= sum2Left[kk];
            }
            float sumRightSum = 0, sum2RightSum = 0;
            for (int kk = k+1; kk < nSplitVals+1; kk++){
              sumRightSum+= sumLeft[kk]; sum2RightSum+= sum2Left[kk];
            }
            float leftError = sum2LeftSum - sumLeftSum*sumLeftSum/wSumLeft,
                  rightError = sum2RightSum - sumRightSum*sumRightSum/wSumRight,
                  splitGain = node.error - leftError - rightError; 
            if (splitGain > maxSplitGain){
              maxSplitGain = splitGain;
              bestSplitFeature = splitFeature;
              bestSplitVal = splitVals[k];
            }
          }
        }
        if (bestSplitFeature >= 0){
          int wSumLeft = 0, wSumRight = 0;
          float sumLeft = 0, sum2Left = 0, sumRight = 0, sum2Right = 0;
          int endLeft = node.endRow;
          float[] featuresSplitFeature = features[bestSplitFeature];
          for (int r = node.startRow; r <= endLeft; r++){
            int row = selectedRows[r];
            float v = values[row];
            if (featuresSplitFeature[row] < bestSplitVal){
              wSumLeft++; sumLeft+= v; sum2Left+= v*v;
            }else{
              wSumRight++; sumRight+= v; sum2Right+= v*v;
              selectedRows[r--] = selectedRows[endLeft];
              selectedRows[endLeft--] = row;
            }
          }
          node.left = new RegressionNode(wSumLeft, sumLeft/wSumLeft,
               sum2Left - sumLeft*sumLeft/wSumLeft, node.startRow, endLeft);
          node.right = new RegressionNode(wSumRight, sumRight/wSumRight, 
               sum2Right - sumRight*sumRight/wSumRight, endLeft+1, node.endRow);
          nodes[nodeCnt++] = node.left;
          nodes[nodeCnt++] = node.right;
          node.splitVal = bestSplitVal;
          node.splitFeature = bestSplitFeature;
        }
      }
      if (displayProgress)
        System.out.println("Tree " + index + ": " + nodeCnt + " nodes");
    }
    
    private int randomNodeRow(int[] rows, RegressionNode node, 
         SplittableRandom rnd){
      return rows[rnd.nextInt(node.endRow - node.startRow + 1) + node.startRow];
    }
  }
  
  public class RegressionNode{
    RegressionNode left, right;
    float splitVal, average, error;
    int splitFeature, startRow, endRow, nRows;
  
    public RegressionNode(int nRows, float average, float error, 
         int startRow, int endRow){
      this.nRows = nRows;
      this.average = average;
      this.error = error;
      this.startRow = startRow;
      this.endRow = endRow;
    }
  }
}
