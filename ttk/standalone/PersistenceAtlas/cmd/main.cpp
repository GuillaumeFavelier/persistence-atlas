#include"main.h"

int main(int argc, char **argv) {
  // register these arguments to the command line parser

  // Super parameters:
  int numberOfClusters=-1;

  // Pipeline parameters:
  // PersistenceMap parameters:
  int criticalType = 1;
  double persistenceThreshold = 0.01,
         persistenceScale = 0.5,
         persistenceSpread = 0.1,
         defaultDeviation = 0,
         defaultAmplitude = 0;

  // LDistanceMatrix parameters:
  string distanceType="2";
  int transposeMatrix=false;

  // SpectralEmbedding parameters:
  string seModulePath="default";
  string seModuleName="spectralEmbedding";
  string functionName="doIt";
  int embeddingDimension=2;
  int minimumNumberOfComponents=1;
  int maximumNumberOfComponents=20;
  int numberOfNeighbors=5;
  double sigma=-1;
  string metaData="ScalarField";

  // Clustering parameters:
  string clModulePath="default";
  string clModuleName="clustering";
  string eigenvalues="Eigenvalues";
  int minimumNumberOfClusters=1;
  int maximumNumberOfClusters=20;
  int currentGap=1;

  // SphereFromPoints
  double radius=0.04;
  int thetaResolution=32;
  int startTheta=0;
  int endTheta=360;
  int phiResolution=32;
  int startPhi=0;
  int endPhi=360;

  // UncertainDataEstimator parameters:
  int computeLowerBound=true;
  int computeUpperBound=true;
  string clusteringData="ClusterIdentifier";
  int binCount=5;

  // MandatoryCriticalPoints parameters:
  double mcpThreshold=0;

  // ttkPersistenceMapBatch pmb;
  MyProgram<ttkSortScalarField> ssf;

  // Super parameters:
  ssf.parser_.setArgument("numberOfClusters", &numberOfClusters,
      "Force number of clusters (-1 for auto-detection)", true);

  // Pipeline parameters:
  ssf.parser_.setArgument("critical", &criticalType,
      "Critical point type (0: min, 1: max, 2: both)", true);
  ssf.parser_.setArgument("threshold", &persistenceThreshold,
      "Persistence threshold", true);
  // ssf.parser_.setArgument("scale", &persistenceScale,
  //     "Persistence scale (amplitude)", true);
  ssf.parser_.setArgument("spread", &persistenceSpread,
      "Persistence spread (standard deviation)", true);
  // ssf.parser_.setArgument("amplitude", &defaultAmplitude,
  //     "Default amplitude", true);
  // ssf.parser_.setArgument("deviation", &defaultDeviation,
  //     "Default standard deviation", true);

  // ssf.parser_.setArgument("distance", &distanceType,
  //     "Distance type", true);
  // ssf.parser_.setArgument("transpose", &transposeMatrix,
  //     "Transpose matrix", true);

  // ssf.parser_.setArgument("seModulePath", &seModulePath,
  //     "Spectral Embedding Module Path", true);
  // ssf.parser_.setArgument("seModuleName", &seModuleName,
  //     "Spectral Embedding Module Name", true);
  // ssf.parser_.setArgument("functionName", &functionName,
  //     "Function Name", true);
  // ssf.parser_.setArgument("embeddingDimension", &embeddingDimension,
  //     "Embedding Dimension", true);
  // ssf.parser_.setArgument("minimumNumberOfComponents", &minimumNumberOfComponents,
  //     "Minimum Number Of Components", true);
  // ssf.parser_.setArgument("maximumNumberOfComponents", &maximumNumberOfComponents,
  //     "Maximum Number Of Components", true);
  ssf.parser_.setArgument("numberOfNeighbors", &numberOfNeighbors,
       "Number Of Neighbors", true);
  // ssf.parser_.setArgument("sigma", &sigma,
  //     "Sigma", true);
  // ssf.parser_.setArgument("metaData", &metaData,
  //     "Meta Data", true);

  // ssf.parser_.setArgument("clModulePath", &clModulePath,
  //     "Clustering Module Path", true);
  // ssf.parser_.setArgument("clModuleName", &clModuleName,
  //     "Clustering Module Name", true);
  // ssf.parser_.setArgument("eigenvalues", &eigenvalues,
  //     "Eigenvalues", true);
  // ssf.parser_.setArgument("minimumNumberOfClusters", &minimumNumberOfClusters,
  //     "Minimum Number Of Clusters", true);
  // ssf.parser_.setArgument("MaximumNumberOfClusters", &maximumNumberOfClusters,
  //     "Maximum Number Of Clusters", true);
  // ssf.parser_.setArgument("currentGap", &currentGap,
  //     "Current Gap", true);

  // ssf.parser_.setArgument("computeLowerBound", &computeLowerBound,
  //     "Compute Lower Bound", true);
  // ssf.parser_.setArgument("computeUpperBound", &computeUpperBound,
  //     "Compute Upper Bound", true);
  // ssf.parser_.setArgument("clusteringData", &clusteringData,
  //     "Clustering Data", true);
  // ssf.parser_.setArgument("binCount", &binCount,
  //     "Bin Count", true);

  ssf.parser_.setArgument("mcpThreshold", &mcpThreshold,
      "MCP Persistence threshold", true);

  int ret = 0;
  ret = ssf.init(argc, argv);
  if(ret != 0)
    return ret;

  int debugLevel = ssf.getDebugLevel();
  int threadNumber = ssf.getThreadNumber();

  ret = ssf.run();
  if(ret != 0)
    return ret;

  vtkDataSet* ssfOutput = ssf.getOutput();

  MyProgram<ttkPersistenceMap> pmb;
  ret = pmb.init(argc, argv);
  if(ret != 0)
    return ret;

  pmb.ttkObject_->SetdebugLevel_(debugLevel);
  pmb.ttkObject_->SetThreadNumber(threadNumber);
  pmb.ttkObject_->SetCriticalType(criticalType);
  pmb.ttkObject_->SetPersistenceThreshold(persistenceThreshold);
  pmb.ttkObject_->SetPersistenceScale(persistenceScale);
  pmb.ttkObject_->SetPersistenceSpread(persistenceSpread);
  pmb.ttkObject_->SetDefaultAmplitude(defaultAmplitude);
  pmb.ttkObject_->SetDefaultDeviation(defaultDeviation);
  pmb.ttkObject_->SetReplaceScalarField(true);
  pmb.ttkObject_->SetBatchMode(true);

  pmb.setInput(ssfOutput, 0);

  // execute PersistenceMap algorithm
  ret = pmb.run();
  if(ret != 0)
    return ret;

  vtkDataSet* pmbOutput = pmb.getOutput();

  // save PersistenceMap output
  string pmbFileName="pmb";
  saveVTK(pmbOutput, pmbFileName);

  MyProgram<ttkLDistanceMatrix> ldm;
  ret = ldm.init(argc, argv);
  if(ret != 0)
    return ret;

  ldm.ttkObject_->SetdebugLevel_(debugLevel);
  ldm.ttkObject_->SetThreadNumber(threadNumber);
  ldm.ttkObject_->SetDistanceType(distanceType);
  ldm.ttkObject_->SetTransposeMatrix(transposeMatrix);

  ldm.setInput(pmbOutput, 0);

  if(numberOfClusters>1){
    sigma = numberOfClusters - 1;
  }

  // execute LDistanceMatrix algorithm
  ret = ldm.run();
  if(ret != 0)
    return ret;

  vtkDataSet* ldmOutput = ldm.getOutput();

  // save LDistanceMatrix output
  string ldmFileName="ldm";
  saveVTK(ldmOutput, ldmFileName);

  MyProgram<ttkSpectralEmbedding> se;
  ret = se.init(argc, argv);
  if(ret != 0)
    return ret;

  se.ttkObject_->SetdebugLevel_(debugLevel);
  se.ttkObject_->SetThreadNumber(threadNumber);
  se.ttkObject_->SetModulePath(seModulePath);
  se.ttkObject_->SetModuleName(seModuleName);
  se.ttkObject_->SetFunctionName(functionName);
  se.ttkObject_->SetEmbeddingDimension(embeddingDimension);
  se.ttkObject_->SetMinimumNumberOfComponents(minimumNumberOfComponents);
  se.ttkObject_->SetMaximumNumberOfComponents(maximumNumberOfComponents);
  se.ttkObject_->SetNumberOfNeighbors(numberOfNeighbors);
  se.ttkObject_->SetSigma(sigma);
  se.ttkObject_->SetMetaData(metaData);

  se.setInput(ldmOutput, 0);

  // execute SpectralEmbedding algorithm
  ret = se.run();
  if(ret != 0)
    return ret;

  vtkDataSet* seOutput=se.getOutput();

  MyProgram<ttkClustering> cl;
  ret = cl.init(argc, argv);
  if(ret != 0)
    return ret;

  cl.setInput(seOutput, 0);

  // check super parameter
  if(numberOfClusters>1){
    minimumNumberOfClusters=numberOfClusters;
    maximumNumberOfClusters=numberOfClusters;
  }

  cl.ttkObject_->SetdebugLevel_(debugLevel);
  cl.ttkObject_->SetThreadNumber(threadNumber);
  cl.ttkObject_->SetModulePath(clModulePath);
  cl.ttkObject_->SetModuleName(clModuleName);
  cl.ttkObject_->SetFunctionName(functionName);
  cl.ttkObject_->SetEigenvalues(eigenvalues);
  cl.ttkObject_->SetMinimumNumberOfClusters(minimumNumberOfClusters);
  cl.ttkObject_->SetMaximumNumberOfClusters(maximumNumberOfClusters);
  cl.ttkObject_->SetCurrentGap(currentGap);

  // execute Clustering algorithm
  ret = cl.run();
  if(ret != 0)
    return ret;

  vtkDataSet* clOutput=cl.getOutput();

  MyProgram<ttkSphereFromPoint> sfp;
  ret = sfp.init(argc, argv);
  if(ret != 0)
    return ret;

  sfp.setInput(clOutput, 0);

  sfp.ttkObject_->SetdebugLevel_(debugLevel);
  sfp.ttkObject_->SetThreadNumber(threadNumber);
  sfp.ttkObject_->SetRadius(radius);
  sfp.ttkObject_->SetThetaResolution(thetaResolution);
  sfp.ttkObject_->SetStartTheta(startTheta);
  sfp.ttkObject_->SetEndTheta(endTheta);
  sfp.ttkObject_->SetPhiResolution(phiResolution);
  sfp.ttkObject_->SetStartPhi(startPhi);
  sfp.ttkObject_->SetEndPhi(endPhi);

  // execute SphereFromPoints algorithm
  ret = sfp.run();
  if(ret != 0)
    return ret;

  vtkDataSet* sfpOutput=sfp.getOutput();

  // save SphereFromPoints output
  string sfpFileName="sfp";
  saveVTK(sfpOutput, sfpFileName);

  MyProgram<ttkUncertainDataEstimator> ude;
  ret = ude.init(argc, argv);
  if(ret != 0)
    return ret;

  // use addInput() instead of setInput()
  // because it has multiple inputs.
  ude.addInput(ssfOutput);
  ude.addInput(clOutput);

  ude.ttkObject_->SetdebugLevel_(debugLevel);
  ude.ttkObject_->SetThreadNumber(threadNumber);
  ude.ttkObject_->ComputeLowerBound(computeLowerBound);
  ude.ttkObject_->ComputeUpperBound(computeUpperBound);
  ude.ttkObject_->SetClusteringData(clusteringData);
  ude.ttkObject_->SetMetaData(metaData);
  ude.ttkObject_->SetBinCount(binCount);

  // execute UncertainDataEstimator algorithm:
  // call run2() instead of run() because it has
  // multiple inputs.
  ret = ude.run2();
  if(ret != 0)
    return ret;

  vtkDataSet* udeOutput=ude.getOutput();

  MyProgram<ttkMandatoryCriticalPoints> mcp;
  ret = mcp.init(argc, argv);
  if(ret != 0)
    return ret;

  mcp.setInput(udeOutput, 0);

  mcp.ttkObject_->SetDebugLevel(debugLevel);
  mcp.ttkObject_->SetThreadNumber(threadNumber);
  mcp.ttkObject_->SetClusteringData(clusteringData);
  mcp.ttkObject_->SetMetaData(metaData);
  mcp.ttkObject_->SetSimplificationThreshold(mcpThreshold);
  mcp.ttkObject_->SetOutputMinimumComponentId(false);
  mcp.ttkObject_->SetOutputJoinSaddleComponentId(false);
  mcp.ttkObject_->SetOutputSplitSaddleComponentId(false);
  mcp.ttkObject_->SetOutputMaximumComponentId(false);
  mcp.ttkObject_->setOutputAllMinimumComponents(true);
  mcp.ttkObject_->setOutputAllJoinSaddleComponents(false);
  mcp.ttkObject_->setOutputAllSplitSaddleComponents(false);
  mcp.ttkObject_->setOutputAllMaximumComponents(true);

  // execute MandatoryCriticalPoints algorithm
  ret = mcp.run();
  if(ret != 0)
    return ret;

  vtkDataSet* mandatoryMinima=mcp.getOutput(0);
  vtkDataSet* mandatoryMaxima=mcp.getOutput(3);

  // save MandatoryCriticalPoints output
  string mcpFileName="mcp";
  switch(criticalType){
    case 0:
      saveVTK(mandatoryMinima, mcpFileName+string("_min"));
      break;

    case 1:
      saveVTK(mandatoryMaxima, mcpFileName+string("_max"));
      break;

    case 2:
      saveVTK(mandatoryMinima, mcpFileName+string("_min"));
      saveVTK(mandatoryMaxima, mcpFileName+string("_max"));
      break;
  }

  return 0;
}
