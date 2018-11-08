#include                  <ttkMandatoryCriticalPoints.h>

vtkStandardNewMacro(ttkMandatoryCriticalPoints)

ttkMandatoryCriticalPoints::ttkMandatoryCriticalPoints(){
  UseAllCores = false;

  computeMinimumOutput_ = false;
  computeJoinSaddleOutput_ = false;
  computeSplitSaddleOutput_ = false;
  computeMaximumOutput_ = false;

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(6);

  simplificationThreshold_ = 0.0;
  simplify_ = true;

  memoryUsage_ = 0.0;

  lowerBoundId = 0;
  upperBoundId = 1;

  outputAllMinimumComponents_ = true;
  outputAllJoinSaddleComponents_ = true;
  outputAllSplitSaddleComponents_ = true;
  outputAllMaximumComponents_ = true;

  computeMinimumOutput_ = true;
  computeJoinSaddleOutput_ = true;
  computeSplitSaddleOutput_ = true;
  computeMaximumOutput_ = true;
}

ttkMandatoryCriticalPoints::~ttkMandatoryCriticalPoints(){
}

int ttkMandatoryCriticalPoints::buildVtkTree(
    MandatoryCriticalPoints &mcp,
    vtkSmartPointer<vtkAppendFilter> app,
    MandatoryCriticalPoints::TreeType treeType,
    int clusterId,
    vtkAbstractArray* clusteringData){
  vtkSmartPointer<vtkUnstructuredGrid> outputTree=vtkSmartPointer<vtkUnstructuredGrid>::New();

  /* Graph Informations */
  const Graph *graph =
    (treeType == MandatoryCriticalPoints::TreeType::JoinTree) ?
    mcp.getJoinTreeGraph()
    : mcp.getSplitTreeGraph();
  const vector<double> *xCoord =
    (treeType == MandatoryCriticalPoints::TreeType::JoinTree) ?
    mcp.getJoinTreeXLayout()
    : mcp.getSplitTreeXLayout();
  const vector<double> *yCoord =
    (treeType == MandatoryCriticalPoints::TreeType::JoinTree) ?
    mcp.getJoinTreeYLayout()
    : mcp.getSplitTreeYLayout();
  const vector<int> *mdtTreePointComponentId =
    (treeType == MandatoryCriticalPoints::TreeType::JoinTree) ?
    mcp.getMdtJoinTreePointComponentId()
    : mcp.getMdtSplitTreePointComponentId();
  const vector<MandatoryCriticalPoints::PointType> *mdtTreePointType =
    (treeType == MandatoryCriticalPoints::TreeType::JoinTree) ?
    mcp.getMdtJoinTreePointType()
    : mcp.getMdtSplitTreePointType();
  const vector<double> *mdtTreePointLowInterval =
    (treeType == MandatoryCriticalPoints::TreeType::JoinTree) ?
    mcp.getMdtJoinTreePointLowInterval()
    : mcp.getMdtSplitTreePointLowInterval();
  const vector<double> *mdtTreePointUpInterval =
    (treeType == MandatoryCriticalPoints::TreeType::JoinTree) ?
    mcp.getMdtJoinTreePointUpInterval()
    : mcp.getMdtSplitTreePointUpInterval();
  const vector<int> *mdtTreeEdgeSwitchable =
    (treeType == MandatoryCriticalPoints::TreeType::JoinTree) ?
    mcp.getMdtJoinTreeEdgeSwitchable()
    : mcp.getMdtSplitTreeEdgeSwitchable();
  const int numberOfPoints = graph->getNumberOfVertices();
  const int numberOfEdges = graph->getNumberOfEdges();

  /* VTK Objects */
  vtkSmartPointer<vtkPoints> mdtTreePoints=vtkSmartPointer<vtkPoints>::New();
  vector<vtkSmartPointer<vtkIdList>> mdtTreeEdge;
  vtkSmartPointer<vtkIntArray> outputTreePointType=vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkDoubleArray> outputTreePointLowInterval=vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> outputTreePointUpInterval=vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> outputTreePointComponentId=vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> outputTreeEdgeSwitchable=vtkSmartPointer<vtkIntArray>::New();

  /* Point datas */
  vtkPointData* outputPointData=outputTree->GetPointData();

  outputTreePointType->SetName("Type");
  outputTreePointType->SetNumberOfTuples(numberOfPoints);
  outputPointData->AddArray(outputTreePointType);

  outputTreePointLowInterval->SetName("LowInterval");
  outputTreePointLowInterval->SetNumberOfTuples(numberOfPoints);
  outputPointData->AddArray(outputTreePointLowInterval);

  outputTreePointUpInterval->SetName("UpInterval");
  outputTreePointUpInterval->SetNumberOfTuples(numberOfPoints);
  outputPointData->AddArray(outputTreePointUpInterval);

  outputTreePointComponentId->SetName("ComponentId");
  outputTreePointComponentId->SetNumberOfTuples(numberOfPoints);
  outputPointData->AddArray(outputTreePointComponentId);

  for(int i=0 ; i<numberOfPoints ; i++){
    outputTreePointType->SetTuple1(i, (*mdtTreePointType)[i]);
    outputTreePointLowInterval->SetTuple1(i, (*mdtTreePointLowInterval)[i]);
    outputTreePointUpInterval->SetTuple1(i, (*mdtTreePointUpInterval)[i]);
    outputTreePointComponentId->SetTuple1(i, (*mdtTreePointComponentId)[i]);
  }

  /* Cell datas */
  vtkCellData* outputCellData=outputTree->GetCellData();

  outputTreeEdgeSwitchable->SetName("Switchable");
  outputTreeEdgeSwitchable->SetNumberOfTuples(numberOfEdges);
  outputCellData->AddArray(outputTreeEdgeSwitchable);

  vtkSmartPointer<vtkIntArray> outputTreeEdgeClusterId=vtkSmartPointer<vtkIntArray>::New();
  if(clusteringData){
    outputTreeEdgeClusterId->SetName("ClusterIdentifier");
    outputTreeEdgeClusterId->SetNumberOfTuples(numberOfEdges);
    outputCellData->AddArray(outputTreeEdgeClusterId);
  }

  for(int i=0 ; i<numberOfEdges ; i++) {
    outputTreeEdgeSwitchable->SetTuple1(i, (*mdtTreeEdgeSwitchable)[i]);

    if(clusteringData){
      outputTreeEdgeClusterId->SetTuple1(i, clusterId);
    }
  }

  /* New Objects */
  // Points
  outputTree->SetPoints(mdtTreePoints);

  // Edges
  mdtTreeEdge.resize(graph->getNumberOfEdges(),NULL);
  for(size_t i=0 ; i<mdtTreeEdge.size() ; i++)
    mdtTreeEdge[i] = vtkSmartPointer<vtkIdList>::New();
  outputTree->Allocate(numberOfEdges);

  // Points
  for(int i=0 ; i<numberOfPoints ; i++) {
    mdtTreePoints->InsertNextPoint((*xCoord)[i], (*yCoord)[i], 0.0);
  }

  // Edges (cells)
  for(int i=0 ; i<numberOfEdges ; i++) {
    mdtTreeEdge[i]->InsertNextId(graph->getEdge(i)->getVertexIdx().first);
    mdtTreeEdge[i]->InsertNextId(graph->getEdge(i)->getVertexIdx().second);
    outputTree->InsertNextCell(VTK_LINE, mdtTreeEdge[i]);
  }

  app->AddInputData(outputTree);
  return 0;
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkMandatoryCriticalPoints::needsToAbort(){
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkMandatoryCriticalPoints::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkMandatoryCriticalPoints] " << progress*100
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkMandatoryCriticalPoints::doIt(vector<vtkDataSet *> &inputs,
    vector<vtkDataSet *> &outputs){
  Timer t;

  const int numInputs=inputs.size();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numInputs<=0){
    cerr << "[ttkMandatoryCriticalPoints] Error: not enough input data is provided." << endl;
    return -1;
  }

  if(numInputs>1){
    cerr << "[ttkMandatoryCriticalPoints] Error: too much input data is provided." << endl;
    return -1;
  }
#endif

  vtkDataSet* input=inputs[0];
  const int numberOfPoints=input->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
    if(numberOfPoints<=0){
      cerr << "[ttkMandatoryCriticalPoints] Error: input has no point." << endl;
      return -1;
    }
#endif

  vtkFieldData* inputFieldData=input->GetFieldData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!inputFieldData){
      cerr << "[ttkMandatoryCriticalPoints] Error: input has no field data." << endl;
      return -1;
    }
#endif
  const bool hasMetaData=(inputFieldData->GetNumberOfArrays()>0);

  vtkPointData* inputPointData=input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!inputPointData){
      cerr << "[ttkMandatoryCriticalPoints] Error: input has no point data." << endl;
      return -1;
    }
#endif

  // read meta data infos
  vtkAbstractArray* clusteringData=nullptr;
  int numberOfObjects=0;
  int numberOfClusters=0;
  if(hasMetaData){
    if(ClusteringData.length())
      clusteringData=inputFieldData->GetAbstractArray(ClusteringData.data());
    else
      clusteringData=inputFieldData->GetAbstractArray(0);
#ifndef TTK_ENABLE_KAMIKAZE
    if(!clusteringData){
      cerr << "[ttkMandatoryCriticalPoints] Error: no valid clustering data provided." << endl;
      return -1;
    }
#endif
    numberOfObjects=clusteringData->GetNumberOfTuples();
#ifndef TTK_ENABLE_KAMIKAZE
    if(numberOfObjects<=0){
      cerr << "[ttkMandatoryCriticalPoints] Error: invalid number of objects: " << numberOfObjects << endl;
      return -1;
    }
#endif

    {
      stringstream msg;
      msg << "[ttkMandatoryCriticalPoints] Number of objects: " << numberOfObjects << endl;
      dMsg(cout, msg.str(), infoMsg);
    }

    int* clusterIds=static_cast<int*>(clusteringData->GetVoidPointer(0));
    numberOfClusters=*std::max_element(clusterIds, clusterIds+numberOfObjects)+1;
#ifndef TTK_ENABLE_KAMIKAZE
    if(numberOfClusters<=0){
      cerr << "[ttkMandatoryCriticalPoints] Error: invalid number of clusters: " << numberOfClusters << endl;
      return -1;
    }
#endif

    {
      stringstream msg;
      msg << "[ttkMandatoryCriticalPoints] Number of clusters: " << numberOfClusters << endl;
      dMsg(cout, msg.str(), infoMsg);
    }
  }

  int numberOfBins=0;
  {
    {
      vtkDataArray* arr=nullptr;
      do{
        string s="probability_bin"+std::to_string(numberOfBins);
        if(clusteringData)
          s="probability_cluster0_bin"+std::to_string(numberOfBins);
        arr=inputPointData->GetArray(s.data());
      }while(arr and ++numberOfBins);
    }
#ifndef TTK_ENABLE_KAMIKAZE
    if(numberOfBins<=0){
      cerr << "[ttkMandatoryCriticalPoints] Error: invalid number of bins: " << numberOfBins << endl;
      return -1;
    }
#endif

    {
      stringstream msg;
      msg << "[ttkMandatoryCriticalPoints] Number of bins: " << numberOfBins << endl;
      dMsg(cout, msg.str(), infoMsg);
    }
  }

  Memory m;
  
  vtkDataSet *outputMinimum = outputs[0];
  vtkDataSet *outputJoinSaddle = outputs[1];
  vtkDataSet *outputSplitSaddle = outputs[2];
  vtkDataSet *outputMaximum = outputs[3];
  vtkUnstructuredGrid *outputJoinTree = 
    vtkUnstructuredGrid::SafeDownCast(outputs[4]);
  vtkUnstructuredGrid *outputSplitTree = 
    vtkUnstructuredGrid::SafeDownCast(outputs[5]);
  
  // Use a pointer-base copy for the input data
  outputMinimum->ShallowCopy(input);
  outputJoinSaddle->ShallowCopy(input);
  outputSplitSaddle->ShallowCopy(input);
  outputMaximum->ShallowCopy(input);

  // Initialize the triangulation object with the input data set
  Triangulation* triangulation = ttkTriangulation::getTriangulation(input);
#ifdef TTK_ENABLE_KAMIKAZE
  if(!triangulation){
    cerr << "[ttkMandatoryCriticalPoints] Error: invalid triangulation. " << endl;
    return -1;
  }
#endif

  triangulation->setWrapper(this);

  const int numberOfOutputs=(clusteringData)?numberOfClusters:1;

  // Allocate the memory for the output scalar field
  using spi=vtkSmartPointer<vtkIntArray>;

  vector<spi> outputMandatoryMinimum;
  vector<spi> outputMandatoryJoinSaddle;
  vector<spi> outputMandatorySplitSaddle;
  vector<spi> outputMandatoryMaximum;

  outputMandatoryMinimum.resize(numberOfOutputs, nullptr);
  outputMandatoryJoinSaddle.resize(numberOfOutputs, nullptr);
  outputMandatorySplitSaddle.resize(numberOfOutputs, nullptr);
  outputMandatoryMaximum.resize(numberOfOutputs, nullptr);
  for(int i=0; i<numberOfOutputs; ++i){
    string strmin="MinimumComponents";
    if(clusteringData)
      strmin="minimumComponents_cluster"+std::to_string(i);

    outputMandatoryMinimum[i]=spi::New();
    outputMandatoryMinimum[i]->SetNumberOfTuples(numberOfPoints);
    outputMandatoryMinimum[i]->SetName(strmin.data());
    outputMinimum->GetPointData()->AddArray(outputMandatoryMinimum[i]);

    string strj="JoinSaddleComponents";
    if(clusteringData)
      strj="joinSaddleComponents_cluster"+std::to_string(i);

    outputMandatoryJoinSaddle[i]=spi::New();
    outputMandatoryJoinSaddle[i]->SetNumberOfTuples(numberOfPoints);
    outputMandatoryJoinSaddle[i]->SetName(strj.data());
    outputJoinSaddle->GetPointData()->AddArray(outputMandatoryJoinSaddle[i]);

    string strs="SplitSaddleComponents";
    if(clusteringData)
      strs="splitSaddleComponents_cluster"+std::to_string(i);

    outputMandatorySplitSaddle[i]=spi::New();
    outputMandatorySplitSaddle[i]->SetNumberOfTuples(numberOfPoints);
    outputMandatorySplitSaddle[i]->SetName(strs.data());
    outputSplitSaddle->GetPointData()->AddArray(outputMandatorySplitSaddle[i]);

    string smax="MaximumComponents";
    if(clusteringData)
      smax="maximumComponents_cluster"+std::to_string(i);

    outputMandatoryMaximum[i]=spi::New();
    outputMandatoryMaximum[i]->SetNumberOfTuples(numberOfPoints);
    outputMandatoryMaximum[i]->SetName(smax.data());
    outputMaximum->GetPointData()->AddArray(outputMandatoryMaximum[i]);
  }

  using spd=vtkSmartPointer<vtkDoubleArray>;

  vector<vector<void*>> inputProbabilityPointers;

  vector<spd> outputMinimumLocationLikelihood;
  vector<void*> outputMinimumLocationLikelihoodPointers;

  vector<spd> outputMaximumLocationLikelihood;
  vector<void*> outputMaximumLocationLikelihoodPointers;

  vector<spd> outputCriticalPointProbability;
  vector<void*> outputCriticalPointProbabilityPointers;

  // reserve memory for input probability bins
  inputProbabilityPointers.resize(numberOfOutputs);

  // reserve memory for output location likelyhood
  outputMinimumLocationLikelihood.resize(numberOfOutputs, nullptr);
  outputMinimumLocationLikelihoodPointers.resize(numberOfOutputs, nullptr);
  outputMaximumLocationLikelihood.resize(numberOfOutputs, nullptr);
  outputMaximumLocationLikelihoodPointers.resize(numberOfOutputs, nullptr);

  for(int i=0; i<numberOfOutputs; ++i){
    // input probability bins
    inputProbabilityPointers[i].resize(numberOfBins, nullptr);
    for(int j=0; j<numberOfBins; ++j){
      string spb="probability_bin"+std::to_string(j);
      if(clusteringData)
        spb="probability_cluster"+std::to_string(i)+"_bin"+std::to_string(j);
      vtkDataArray* arr=inputPointData->GetArray(spb.data());
      inputProbabilityPointers[i][j]=arr->GetVoidPointer(0);
    }

    // output location likelyhood
    string smin="MinimumLocationLikelihood";
    if(clusteringData)
      smin="minimumLocationLikelihood_cluster"+std::to_string(i);
    outputMinimumLocationLikelihood[i]=spd::New();
    outputMinimumLocationLikelihood[i]->SetNumberOfTuples(numberOfPoints);
    outputMinimumLocationLikelihood[i]->SetName(smin.data());
    outputMinimumLocationLikelihoodPointers[i]=outputMinimumLocationLikelihood[i]->GetVoidPointer(0);

    outputMinimum->GetPointData()->AddArray(outputMinimumLocationLikelihood[i]);

    string smax="MaximumLocationLikelihood";
    if(clusteringData)
      smax="maximumLocationLikelihood_cluster"+std::to_string(i);
    outputMaximumLocationLikelihood[i]=spd::New();
    outputMaximumLocationLikelihood[i]->SetNumberOfTuples(numberOfPoints);
    outputMaximumLocationLikelihood[i]->SetName(smax.data());
    outputMaximumLocationLikelihoodPointers[i]=outputMaximumLocationLikelihood[i]->GetVoidPointer(0);

    outputMaximum->GetPointData()->AddArray(outputMaximumLocationLikelihood[i]);
  }

  if(clusteringData){
    // reserve memory for output critical point probability
    outputCriticalPointProbability.resize(numberOfClusters, nullptr);
    outputCriticalPointProbabilityPointers.resize(numberOfClusters, nullptr);

    for(int i=0; i<numberOfClusters; ++i){
      // output critical point probability
      string scp="criticalPointProbability_cluster"+std::to_string(i);
      outputCriticalPointProbability[i]=spd::New();
      outputCriticalPointProbability[i]->SetNumberOfTuples(numberOfPoints);
      outputCriticalPointProbability[i]->SetName(scp.data());
      outputCriticalPointProbabilityPointers[i]=outputCriticalPointProbability[i]->GetVoidPointer(0);

      outputMinimum->GetPointData()->AddArray(outputCriticalPointProbability[i]);
      outputJoinSaddle->GetPointData()->AddArray(outputCriticalPointProbability[i]);
      outputSplitSaddle->GetPointData()->AddArray(outputCriticalPointProbability[i]);
      outputMaximum->GetPointData()->AddArray(outputCriticalPointProbability[i]);
    }
  }

  vtkSmartPointer<vtkAppendFilter> joinTrees=vtkSmartPointer<vtkAppendFilter>::New();
  vtkSmartPointer<vtkAppendFilter> splitTrees=vtkSmartPointer<vtkAppendFilter>::New();

  // preprocess beforehand to avoid race conditions
  MandatoryCriticalPoints::preprocessTriangulation(triangulation);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_) schedule(dynamic)
#endif
  for(int i=0; i<numberOfOutputs; ++i){
    int clusterId=0;

    string upperBoundFieldName=upperBoundFieldName_;
    string lowerBoundFieldName=lowerBoundFieldName_;
    if(clusteringData){
      clusterId=i;
      upperBoundFieldName="upperBound_cluster"+std::to_string(clusterId);
      lowerBoundFieldName="lowerBound_cluster"+std::to_string(clusterId);
    }

    // Input data arrays
    vtkDataArray *inputLowerBoundField = NULL;
    vtkDataArray *inputUpperBoundField = NULL;

    // Get the upper bound field array in the input data set
    if(upperBoundFieldName.length()){
      inputUpperBoundField =
        inputPointData->GetArray(upperBoundFieldName.c_str());
    } else {
      inputUpperBoundField = inputPointData->GetArray(upperBoundId);
    }
#ifdef TTK_ENABLE_KAMIKAZE
    if(!inputUpperBoundField){
      cerr << "[ttkMandatoryCriticalPoints] Error: input upper bound field not found. " << endl;
    }
#endif

    // Get the lower bound field array in the input data set
    if(lowerBoundFieldName.length()){
      inputLowerBoundField =
        inputPointData->GetArray(lowerBoundFieldName.c_str());
    } else {
      inputLowerBoundField = inputPointData->GetArray(lowerBoundId);
    }
#ifdef TTK_ENABLE_KAMIKAZE
    if(!inputLowerBoundField){
      cerr << "[ttkMandatoryCriticalPoints] Error: input lower bound field not found. " << endl;
    }
#endif

    {
      stringstream msg;
      msg << "[ttkMandatoryCriticalPoints] Using `"
        << inputLowerBoundField->GetName() << "' as lower bound..." << endl;
      msg << "[ttkMandatoryCriticalPoints] Using `"
        << inputUpperBoundField->GetName() << "' as upper bound..." << endl;
      dMsg(cout, msg.str(), infoMsg);
    }

    MandatoryCriticalPoints mandatoryCriticalPoints;
    mandatoryCriticalPoints.setupTriangulation(triangulation);
    if(clusteringData)
      mandatoryCriticalPoints.setParallel(false);

    // Calling the executing package
    switch(inputUpperBoundField->GetDataType()){

      vtkTemplateMacro({
          // Reset the baseCode object
          mandatoryCriticalPoints.flush();

          // Wrapper
          mandatoryCriticalPoints.setWrapper(this);
          // Set the number of vertex
          mandatoryCriticalPoints.setVertexNumber(input->GetNumberOfPoints());
          // Set the coordinates of each vertex
          double point[3];
          for (int i = 0; i<numberOfPoints; i++) {
          input->GetPoint(i, point);
          mandatoryCriticalPoints.setVertexPosition(i,point);
          }
          // Set the void pointers to the upper and lower bound fields
          mandatoryCriticalPoints.setLowerBoundFieldPointer(
              inputLowerBoundField->GetVoidPointer(0));
          mandatoryCriticalPoints.setUpperBoundFieldPointer(
              inputUpperBoundField->GetVoidPointer(0));

          // Set the output data pointers
          mandatoryCriticalPoints.setOutputMinimumDataPointer(
              outputMandatoryMinimum[clusterId]->GetVoidPointer(0));
          mandatoryCriticalPoints.setOutputJoinSaddleDataPointer(
              outputMandatoryJoinSaddle[clusterId]->GetVoidPointer(0));
          mandatoryCriticalPoints.setOutputSplitSaddleDataPointer(
              outputMandatorySplitSaddle[clusterId]->GetVoidPointer(0));
          mandatoryCriticalPoints.setOutputMaximumDataPointer(
              outputMandatoryMaximum[clusterId]->GetVoidPointer(0));

          mandatoryCriticalPoints.setClusterId(clusterId);
          mandatoryCriticalPoints.setNumberOfBins(numberOfBins);
          mandatoryCriticalPoints.setInputProbabilityPointer(
                &inputProbabilityPointers[clusterId]);
          mandatoryCriticalPoints.setOutputMaximumLocationLikelihoodPointer(
              outputMaximumLocationLikelihood[clusterId]->GetVoidPointer(0));
          mandatoryCriticalPoints.setOutputMinimumLocationLikelihoodPointer(
              outputMinimumLocationLikelihood[clusterId]->GetVoidPointer(0));
          if(clusteringData){
            mandatoryCriticalPoints.setOutputCriticalPointProbabilityPointer(
                outputCriticalPointProbability[clusterId]->GetVoidPointer(0));
            mandatoryCriticalPoints.setClusteringData(
                clusteringData->GetVoidPointer(0));
            mandatoryCriticalPoints.setNumberOfClusters(numberOfClusters);
            mandatoryCriticalPoints.setNumberOfObjects(numberOfObjects);
          }
          // Set the triangulation object
          // Set the offsets
          mandatoryCriticalPoints.setSoSoffsets();
          // Simplification threshold

          mandatoryCriticalPoints.setSimplificationThreshold(
              simplificationThreshold_);
          // Execute process
          mandatoryCriticalPoints.execute<VTK_TT>();
      });
    }

    // Simplification
    if(simplify_) {
      // Simplify
      mandatoryCriticalPoints.simplifyJoinTree();
      mandatoryCriticalPoints.buildJoinTreePlanarLayout();
      mandatoryCriticalPoints.simplifySplitTree();
      mandatoryCriticalPoints.buildSplitTreePlanarLayout();
    }
    // Outputs
    if(simplify_ or computeMinimumOutput_) {
      if(outputAllMinimumComponents_){
        mandatoryCriticalPoints.outputAllMinima();
        mandatoryCriticalPoints.outputMinimumLocationLikelihood();
      }
      else{
        mandatoryCriticalPoints.outputMinimum(outputMinimumComponentId_);
        mandatoryCriticalPoints.outputMinimumLocationLikelihood(outputMinimumComponentId_);
      }
    }
    if(simplify_ or computeJoinSaddleOutput_) {
      if(outputAllJoinSaddleComponents_)
        mandatoryCriticalPoints.outputAllJoinSaddle();
      else
        mandatoryCriticalPoints.outputJoinSaddle(outputJoinSaddleComponentId_);
    }
    if(simplify_ or computeSplitSaddleOutput_) {
      if(outputAllSplitSaddleComponents_)
        mandatoryCriticalPoints.outputAllSplitSaddle();
      else
        mandatoryCriticalPoints.outputSplitSaddle(outputSplitSaddleComponentId_);
    }
    if(simplify_ or computeMaximumOutput_) {
      if(outputAllMaximumComponents_){
        mandatoryCriticalPoints.outputAllMaxima();
        mandatoryCriticalPoints.outputMaximumLocationLikelihood();
      }
      else{
        mandatoryCriticalPoints.outputMaximum(outputMaximumComponentId_);
        mandatoryCriticalPoints.outputMaximumLocationLikelihood(outputMaximumComponentId_);
      }
    }

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
    {
#endif
    buildVtkTree(mandatoryCriticalPoints,
        joinTrees,
        MandatoryCriticalPoints::TreeType::JoinTree,
        clusterId,
        clusteringData);
    buildVtkTree(mandatoryCriticalPoints,
        splitTrees,
        MandatoryCriticalPoints::TreeType::SplitTree,
        clusterId,
        clusteringData);
#ifdef TTK_ENABLE_OPENMP
    }
#endif
  }

  joinTrees->Update();
  splitTrees->Update();
  outputJoinTree->ShallowCopy(joinTrees->GetOutput());
  outputSplitTree->ShallowCopy(splitTrees->GetOutput());

  computeMinimumOutput_ = false;
  computeJoinSaddleOutput_ = false;
  computeSplitSaddleOutput_ = false;
  computeMaximumOutput_ = false;

  {
    stringstream msg;
    memoryUsage_ += m.getElapsedUsage();
    msg << "[ttkMandatoryCriticalPoints] Memory usage: " << memoryUsage_
        << " MB, elapsed time: " << t.getElapsedTime() << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}

int ttkMandatoryCriticalPoints::FillInputPortInformation(int port, vtkInformation *info){
  if(!this->Superclass::FillInputPortInformation(port, info))
    return 0;
  if(port==0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}

int ttkMandatoryCriticalPoints::FillOutputPortInformation(int port, vtkInformation *info){
  if(!this->Superclass::FillOutputPortInformation(port, info))
    return 0;
  if(port >= 0 && port < 4)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  if((port == 4) || (port == 5))
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  return 1;
}
