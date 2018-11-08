#include                  <ttkUncertainDataEstimator.h>

vtkStandardNewMacro(ttkUncertainDataEstimator)

  ttkUncertainDataEstimator::ttkUncertainDataEstimator(){
    UseAllCores = false;
    computeLowerBound_ = true;
    computeUpperBound_ = true;
    BinCount = 5;

    SetNumberOfInputPorts(1);
    SetNumberOfOutputPorts(1);
  }

ttkUncertainDataEstimator::~ttkUncertainDataEstimator(){
}

bool ttkUncertainDataEstimator::needsToAbort(){
  return GetAbortExecute();
}

int ttkUncertainDataEstimator::updateProgress(const float &progress){
  {
    stringstream msg;
    msg << "[ttkUncertainDataEstimator] " << progress*100
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkUncertainDataEstimator::doIt(vector<vtkDataSet*>& input, vtkDataSet *output){
  Timer t;

  const int numInputs=input.size();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numInputs<=0){
    cerr << "[ttkUncertainDataEstimator] Error: not enough input data is provided." << endl;
    return -1;
  }
#endif

  vtkDataSet* representativeDataSet=nullptr;
  vtkDataSet* clusteringDataSet=nullptr;
  vtkDataArray* clusteringData=nullptr;
  vtkAbstractArray* metaData=nullptr;
  int numberOfClusters=0;

  // Get arrays from input data
  vector<vtkDataArray*> inputScalarFields;
  if(numInputs>2){
    inputScalarFields.resize(numInputs, nullptr);
    for(int i=0; i<numInputs; i++){
#ifndef TTK_ENABLE_KAMIKAZE
      if(!input[i]){
        cerr << "[ttkUncertainDataEstimator] Error: input dataset [" << i << "] is null." << endl;
        return -1;
      }
#endif

      vtkPointData* pointData=input[i]->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
      if(!pointData){
        cerr << "[ttkUncertainDataEstimator] Error: input dataset ["<< i << "] has no point data." << endl;
        return -1;
      }
#endif

      if(ScalarField.length())
        inputScalarFields[i]=pointData->GetArray(ScalarField.data());
      else
        inputScalarFields[i]=pointData->GetArray(0);
    }

    representativeDataSet=input[0];
  }
  else if(numInputs>0){
    if(numInputs==1)
      representativeDataSet=input[0];
    else{
      vtkPointData* pointData=nullptr;
      vtkPointData* pointData0=input[0]->GetPointData();
      vtkPointData* pointData1=input[1]->GetPointData();

      vtkFieldData* fieldData=nullptr;
      vtkFieldData* fieldData0=input[0]->GetFieldData();
      vtkFieldData* fieldData1=input[1]->GetFieldData();

#ifndef TTK_ENABLE_KAMIKAZE
      if(!pointData0 or !pointData1){
        cerr << "[ttkUncertainDataEstimator] Error: no point data information provided." << endl;
        return -1;
      }
#endif

#ifndef TTK_ENABLE_KAMIKAZE
      if(!fieldData0 or !fieldData1){
        cerr << "[ttkUncertainDataEstimator] Error: no field data information provided." << endl;
        return -1;
      }
#endif

      if(fieldData0->GetNumberOfArrays() and !fieldData1->GetNumberOfArrays()){
        representativeDataSet=input[1];
        clusteringDataSet=input[0];
        fieldData=fieldData0;
        pointData=pointData0;
      }
      else if(!fieldData0->GetNumberOfArrays() and fieldData1->GetNumberOfArrays()){
        representativeDataSet=input[0];
        clusteringDataSet=input[1];
        fieldData=fieldData1;
        pointData=pointData1;
      }
#ifndef TTK_ENABLE_KAMIKAZE
      else{
        cerr << "[ttkUncertainDataEstimator] Error: field data does not disambiguate representative dataset." << endl;
        return -1;
      }
#endif

      if(ClusteringData.length())
        clusteringData=pointData->GetArray(ClusteringData.data());
      else
        clusteringData=pointData->GetArray(0);
#ifndef TTK_ENABLE_KAMIKAZE
      if(!clusteringData){
        cerr << "[ttkUncertainDataEstimator] Error: no valid clustering data provided." << endl;
        return -1;
      }
#endif

      if(MetaData.length())
        metaData=fieldData->GetAbstractArray(MetaData.data());
      else
        metaData=fieldData->GetAbstractArray(0);
#ifndef TTK_ENABLE_KAMIKAZE
      if(!metaData){
        cerr << "[ttkUncertainDataEstimator] Error: no valid meta data provided." << endl;
        return -1;
      }
#endif

#ifndef TTK_ENABLE_KAMIKAZE
      if(clusteringData->GetNumberOfTuples() != metaData->GetNumberOfTuples()){
        cerr << "[ttkUncertainDataEstimator] Error: inconsistency between clustering and meta data." << endl;
        return -1;
      }
#endif

      const int numberOfElements=clusteringData->GetNumberOfTuples();
      int* clusterIds=static_cast<int*>(clusteringData->GetVoidPointer(0));
      numberOfClusters=*std::max_element(clusterIds, clusterIds+numberOfElements)+1;

#ifndef TTK_ENABLE_KAMIKAZE
      if(numberOfClusters<=0){
        cerr << "[ttkUncertainDataEstimator] Error: invalid number of clusters: " << numberOfClusters << endl;
        return -1;
      }
#endif

      stringstream msg;
      msg << "[ttkUncertainDataEstimator] Number of clusters: " << numberOfClusters << endl;
      dMsg(cout, msg.str(), infoMsg);
    }

    vtkPointData* pointData=representativeDataSet->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData){
      cerr << "[ttkUncertainDataEstimator] Error: representative dataset has no point data." << endl;
      return -1;
    }
#endif

    const int numberOfArrays=pointData->GetNumberOfArrays();
    inputScalarFields.resize(numberOfArrays, nullptr);
    for(int i=0; i<numberOfArrays; i++)
      inputScalarFields[i]=pointData->GetArray(i);
  }

  const int numberOfPoints=representativeDataSet->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
    if(numberOfPoints<=0){
      cerr << "[ttkUncertainDataEstimator] Error: representative dataset has no point." << endl;
      return -1;
    }
#endif

  // Use a pointer-base copy for the input data
  output->ShallowCopy(representativeDataSet);

  vtkDataArray* representativeDataArray=nullptr;
  const int numberOfScalarFields=inputScalarFields.size();
  vector<void*> inputScalarFieldPointers(numberOfScalarFields, nullptr);
  {
    stringstream msg;
    msg << "[ttkUncertainDataEstimator] Number of scalar fields: " << numberOfScalarFields << endl;
    dMsg(cout, msg.str(), infoMsg);

#ifndef TTK_ENABLE_KAMIKAZE
    if(numberOfScalarFields<=0){
      cerr << "[ttkUncertainDataEstimator] Error: no valid input scalar field provided." << endl;
      return -1;
    }
#endif

    representativeDataArray=inputScalarFields[0];

    for(int i=0; i<numberOfScalarFields; i++){
#ifndef TTK_ENABLE_KAMIKAZE
      // Check if all the inputs are here
      if(!inputScalarFields[i]){
        cerr << "[ttkUncertainDataEstimator] Error: input scalar field ["<< i<<"] is null." << endl;
        return -1;
      }
#endif
      inputScalarFieldPointers[i]=inputScalarFields[i]->GetVoidPointer(0);

#ifndef TTK_ENABLE_KAMIKAZE
      // Check if inputs have the same data type and the same number of points
      if(inputScalarFields[i]->GetDataType() != representativeDataArray->GetDataType()){
        stringstream msg;
        msg << "[ttkUncertainDataEstimator] Error: inputs of different data types." << endl;
        dMsg(cerr, msg.str(), fatalMsg);
        return -1;
      }
      if(inputScalarFields[i]->GetNumberOfTuples() != representativeDataArray->GetNumberOfTuples()){
        stringstream msg;
        msg << "[ttkUncertainDataEstimator] Error: inputs with different number of points." << endl;
        dMsg(cerr, msg.str(), fatalMsg);
        return -1;
      }
#endif
    }
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(BinCount<=0){
    cerr << "[ttkUncertainDataEstimator] Error: input number of bins is invalid." << endl;
    return -1;
  }
#endif

  // Allocate the memory for the output scalar fields
  vector<vtkDataArray*> outputLowerBoundScalarFields;
  vector<vtkDataArray*> outputUpperBoundScalarFields;
  vector<vtkDataArray*> outputMeanScalarFields;
  vector<vector<vtkDataArray*>> outputProbabilityFields;
  vector<void*> outputLowerBoundScalarFieldPointers;
  vector<void*> outputUpperBoundScalarFieldPointers;
  vector<void*> outputMeanScalarFieldPointers;
  vector<vector<void*>> outputProbabilityFieldPointers;

  vtkSmartPointer<vtkPointData> outputPointData=vtkSmartPointer<vtkPointData>::New();

  if(clusteringData and metaData){
    if(computeLowerBound_){
      outputLowerBoundScalarFields.resize(numberOfClusters, nullptr);
      outputLowerBoundScalarFieldPointers.resize(numberOfClusters, nullptr);
    }

    if(computeUpperBound_){
      outputUpperBoundScalarFields.resize(numberOfClusters, nullptr);
      outputUpperBoundScalarFieldPointers.resize(numberOfClusters, nullptr);
    }

    if(computeLowerBound_ and computeUpperBound_){
      outputProbabilityFields.resize(numberOfClusters);
      outputProbabilityFieldPointers.resize(numberOfClusters);
    }

    outputMeanScalarFields.resize(numberOfClusters, nullptr);
    outputMeanScalarFieldPointers.resize(numberOfClusters, nullptr);

    for(int i=0; i<numberOfClusters; ++i){
      const string s=std::to_string(i);

      if(computeLowerBound_){
        const string lowerBound="lowerBound_cluster"+s;
        outputLowerBoundScalarFields[i]=vtkDoubleArray::New();
        outputLowerBoundScalarFields[i]->SetNumberOfTuples(numberOfPoints);
        outputLowerBoundScalarFields[i]->SetName(lowerBound.data());
        outputLowerBoundScalarFieldPointers[i]=outputLowerBoundScalarFields[i]->GetVoidPointer(0);
        outputPointData->AddArray(outputLowerBoundScalarFields[i]);
      }

      if(computeUpperBound_){
        const string upperBound="upperBound_cluster"+s;
        outputUpperBoundScalarFields[i]=vtkDoubleArray::New();
        outputUpperBoundScalarFields[i]->SetNumberOfTuples(numberOfPoints);
        outputUpperBoundScalarFields[i]->SetName(upperBound.data());
        outputUpperBoundScalarFieldPointers[i]=outputUpperBoundScalarFields[i]->GetVoidPointer(0);
        outputPointData->AddArray(outputUpperBoundScalarFields[i]);
      }

      if(computeLowerBound_ and computeUpperBound_){
        outputProbabilityFields[i].resize(BinCount, nullptr);
        outputProbabilityFieldPointers[i].resize(BinCount, nullptr);
        for(int j=0; j<BinCount; ++j){
          const string t=std::to_string(j);
          const string probability="probability_cluster"+s+"_bin"+t;

          outputProbabilityFields[i][j]=vtkDoubleArray::New();
          outputProbabilityFields[i][j]->SetNumberOfTuples(numberOfPoints);
          outputProbabilityFields[i][j]->SetName(probability.data());
          outputProbabilityFieldPointers[i][j]=outputProbabilityFields[i][j]->GetVoidPointer(0);
          outputPointData->AddArray(outputProbabilityFields[i][j]);
        }
      }

      const string mean="mean_cluster"+s;
      outputMeanScalarFields[i]=vtkDoubleArray::New();
      outputMeanScalarFields[i]->SetNumberOfTuples(numberOfPoints);
      outputMeanScalarFields[i]->SetName(mean.data());
      outputMeanScalarFieldPointers[i]=outputMeanScalarFields[i]->GetVoidPointer(0);
      outputPointData->AddArray(outputMeanScalarFields[i]);
    }
  }
  else{
    if(computeLowerBound_){
      outputLowerBoundScalarFields.resize(1, nullptr);
      outputLowerBoundScalarFieldPointers.resize(1, nullptr);
      outputLowerBoundScalarFields[0]=vtkDoubleArray::New();
      outputLowerBoundScalarFields[0]->SetNumberOfTuples(numberOfPoints);
      outputLowerBoundScalarFields[0]->SetName("LowerBoundField");
      outputLowerBoundScalarFieldPointers[0]=outputLowerBoundScalarFields[0]->GetVoidPointer(0);
      outputPointData->AddArray(outputLowerBoundScalarFields[0]);
    }

    if(computeUpperBound_){
      outputUpperBoundScalarFields.resize(1, nullptr);
      outputUpperBoundScalarFieldPointers.resize(1, nullptr);
      outputUpperBoundScalarFields[0]=vtkDoubleArray::New();
      outputUpperBoundScalarFields[0]->SetNumberOfTuples(numberOfPoints);
      outputUpperBoundScalarFields[0]->SetName("UpperBoundField");
      outputUpperBoundScalarFieldPointers[0]=outputUpperBoundScalarFields[0]->GetVoidPointer(0);
      outputPointData->AddArray(outputUpperBoundScalarFields[0]);
    }

    if(computeLowerBound_ and computeUpperBound_){
      outputProbabilityFields.resize(1);
      outputProbabilityFieldPointers.resize(1);
      outputProbabilityFields[0].resize(BinCount, nullptr);
      outputProbabilityFieldPointers[0].resize(BinCount, nullptr);
      for(int j=0; j<BinCount; ++j){
        const string t=std::to_string(j);
        const string probability="probability_bin"+t;

        outputProbabilityFields[0][j]=vtkDoubleArray::New();
        outputProbabilityFields[0][j]->SetNumberOfTuples(numberOfPoints);
        outputProbabilityFields[0][j]->SetName(probability.data());
        outputProbabilityFieldPointers[0][j]=outputProbabilityFields[0][j]->GetVoidPointer(0);
        outputPointData->AddArray(outputProbabilityFields[0][j]);
      }
    }

    outputMeanScalarFields.resize(1, nullptr);
    outputMeanScalarFieldPointers.resize(1, nullptr);
    outputMeanScalarFields[0]=vtkDoubleArray::New();
    outputMeanScalarFields[0]->SetNumberOfTuples(numberOfPoints);
    outputMeanScalarFields[0]->SetName("MeanField");
    outputMeanScalarFieldPointers[0]=outputMeanScalarFields[0]->GetVoidPointer(0);
    outputPointData->AddArray(outputMeanScalarFields[0]);
  }

  // Calling the executing package
  {
    UncertainDataEstimator uncertainDataEstimator;
    uncertainDataEstimator.setWrapper(this);

    uncertainDataEstimator.setVertexNumber(numberOfPoints);
    uncertainDataEstimator.setBinCount(BinCount);

    uncertainDataEstimator.setNumberOfInputFields(numberOfScalarFields);
    uncertainDataEstimator.setInputScalarFields(&inputScalarFieldPointers);

    uncertainDataEstimator.setComputeLowerBound(computeLowerBound_);
    uncertainDataEstimator.setComputeUpperBound(computeUpperBound_);

    uncertainDataEstimator.setNumberOfClusters(numberOfClusters);
    if(clusteringData)
      uncertainDataEstimator.setClusteringData(clusteringData->GetVoidPointer(0));
    uncertainDataEstimator.setOutputLowerBoundFields(&outputLowerBoundScalarFieldPointers);
    uncertainDataEstimator.setOutputUpperBoundFields(&outputUpperBoundScalarFieldPointers);
    uncertainDataEstimator.setOutputMeanFields(&outputMeanScalarFieldPointers);
    uncertainDataEstimator.setOutputProbabilityFields(&outputProbabilityFieldPointers);

    switch(representativeDataArray->GetDataType()){
      vtkTemplateMacro((
            {
            uncertainDataEstimator.execute<VTK_TT>();
            }
            ));
    }
  }

  output->GetPointData()->ShallowCopy(outputPointData);

  if(clusteringData and metaData)
    output->GetFieldData()->ShallowCopy(clusteringDataSet->GetFieldData());

  return 0;
}

int ttkUncertainDataEstimator::FillInputPortInformation(int port, vtkInformation *info){
  if(!this->Superclass::FillInputPortInformation(port, info)){
    return 0;
  }
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}

int ttkUncertainDataEstimator::FillOutputPortInformation(int port, vtkInformation *info){
  if(!this->Superclass::FillOutputPortInformation(port, info)){
    return 0;
  }
  if(port==0 || port==1 || port==3)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}

int ttkUncertainDataEstimator::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector, vtkInformationVector *outputVector){
  Memory m;

  // Output pointers informations
  vtkInformation* outInfo;
  // Unified bound fields
  outInfo = outputVector->GetInformationObject(0);
  vtkDataSet *output = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Number of input files
  const int numInputs = inputVector[0]->GetNumberOfInformationObjects();
  {
    stringstream msg;
    msg << "[ttkUncertainDataEstimator] Number of inputs: " << numInputs << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  // Get input datas
  vector<vtkDataSet*> input(numInputs, nullptr);
  for(int i=0 ; i<numInputs ; i++)
    input[i] = vtkDataSet::GetData(inputVector[0], i);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkUncertainDataEstimator] Memory usage: " << m.getElapsedUsage()
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
