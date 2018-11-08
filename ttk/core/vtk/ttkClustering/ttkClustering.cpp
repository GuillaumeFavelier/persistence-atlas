#include<ttkClustering.h>

vtkStandardNewMacro(ttkClustering)

  bool ttkClustering::needsToAbort(){
    return GetAbortExecute();
  }

int ttkClustering::updateProgress(const float &progress){
  {
    stringstream msg;
    msg << "[ttkClustering] " << progress*100
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkClustering::doIt(vtkDataSet* input, vtkDataSet* output){
  Memory m;

  const int numberOfObjects=input->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfObjects){
    cerr << "[ttkClustering] Error: input has no point." << endl;
    return -1;
  }
#endif

  vtkSmartPointer<vtkIntArray> assignation=vtkSmartPointer<vtkIntArray>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!assignation){
    cerr << "[ttkClustering] Error: vtkIntArray memory allocation problem." << endl;
    return -1;
  }
#endif
  assignation->SetNumberOfComponents(1);
  assignation->SetNumberOfTuples(numberOfObjects);
  assignation->SetName("ClusterIdentifier");

  vtkPointData* inputPointData=input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputPointData){
    cerr << "[ttkClustering] Error: input has no point data." << endl;
    return -1;
  }
#endif

  const int numberOfArrays=inputPointData->GetNumberOfArrays();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfArrays){
    cerr << "[ttkClustering] Error: no point data arrays found." << endl;
    return -1;
  }
#endif

  int numberOfComponents=0;
  {
    vtkDataArray* arr=nullptr;
    do{
      string s="component_"+std::to_string(numberOfComponents);
      arr=inputPointData->GetArray(s.data());
    }while(arr and ++numberOfComponents);
  }

  if(numberOfComponents>=numberOfObjects-1)
    numberOfComponents=numberOfObjects-2;

  vector<double> maps(numberOfObjects*numberOfComponents);
  for(int i=0; i<numberOfComponents; ++i){
    vtkDataArray* arr=inputPointData->GetArray(i);
#ifndef TTK_ENABLE_KAMIKAZE
    if(!arr){
      cerr << "[ttkClustering] Error: input point data array [" << i << "] is null." << endl;
      return -1;
    }
#endif

    double* data=static_cast<double*>(arr->GetVoidPointer(0));

    for(int j=0; j<numberOfObjects; ++j)
      maps[j*numberOfComponents+i]=data[j];
  }

  vtkDataArray* eigenvalues=inputPointData->GetArray(Eigenvalues.data());
#ifndef TTK_ENABLE_KAMIKAZE
  if(!eigenvalues){
    cerr << "[ttkClustering] Error: no eigenvalues provided." << endl;
    return -1;
  }
#endif

  clustering_.setWrapper(this);
  clustering_.setInputModulePath(ModulePath);
  clustering_.setInputModuleName(ModuleName);
  clustering_.setInputFunctionName(FunctionName);
  clustering_.setInputNumberOfObjects(numberOfObjects);
  clustering_.setInputMaximumNumberOfClusters(MaximumNumberOfClusters);
  clustering_.setInputMinimumNumberOfClusters(MinimumNumberOfClusters);
  clustering_.setInputCurrentGap(CurrentGap);
  clustering_.setInputNumberOfComponents(numberOfComponents);
  clustering_.setInputMaps(&maps);
  clustering_.setInputEigenvalues(eigenvalues->GetVoidPointer(0));
  clustering_.setOutputAssignation(assignation->GetVoidPointer(0));
  clustering_.execute();

  output->ShallowCopy(input);

  vtkPointData* outputPointData=output->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputPointData){
    cerr << "[ttkClustering] Error: output has no point data." << endl;
    return -1;
  }
#endif

  outputPointData->AddArray(assignation);

  vtkFieldData* outputFieldData=output->GetFieldData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputFieldData){
    cerr << "[ttkLDistanceMatrix] Error: output vtkImageData has no field data." << endl;
    return -1;
  }
#endif
 outputFieldData->AddArray(assignation);

  {
    stringstream msg;
    msg << "[ttkClustering] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkClustering::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector){
  vtkDataSet* input=vtkDataSet::GetData(inputVector[0]);
  vtkInformation* outInfo=outputVector->GetInformationObject(0);
  vtkDataSet* output=vtkDataSet::GetData(outInfo);

  doIt(input, output);

  return 1;
}
