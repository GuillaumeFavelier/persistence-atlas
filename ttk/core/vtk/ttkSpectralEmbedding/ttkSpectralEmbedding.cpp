#include<ttkSpectralEmbedding.h>

vtkStandardNewMacro(ttkSpectralEmbedding)

int ttkSpectralEmbedding::doIt(vector<vtkDataSet*>& inputs, vector<vtkDataSet*>& outputs){
  Memory m;

  vtkUnstructuredGrid* input=vtkUnstructuredGrid::SafeDownCast(inputs[0]);
  vtkUnstructuredGrid* output=vtkUnstructuredGrid::SafeDownCast(outputs[0]);

  // get matrix values
  vtkDataArray* distanceMatrix=nullptr;
  {
    vtkPointData* pointData=input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData){
      cerr << "[vtkSpectralEmbedding] Error: input has no point data." << endl;
      return -1;
    }
#endif

    if(ScalarField.length()>0)
      distanceMatrix=pointData->GetArray(ScalarField.data());
    else
      distanceMatrix=pointData->GetArray(0);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(!distanceMatrix){
    cerr << "[vtkSpectralEmbedding] Error: no scalar field found." << endl;
    return -1;
  }
#endif

#ifndef TTK_ENABLE_KAMIKAZE
  if(distanceMatrix->GetDataType()!=VTK_DOUBLE){
    cerr << "[vtkSpectralEmbedding] Error: input matrix has unexpected type." << endl;
    return -1;
  }
#endif

#ifndef TTK_ENABLE_KAMIKAZE
  if(MaximumNumberOfComponents<1){
    cerr << "[vtkSpectralEmbedding] Error: maximum number of components is invalid." << endl;
    return -1;
  }
#endif

  const int numberOfPoints=input->GetNumberOfPoints();
#ifndef TTK_ENABLE_KAMIKAZE
  if(numberOfPoints<=0){
    cerr << "[vtkSpectralEmbedding] Error: input has no point." << endl;
    return -1;
  }

  {
    const float tmp=std::sqrt(numberOfPoints);
    const float dec=tmp-std::floor(tmp);
    if(dec>0.0f){
      cerr << "[vtkSpectralEmbedding] Error: input matrix dimensions are not valid." << endl;
      return -1;
    }
  }
#endif
  const int numberOfObjects=std::sqrt(numberOfPoints);
  int numberOfComponents=MaximumNumberOfComponents;

  vector<vtkSmartPointer<vtkDoubleArray>> components(numberOfComponents);
  vector<void*> componentPointers(numberOfComponents);
  for(int i=0; i<numberOfComponents; ++i){
    string s="component_"+std::to_string(i);
    components[i]=vtkSmartPointer<vtkDoubleArray>::New();
    components[i]->SetNumberOfComponents(1);
    components[i]->SetNumberOfTuples(numberOfObjects);
    components[i]->SetName(s.data());
    componentPointers[i]=components[i]->GetVoidPointer(0);
  }

  vtkSmartPointer<vtkDoubleArray> eigenvalues=vtkSmartPointer<vtkDoubleArray>::New();
  eigenvalues->SetNumberOfComponents(1);
  eigenvalues->SetNumberOfTuples(numberOfObjects);
  eigenvalues->SetName("Eigenvalues");

  vtkSmartPointer<vtkDoubleArray> embedding=vtkSmartPointer<vtkDoubleArray>::New();
  embedding->SetNumberOfComponents(3);
  embedding->SetNumberOfTuples(numberOfObjects);
  embedding->SetName("Embedding");

  spectralEmbedding_.setWrapper(this);
  spectralEmbedding_.setInputModulePath(ModulePath);
  spectralEmbedding_.setInputModuleName(ModuleName);
  spectralEmbedding_.setInputFunctionName(FunctionName);
  spectralEmbedding_.setInputMatrixDimension(numberOfObjects);
  spectralEmbedding_.setInputDistanceMatrix(distanceMatrix->GetVoidPointer(0));
  spectralEmbedding_.setInputEmbeddingDimension(EmbeddingDimension);
  spectralEmbedding_.setInputMinimumNumberOfComponents(MinimumNumberOfComponents);
  spectralEmbedding_.setInputMaximumNumberOfComponents(MaximumNumberOfComponents);
  spectralEmbedding_.setInputNumberOfNeighbors(NumberOfNeighbors);
  spectralEmbedding_.setInputSigma(Sigma);
  spectralEmbedding_.setOutputEmbedding(embedding->GetVoidPointer(0));
  spectralEmbedding_.setOutputEigenvalues(eigenvalues->GetVoidPointer(0));
  spectralEmbedding_.setOutputComponents(&componentPointers);
  spectralEmbedding_.execute();

  {
    vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
    points->SetData(embedding);
    output->SetPoints(points);

    vtkPointData* pointData=output->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData){
      cerr << "[vtkSpectralEmbedding] Error: output has no point data." << endl;
      return -1;
    }
#endif

    numberOfComponents=componentPointers.size();
    for(int i=0; i<numberOfComponents; ++i)
      pointData->AddArray(components[i]);
    pointData->AddArray(eigenvalues);

    vtkFieldData* fieldData=input->GetFieldData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!fieldData){
      cerr << "[vtkSpectralEmbedding] Error: input has no field data." << endl;
      return -1;
    }

    if(fieldData->GetNumberOfArrays()<=0){
      cerr << "[vtkSpectralEmbedding] Error: no field data arrays provided." << endl;
      return -1;
    }
#endif

    vtkAbstractArray* metaData=nullptr;
    if(MetaData.length())
      metaData=fieldData->GetAbstractArray(MetaData.data());
    else
      metaData=fieldData->GetAbstractArray(0);
#ifndef TTK_ENABLE_KAMIKAZE
    if(!metaData){
      cerr << "[vtkSpectralEmbedding] Error: no valid field data array provided." << endl;
      return -1;
    }
#endif

    output->GetFieldData()->AddArray(metaData);
  }

  {
    stringstream msg;
    msg << "[ttkSpectralEmbedding] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}
