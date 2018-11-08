#include                  <ttkPersistenceMap.h>

vtkStandardNewMacro(ttkPersistenceMap)

int ttkPersistenceMap::doIt(
  vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

  Memory m;
  Timer t;
  
  vtkDataSet *input = inputs[0];
  vtkDataSet *output = outputs[0];
  
  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
 
  if(!triangulation)
    return -1;
  
  triangulation->setWrapper(this);
  persistenceMap_.setupTriangulation(triangulation);
  persistenceMap_.setWrapper(this);
 
  // use a pointer-base copy for the input data -- to adapt if your wrapper does
  // not produce an output of the type of the input.
  output->ShallowCopy(input);

  vtkPointData* inputPointData=input->GetPointData();
  const int numberOfArrays=inputPointData->GetNumberOfArrays();

  for(int i=0; i<numberOfArrays; ++i){
    // in the following, the target scalar field of the input is replaced in the 
    // variable 'output' with the result of the computation.
    // if your wrapper produces an output of the same type of the input, you 
    // should proceed in the same way.
    vtkDataArray *inputScalarField = NULL;

    if(BatchMode){
      inputScalarField = inputPointData->GetArray(i);
    }
    else{
      if(ScalarField.length()){
        inputScalarField = inputPointData->GetArray(ScalarField.data());
      }
      else{
        inputScalarField = inputPointData->GetArray(0);
      }
    }

    {
      stringstream msg;
      msg << "[ttkPersistenceMap] Processing " << inputScalarField->GetName()
       << " scalar field..." << endl;
      dMsg(cout, msg.str(), memoryMsg);
    }

    if(!inputScalarField)
      return -2;

    outputScalarField_ = vtkSmartPointer<vtkDoubleArray>::New();
    outputScalarField_->SetNumberOfTuples(input->GetNumberOfPoints());
    if(ReplaceScalarField)
      outputScalarField_->SetName(inputScalarField->GetName());
    else if(BatchMode)
      outputScalarField_->SetName(string("PersistenceMap_" + string(inputScalarField->GetName())).data());
    else{
      if(ScalarField.length())
        outputScalarField_->SetName(string("PersistenceMap_" + ScalarField).data());
      else
        outputScalarField_->SetName("PersistenceMap");
    }

    // calling the executing package
    switch(inputScalarField->GetDataType()){

      vtkTemplateMacro(
          {
          persistenceMap_.setInputDataPointer(inputScalarField->GetVoidPointer(0));
          persistenceMap_.setOutputDataPointer(outputScalarField_->GetVoidPointer(0));
          persistenceMap_.execute<VTK_TT>(
              PersistenceThreshold,
              PersistenceScale, PersistenceSpread, 
              DefaultAmplitude, DefaultDeviation,
              CriticalType);
          }
          );
    }
    if(ReplaceScalarField)
      output->GetPointData()->RemoveArray(inputScalarField->GetName());
    output->GetPointData()->AddArray(outputScalarField_);

    if(!BatchMode) break;
  }
  
  {
    stringstream msg;
    msg << "[ttkPersistenceMap] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    if(BatchMode)
      msg << "[ttkPersistenceMap] Batchmode processed overall in "
        << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 0;
}
