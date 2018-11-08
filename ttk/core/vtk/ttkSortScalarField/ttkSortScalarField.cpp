#include                  <ttkSortScalarField.h>

vtkStandardNewMacro(ttkSortScalarField)

  bool ttkSortScalarField::needsToAbort(){
    return GetAbortExecute();
  }

int ttkSortScalarField::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkSortScalarField] " << progress*100 
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkSortScalarField::doIt(vtkDataSet *input, vtkDataSet *output){
  Memory m;

  // point data
  vtkPointData* inputPointData=input->GetPointData();
  vtkSmartPointer<vtkPointData> outputPointData=vtkSmartPointer<vtkPointData>::New();
  {
    const int numberOfArrays=inputPointData->GetNumberOfArrays();
    vector<string> arrays;
    for(int i=0; i<numberOfArrays; ++i){
      vtkDataArray* arr=inputPointData->GetArray(i);
      arrays.push_back(arr->GetName());
    }

    std::sort(arrays.begin(), arrays.end());

    for(int i=0; i<numberOfArrays; ++i){
      vtkDataArray* arr=inputPointData->GetArray(arrays[i].data());
      outputPointData->AddArray(arr);
    }
  }

  // cell data
  vtkCellData* inputCellData=input->GetCellData();
  vtkSmartPointer<vtkCellData> outputCellData=vtkSmartPointer<vtkCellData>::New();
  {
    const int numberOfArrays=inputCellData->GetNumberOfArrays();
    vector<string> arrays;
    for(int i=0; i<numberOfArrays; ++i){
      vtkDataArray* arr=inputCellData->GetArray(i);
      arrays.push_back(arr->GetName());
    }

    std::sort(arrays.begin(), arrays.end());

    for(int i=0; i<numberOfArrays; ++i){
      vtkDataArray* arr=inputCellData->GetArray(arrays[i].data());
      outputCellData->AddArray(arr);
    }
  }

  output->ShallowCopy(input);
  output->GetPointData()->ShallowCopy(outputPointData);
  output->GetCellData()->ShallowCopy(outputCellData);

  {
    stringstream msg;
    msg << "[ttkSortScalarField] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkSortScalarField::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;

  vtkDataSet *input=vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output=vtkDataSet::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkSortScalarField] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}
