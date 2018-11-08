#include                  <ttkSphereFromPoint.h>

vtkStandardNewMacro(ttkSphereFromPoint)

ttkSphereFromPoint::ttkSphereFromPoint(){

  ThetaResolution = 20;
  StartTheta = 0;
  EndTheta = 360;
  
  PhiResolution = 20;
  StartPhi = 0;
  EndPhi = 180;
  
  Radius = 0.5;
 
  masterAppender_ = NULL;
}

ttkSphereFromPoint::~ttkSphereFromPoint(){

  if(masterAppender_)
    masterAppender_->Delete();

  for(int i = 0; i < (int) appenderList_.size(); i++){
    appenderList_[i]->Delete();
  }
  
  for(int i = 0; i < (int) sphereList_.size(); i++){
    sphereList_[i]->Delete();
  }
 
  for(int i = 0; i < (int) dataArrayList_.size(); i++){
    for(int j = 0; j < (int) dataArrayList_[i].size(); j++){
      dataArrayList_[i][j]->Delete();
    }
  }
}


// transmit abort signals -- to copy paste in other wrappers
bool ttkSphereFromPoint::needsToAbort(){
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkSphereFromPoint::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[ttkSphereFromPoint] " << progress*100 
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }
  
  UpdateProgress(progress);
  return 0;
}

int ttkSphereFromPoint::doIt(vtkDataSet *input, vtkPolyData *output){

  Timer t;
  
  if(masterAppender_){
    masterAppender_->Delete();
    masterAppender_ = NULL;
  }

  for(int i = 0; i < (int) appenderList_.size(); i++){
    appenderList_[i]->Delete();
  }
 
  masterAppender_ = vtkAppendPolyData::New();
  appenderList_.resize(threadNumber_);
  for(int i = 0; i < (int) appenderList_.size(); i++){
    appenderList_[i] = vtkAppendPolyData::New();
  }
  
  if((int) sphereList_.size() > input->GetNumberOfPoints()){
    for(int i = input->GetNumberOfPoints(); i < (int) sphereList_.size(); i++){
      sphereList_[i]->Delete();
    }
    
    sphereList_.resize(input->GetNumberOfPoints());
  }
  else if((int) sphereList_.size() < input->GetNumberOfPoints()){
    int oldSize = sphereList_.size();
    
    sphereList_.resize(input->GetNumberOfPoints());
    
    for(int i = oldSize; i < (int) sphereList_.size(); i++){
      sphereList_[i] = vtkSphereSource::New();
    }
  }
  
  if(dataArrayList_.size()){
    for(int i = 0; i < (int) dataArrayList_.size(); i++){
      for(int j = 0; j < (int) dataArrayList_[i].size(); j++){
        dataArrayList_[i][j]->Delete();
      }
      dataArrayList_[i].clear();
    }
  }
  dataArrayList_.resize(threadNumber_);
  
  vector<vector<double> > p(threadNumber_);
  
  for(int i = 0; i < (int) p.size(); i++){
    p[i].resize(3);
  }

  if(!input->GetNumberOfPoints())
    return -1;
  // make this function thread-safe
  input->GetPoint(0, p[0].data());
  
  int count = 0;
  
#ifdef TTK_ENABLE_OPENMP
  omp_lock_t writeLock;
  omp_init_lock(&writeLock);
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < input->GetNumberOfPoints(); i++){
    
    if(!needsToAbort()){
    
      int threadId = 0;
    
#ifdef TTK_ENABLE_OPENMP
      threadId = omp_get_thread_num();
#endif
    
      input->GetPoint(i, p[threadId].data());
      
      sphereList_[i]->SetCenter(
        p[threadId][0], p[threadId][1], p[threadId][2]);
      
      sphereList_[i]->SetRadius(Radius);
      
      sphereList_[i]->SetThetaResolution(ThetaResolution);
      sphereList_[i]->SetStartTheta(StartTheta);
      sphereList_[i]->SetEndTheta(EndTheta);

      sphereList_[i]->SetPhiResolution(PhiResolution);
      sphereList_[i]->SetStartPhi(StartPhi);
      sphereList_[i]->SetEndPhi(EndPhi); 
      
      sphereList_[i]->Update();
    
      vtkPolyData *sphereSurface = sphereList_[i]->GetOutput();
      
      // copy the data values
      for(int j = 0; j < input->GetPointData()->GetNumberOfArrays(); j++){
        
        vtkDataArray *array = input->GetPointData()->GetArray(j);
        
        if(array->GetNumberOfComponents() == 1){
          
          double value = 0;
          
          switch(array->GetDataType()){

            case VTK_CHAR:
              {
                array->GetTuple(i, &value);
                
                vtkCharArray *charArray = vtkCharArray::New();
                
                charArray->SetName(array->GetName());
                charArray->SetNumberOfTuples(
                  sphereSurface->GetNumberOfPoints());
                for(int k = 0; k < sphereSurface->GetNumberOfPoints(); k++){
                  charArray->SetTuple(k, &value);
                }
                sphereSurface->GetPointData()->AddArray(charArray);
                dataArrayList_[threadId].push_back(charArray);
              }
              break;
            
            case VTK_DOUBLE:
              {
                array->GetTuple(i, &value);
                
                vtkDoubleArray *doubleArray = vtkDoubleArray::New();
                
                doubleArray->SetName(array->GetName());
                doubleArray->SetNumberOfTuples(
                  sphereSurface->GetNumberOfPoints());
                for(int k = 0; k < sphereSurface->GetNumberOfPoints(); k++){
                  doubleArray->SetTuple(k, &value);
                }
                sphereSurface->GetPointData()->AddArray(doubleArray);
                dataArrayList_[threadId].push_back(doubleArray);
              }
              break;
            
            case VTK_FLOAT:
              {
                array->GetTuple(i, &value);
                
                vtkFloatArray *floatArray = vtkFloatArray::New();
                
                floatArray->SetName(array->GetName());
                floatArray->SetNumberOfTuples(
                  sphereSurface->GetNumberOfPoints());
                for(int k = 0; k < sphereSurface->GetNumberOfPoints(); k++){
                  floatArray->SetTuple(k, &value);
                }
                sphereSurface->GetPointData()->AddArray(floatArray);
                dataArrayList_[threadId].push_back(floatArray);
              }
              break;

            case VTK_INT:
              {
                array->GetTuple(i, &value);
                
                vtkIntArray *intArray = vtkIntArray::New();
                
                intArray->SetName(array->GetName());
                intArray->SetNumberOfTuples(
                  sphereSurface->GetNumberOfPoints());
                for(int k = 0; k < sphereSurface->GetNumberOfPoints(); k++){
                  intArray->SetTuple(k, &value);
                }
                sphereSurface->GetPointData()->AddArray(intArray);
                dataArrayList_[threadId].push_back(intArray);
              }
              break;
              
            default:
              {
                stringstream msg;
                msg << "[ttkSphereFromPoint] Unsupported data type :(" << endl;
                dMsg(cerr, msg.str(), detailedInfoMsg);
              }
              break;
          }
        }
        else{
          stringstream msg;
          msg << "[ttkSphereFromPoint] Unsupported number of components :(" 
            << endl;
          dMsg(cerr, msg.str(), detailedInfoMsg);
        }
      }
     
      appenderList_[threadId]->AddInputConnection(
        sphereList_[i]->GetOutputPort());
 
      if(debugLevel_ > 3){
#ifdef TTK_ENABLE_OPENMP
        omp_set_lock(&writeLock);
#endif   
        if(!(count % (input->GetNumberOfPoints()/10))){
          updateProgress((count + 1.0)/input->GetNumberOfPoints());
        }

        count++;
#ifdef TTK_ENABLE_OPENMP
        omp_unset_lock(&writeLock);
#endif
      }
    }
  }
 
  for(int i = 0; i < (int) appenderList_.size(); i++){
    if(appenderList_[i]->GetInput())
      masterAppender_->AddInputConnection(
        appenderList_[i]->GetOutputPort());
  }
  masterAppender_->Update();
  
#ifdef TTK_ENABLE_OPENMP
  omp_destroy_lock(&writeLock);
#endif
  
  output->ShallowCopy(masterAppender_->GetOutput());

  {
    stringstream msg;
    msg << "[ttkSphereFromPoint] Spheres computed in " 
      << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkSphereFromPoint::RequestData(vtkInformation *request, 
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;
  
  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);
  
  doIt(input, output);
  
  {
    stringstream msg;
    msg << "[ttkSphereFromPoint] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 1;
}