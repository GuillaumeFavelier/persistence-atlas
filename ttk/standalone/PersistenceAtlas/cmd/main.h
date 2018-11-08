#pragma once

#include<ttkProgramBase.h>
#include<ttkSortScalarField.h>
#include<ttkPersistenceMap.h>
#include<ttkLDistanceMatrix.h>
#include<ttkSpectralEmbedding.h>
#include<ttkClustering.h>
#include<ttkSphereFromPoint.h>
#include<ttkUncertainDataEstimator.h>
#include<ttkMandatoryCriticalPoints.h>

template<typename ttkModule>
class MyProgram : public vtkProgram<ttkModule>{
  public:
    int getNumberOfInputs(){
      return this->inputs_.size();
    }

    int addInput(vtkDataSet* data){
      this->inputs_.push_back(data);
      this->vtkWrapper_->AddInputData(data);
      return 0;
    }

    vtkDataSet* getInput(const int i){
      return this->inputs_[i];
    }

    int setInput(vtkDataSet* data, const int i){
      this->inputs_[i]=data;
      return 0;
    }

    vtkDataSet* getOutput(){
      return this->vtkWrapper_->GetOutput();
    }

    vtkDataSet* getOutput(int i){
      return this->vtkWrapper_->GetOutput(i);
    }

    int run2(){
      this->vtkWrapper_->Update();
      this->vtkWrapper_->Modified();
      return 0;
    }

    string getOutputPath(){
      return this->outputPath_;
    }

    int getDebugLevel(){
      return this->debugLevel_;
    }

    int getThreadNumber(){
      return this->threadNumber_;
    }
};

template<class vtkWriterClass>
int save(vtkDataSet* data, const string& outputFileName){
  vtkSmartPointer<vtkWriterClass> writer=vtkSmartPointer<vtkWriterClass>::New();
  writer->SetFileName(outputFileName.data());
  writer->SetInputData(data);
  writer->Write();

  return 0;
}

int saveVTK(vtkDataSet* data, const string& outputFileName){
  string fileName=outputFileName.c_str();
  switch(data->GetDataObjectType()){
    case VTK_UNSTRUCTURED_GRID:
    case TTK_UNSTRUCTURED_GRID:
      fileName=fileName+string(".vtu");
      save<vtkXMLUnstructuredGridWriter>(data, fileName);
      break;

    case VTK_IMAGE_DATA:
    case TTK_IMAGE_DATA:
      fileName=fileName+string(".vti");
      save<vtkXMLImageDataWriter>(data, fileName);
      break;

    case VTK_POLY_DATA:
    case TTK_POLY_DATA:
      fileName=fileName+string(".vtp");
      save<vtkXMLPolyDataWriter>(data, fileName);
      break;
  }
  return 0;
}
