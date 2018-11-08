/// \ingroup vtk
/// \class ttkClustering
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the clustering processing package.
///
/// VTK wrapping code for the @Clustering package.
/// 
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa ttk::Clustering
#pragma once

// ttk code includes
#include<Clustering.h>
#include<ttkWrapper.h>

// VTK includes
#include<vtkLongLongArray.h>
#include<vtkDataArray.h>
#include<vtkDataSet.h>
#include<vtkDataSetAlgorithm.h>
#include<vtkDoubleArray.h>
#include<vtkFiltersCoreModule.h>
#include<vtkFloatArray.h>
#include<vtkInformation.h>
#include<vtkIntArray.h>
#include<vtkObjectFactory.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>
#include<vtkTable.h>
#include<vtkPoints.h>

class VTKFILTERSCORE_EXPORT ttkClustering 
: public vtkDataSetAlgorithm, public Wrapper{

  public:

    static ttkClustering* New();
    vtkTypeMacro(ttkClustering, vtkDataSetAlgorithm)

      // default ttk setters
      vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = OsCall::getNumberOfCores();
      }
      Modified();
    }

    void SetThreadNumber(int threadNumber){\
      ThreadNumber = threadNumber;\
        SetThreads();\
    }\
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    vtkSetMacro(ModulePath, string);
    vtkGetMacro(ModulePath, string);

    vtkSetMacro(ModuleName, string);
    vtkGetMacro(ModuleName, string);

    vtkSetMacro(FunctionName, string);
    vtkGetMacro(FunctionName, string);

    vtkSetMacro(Eigenvalues, string);
    vtkGetMacro(Eigenvalues, string);

    vtkSetMacro(MaximumNumberOfClusters, int);
    vtkGetMacro(MaximumNumberOfClusters, int);

    vtkSetMacro(MinimumNumberOfClusters, int);
    vtkGetMacro(MinimumNumberOfClusters, int);

    vtkSetMacro(CurrentGap, int);
    vtkGetMacro(CurrentGap, int);

    int FillInputPortInformation(int port, vtkInformation *info){

      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
          break;
        default:
          break;
      }

      return 1;
    }

    int FillOutputPortInformation(int port, vtkInformation *info){

      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
          break;
        default:
          break;
      }

      return 1;
    }

  protected:

    ttkClustering(){
      UseAllCores=false;
      ThreadNumber=1;
      MaximumNumberOfClusters=2;
      MinimumNumberOfClusters=2;
      CurrentGap=1;
      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~ttkClustering(){};

    int RequestData(vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector);

  private:

    bool UseAllCores;
    int ThreadNumber;
    string ModulePath;
    string ModuleName;
    string FunctionName;
    string Eigenvalues;
    int MaximumNumberOfClusters;
    int MinimumNumberOfClusters;
    int CurrentGap;
    Clustering clustering_;

    int doIt(vtkDataSet* input,
        vtkDataSet* output);

    bool needsToAbort();

    int updateProgress(const float &progress);

};
