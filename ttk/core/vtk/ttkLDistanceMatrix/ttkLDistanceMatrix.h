/// \ingroup vtk
/// \class ttkLDistanceMatrix
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the lDistanceMatrix processing package.
///
/// VTK wrapping code for the @LDistanceMatrix package.
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
/// \sa ttk::LDistanceMatrix
#pragma once

#include<list>

// ttk code includes
#include<LDistanceMatrix.h>
#include<ttkWrapper.h>

// VTK includes
#include<vtkCharArray.h>
#include<vtkDataArray.h>
#include<vtkDataSet.h>
#include<vtkDataSetAlgorithm.h>
#include<vtkDoubleArray.h>
#include<vtkFiltersCoreModule.h>
#include<vtkFloatArray.h>
#include<vtkStringArray.h>
#include<vtkInformation.h>
#include<vtkIntArray.h>
#include<vtkObjectFactory.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>
#include<vtkPointData.h>
#include<vtkStreamingDemandDrivenPipeline.h>

class VTKFILTERSCORE_EXPORT ttkLDistanceMatrix
: public vtkDataSetAlgorithm, public Wrapper{

  public:

    static ttkLDistanceMatrix* New();
    vtkTypeMacro(ttkLDistanceMatrix, vtkDataSetAlgorithm)

      // default ttk setters
      vtkSetMacro(debugLevel_, int);

    void SetThreadNumber(int threadNumber){\
      ThreadNumber = threadNumber;\
        SetThreads();\
    }\
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    vtkSetMacro(DistanceType, string);
    vtkGetMacro(DistanceType, string);

    vtkSetMacro(TransposeMatrix, bool);
    vtkGetMacro(TransposeMatrix, bool);

    int FillOutputPortInformation(int port, vtkInformation *info){
      if(port==0)
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      return 1;
    }

  protected:

    ttkLDistanceMatrix(){
      UseAllCores=false;
      ThreadNumber=1;
      TransposeMatrix=false;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~ttkLDistanceMatrix(){};

    TTK_SETUP();

  private:

    string DistanceType;
    bool TransposeMatrix;
    LDistanceMatrix lDistanceMatrix_;
};
