/// \ingroup vtk
/// \class ttkSortScalarField
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the sortScalarField processing package.
///
/// VTK wrapping code for the @SortScalarField package.
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
/// \sa ttk::SortScalarField
#pragma once

// VTK includes -- to adapt
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkCellData.h>
#include                  <vtkSmartPointer.h>

// ttk code includes
#include                  <Wrapper.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkSortScalarField
#else
class ttkSortScalarField
#endif
: public vtkDataSetAlgorithm, public Wrapper{

  public:

    static ttkSortScalarField* New();
    vtkTypeMacro(ttkSortScalarField, vtkDataSetAlgorithm)

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

  protected:

    ttkSortScalarField(){
      UseAllCores = false;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~ttkSortScalarField(){};

    int RequestData(vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector);

  private:

    bool UseAllCores;
    int ThreadNumber;

    int doIt(vtkDataSet *input, vtkDataSet *output);
    bool needsToAbort();
    int updateProgress(const float &progress);
};
