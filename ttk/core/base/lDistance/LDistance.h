/// \ingroup base
/// \class ttk::LDistance 
/// \author Maxime Soler <soler.maxime@total.com>
/// \date 26/02/2017
///
/// \brief TTK %lDistance processing package.
///
/// %LDistance is a TTK processing package that takes a scalar field on the input 
/// and produces a scalar field on the output.
///
/// \sa ttkLDistance.cpp for a usage example.

#ifndef _LDISTANCE_H
#define _LDISTANCE_H

// Standard.
#include                  <cstdlib>
#include                  <string>
#include                  <cmath>

// Base code.
#include                  <Wrapper.h>


namespace ttk {
  
  class LDistance : public Debug {

    public:
        
      LDistance();
      
      ~LDistance();

      template <class dataType>
        int execute(const string &distanceType);

      // computation of Ln on transpose of the matrix
      template <class dataType>
        int computeLn(vector<dataType*> &scalarFields, const int v1, const int v2,
          dataType *output, const int n, const int numberOfScalarFields);
    
      template <class dataType>
        int computeLn(dataType *input1, dataType *input2, dataType *output,
          const int n, const int vertexNumber);

      // computation of Linf on transpose the matrix
      template <class dataType>
        int computeLinf(vector<dataType*> &scalarFields, const int v1, const int v2,
          dataType *output, const int numberOfScalarFields);
      
      template <class dataType>
        int computeLinf(dataType *input1, dataType *input2, dataType *output,
           const int vertexNumber);
    
      /// Pass a pointer to an input array representing a scalarfield.
      /// The expected format for the array is the following:
      /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
      /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
      /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
      /// The array is expected to be correctly allocated. 
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa setVertexNumber() and setDimensionNumber().
      inline int setInputDataPointer1(void *data) {
        inputData1_ = data;
        return 0;
      }
	  
  	  inline int setInputDataPointer2(void *data) {
	    	inputData2_ = data;
    		return 0;
      }

      /// Pass a pointer to an output array representing a scalar field.
      /// The expected format for the array is the following:
      /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
      /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
      /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
      /// The array is expected to be correctly allocated. 
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise. 
      /// \sa setVertexNumber() and setDimensionNumber().
      inline int setOutputDataPointer(void *data) {
        outputData_ = data;
        return 0;
      }

      inline int setNumberOfPoints(long numberOfPoints) {
        numberOfPoints_ = numberOfPoints;
        return 0;
      }
    
      inline double getResult() {
        return result;
      }

      template<typename type>
      static type abs_diff(const type var1, const type var2){
         return (var1 > var2) ? var1 - var2 : var2 - var1;
      }
    
    protected:
    
      void                  *inputData1_,
                            *inputData2_,
                            *outputData_;
      double                result;
      long                  numberOfPoints_;
  };
}

// template functions
template <class dataType> int LDistance::execute(
  const string &distanceType)
{
  
  Timer t;
  int status;
  
  // Check variables consistency
  #ifndef TTK_ENABLE_KAMIKAZE
  if (!inputData1_ || !inputData2_ || !outputData_)
    return -1;
  #endif

  dataType *outputData = (dataType *) outputData_;
  dataType *inputData1 = (dataType *) inputData1_;
  dataType *inputData2 = (dataType *) inputData2_;
  
  int vertexNumber = (int) numberOfPoints_;

  if (distanceType == "inf") {
    status = computeLinf(inputData1, inputData2, outputData, vertexNumber);
  }
  else {
    int n = stoi(distanceType);
    if (n < 1) return -4;
    
    status = computeLn(inputData1, inputData2, outputData, n, vertexNumber);
  }
   
  {
    stringstream msg;
    msg << "[LDistance] Data-set (" << vertexNumber
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return status;
}

template <class dataType> int LDistance::computeLn(
  vector<dataType*> &scalarFields,
  const int v1,
  const int v2,
  dataType *output,
  const int n,
  const int numberOfScalarFields) 
{
  dataType sum = 0;
  
  // Compute difference for each point.
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_) reduction(+:sum)
  #endif
  for (int i = 0; i < numberOfScalarFields; ++i) {
    const dataType diff =  abs_diff<dataType>(scalarFields[i][v1], scalarFields[i][v2]);
    const dataType power = pow(diff, (double)n);
    
    // Careful: huge dataset + huge values
    // may exceed double capacity.
    sum += power;
    
    // Store difference.
    if(output)
      output[i] = power;
  }
  
  sum = pow(sum, 1.0 / (double) n);
  
  // Affect result.
  result = (double) sum;
  {
    stringstream msg;
    msg << "[LDistance] Distance: " << result << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}

template <class dataType> int LDistance::computeLn(
  dataType *input1, dataType *input2, 
  dataType *output, 
  const int n,
  const int vertexNumber) 
{
  dataType sum = 0;
  
  // Compute difference for each point.
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_) reduction(+:sum)
  #endif
  for (int i = 0; i < vertexNumber; ++i) {
    const dataType diff =  abs_diff<dataType>(input1[i], input2[i]);
    const dataType power = pow(diff, (double)n);
    
    // Careful: huge dataset + huge values
    // may exceed double capacity.
    sum += power;
    
    // Store difference.
    if(output)
      output[i] = power;
  }
  
  sum = pow(sum, 1.0 / (double) n);
  
  // Affect result.
  result = (double) sum;
  {
    stringstream msg;
    msg << "[LDistance] Distance: " << result << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}

template <class dataType> int LDistance::computeLinf(
  vector<dataType*> &scalarFields,
  const int v1,
  const int v2,
  dataType *output,
  const int numberOfScalarFields) 
{
  if (numberOfScalarFields< 1) return 0;
  
  dataType maxValue = abs_diff<dataType>(scalarFields[0][v1], scalarFields[0][v2]);
  
  // Compute difference for each point.
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_) reduction(max:maxValue)
  #endif
  for (int i = 1; i < numberOfScalarFields; ++i) {
    const dataType iter = abs_diff<dataType>(scalarFields[i][v1], scalarFields[i][v2]);
    if (iter > maxValue) maxValue = iter;
    
    // Store absolute difference in output.
    if(output)
      output[i] = iter;
  }
  
  // Affect result.
  result = (double) maxValue;
  {
    stringstream msg;
    msg << "[LDistance] Distance: " << result << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}
  
template <class dataType> int LDistance::computeLinf(
  dataType *input1, dataType *input2, 
  dataType *output,
  const int vertexNumber) 
{
  if (vertexNumber < 1) return 0;
  
  dataType maxValue = abs_diff<dataType>(input1[0], input2[0]);
  
  // Compute difference for each point.
  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(threadNumber_) reduction(max:maxValue)
  #endif
  for (int i = 1; i < vertexNumber; ++i) {
    const dataType iter = abs_diff<dataType>(input1[i], input2[i]);
    if (iter > maxValue) maxValue = iter;
    
    // Store absolute difference in output.
    if(output)
      output[i] = iter;
  }
  
  // Affect result.
  result = (double) maxValue;
  {
    stringstream msg;
    msg << "[LDistance] Distance: " << result << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}

#endif // LDISTANCE_H
