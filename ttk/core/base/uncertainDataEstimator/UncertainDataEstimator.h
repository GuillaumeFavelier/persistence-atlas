/// \ingroup base
/// \class ttk::UncertainDataEstimator
/// \author Michael Michaux <michauxmichael89@gmail.com>
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date August 2016.
///
/// \brief TTK processing package that takes an input ensemble data set 
/// (represented by a list of scalar fields) and which computes various 
/// vertexwise statistics (PDF estimation, bounds, moments, etc.)
///
/// \sa ttkUncertainDataEstimator.cpp %for a usage example.

#ifndef _UNCERTAINDATAESTIMATOR_H
#define _UNCERTAINDATAESTIMATOR_H

// base code includes
#include                  <Wrapper.h>

namespace ttk{

  template <class dataType>
    class PDFBounds : public Debug {
      public:

        PDFBounds() {
          numberOfVertices_ = 0;
        }

        ~PDFBounds() {
          flush();
        }

        int evaluateRealization(const void *voidPointer) {
#ifdef TTK_ENABLE_KAMIKAZE
          if(!(numberOfVertices_>0)) {
            return -1; // Number of vertices not defined
          }
#endif
          const dataType *inputData =
            reinterpret_cast<const dataType*>(voidPointer);
          unsigned int numberOfVertices =
            static_cast<unsigned int>(numberOfVertices_);
          /* Initialize if first call since a change */
          if (!(upperBound_.size()==numberOfVertices) || !(lowerBound_.size()==numberOfVertices)) {
            upperBound_.resize(numberOfVertices);
            lowerBound_.resize(numberOfVertices);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
            for (size_t i = 0 ; i < numberOfVertices ; i++) {
              upperBound_[i] = inputData[i];
              lowerBound_[i] = inputData[i];
            }
          } else { /* Update the two fields with the new input */
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
            for (size_t i = 0 ; i < numberOfVertices ; i++) {
              // Upper Bound
              if (inputData[i] > upperBound_[i]) {
                upperBound_[i] = inputData[i];
              }
              // Lower Bound
              if (inputData[i] < lowerBound_[i]) {
                lowerBound_[i] = inputData[i];
              }
            }
          }
          return 0;
        }

        pair<dataType,dataType> getRange() const {
          pair<dataType,dataType> range;
          range.first = getRangeMin();
          range.second = getRangeMax();
          return range;
        }

        dataType getRangeMax() const {
          if(upperBound_.size()) {
            dataType maxValue = upperBound_[0];
            for (size_t i = 1; i < upperBound_.size(); i++) {
              if (upperBound_[i] > maxValue) {
                maxValue = upperBound_[i];
              }
            }
            return maxValue;
          } else {
            return 0;
          }
        }

        dataType getRangeMin() const {
          if(lowerBound_.size()) {
            dataType minValue = lowerBound_[0];
            for (size_t i = 1; i < lowerBound_.size(); i++) {
              if (lowerBound_[i] < minValue) {
                minValue = lowerBound_[i];
              }
            }
            return minValue;
          } else {
            return 0;
          }
        }

        inline int flush() {
          numberOfVertices_ = 0;
          upperBound_.clear();
          lowerBound_.clear();
          return 0;
        }

        inline dataType* getLowerBoundPointer() {
          return lowerBound_.data();
        }

        inline dataType* getUpperBoundPointer() {
          return upperBound_.data();
        }

        inline int setNumberOfVertices(const int number) {
          numberOfVertices_ = number;
          return 0;
        }

      protected:
        int numberOfVertices_;
        vector<dataType> upperBound_;
        vector<dataType> lowerBound_;
    };

  class PDFHistograms : public Debug {
    public:

      PDFHistograms() {
        numberOfBins_ = 0;
        numberOfInputFields_ = 0;
        numberOfVertices_ = 0;
        rangeMax_ = 0;
        rangeMin_ = 0;
      }

      ~PDFHistograms() {
        flush();
      }

      template <class dataType>
        int evaluateRealization(const dataType *inputData) {
#ifdef TTK_ENABLE_KAMIKAZE
          if(!(rangeMin_<rangeMax_)) {
            return -1; // Range error
          }
          if(!(numberOfBins_>0)) {
            return -2; // Number of bins not defined
          }
          if(!(numberOfVertices_>0)) {
            return -3; // Number of vertices not defined
          }
#endif
          if(numberOfInputFields_ == 0) {
            /* Initialize */
            probability_.resize(numberOfBins_);
            double dx = (rangeMax_-rangeMin_) / static_cast<double>(numberOfBins_);
            for (size_t i=0 ; i < numberOfBins_ ; i++) {
              probability_[i].resize(numberOfVertices_);
              binValue_[i] = rangeMin_+(dx/2.0)+(static_cast<double>(i)*dx);
            }
          }
          /* Add input datas */
          for(unsigned int i=0 ; i<numberOfVertices_ ; i++) {
            unsigned int bin = static_cast<unsigned int>(floor((inputData[i]-rangeMin_)*numberOfBins_/(rangeMax_-rangeMin_)));
            bin = (bin == numberOfBins_) ? numberOfBins_-1 : bin;
            probability_[bin][i] += 1.0;
          }
          numberOfInputFields_++;
          return 0;
        }

      inline int flush() {
        binValue_.clear();
        numberOfBins_ = 0;
        numberOfInputFields_ = 0;
        numberOfVertices_ = 0;
        rangeMax_ = 0;
        rangeMin_ = 0;
        // selection_.clear(); // TODO : selection support
        return 0;
      }

      inline double* getBinFieldPointer(const unsigned int binId) {
        if(binId < numberOfBins_) {
          return probability_[binId].data();
        } else {
          return nullptr;
        }
      }

      void getVertexHistogram(const unsigned int vertexId, vector<double> &histogram) const {
        histogram.resize(numberOfBins_);
        if(vertexId < numberOfVertices_) {
#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp parallel for num_threads(threadNumber_)
#else
#pragma omp parallel for num_threads(threadNumber_) \
          schedule(static, numberOfBins_/threadNumber_)
#endif
#endif
          for(int i=0 ; i< (int) numberOfBins_ ; i++) {
            if(probability_[i].size()==numberOfVertices_) {
              histogram[i] = probability_[i][vertexId];
            } else {
              histogram[i] = 0.0;
            }
          }
        } else {
          fill(histogram.begin(), histogram.end(), 0.0);
        }
      }

      int normalize() {
        const double normalization = 1.0 / static_cast<double>(numberOfInputFields_);
#ifdef TTK_ENABLE_OPENMP
#ifdef _WIN32
#pragma omp parallel for num_threads(threadNumber_)
#else
#pragma omp parallel for num_threads(threadNumber_) collapse(2) \
        schedule(static, (numberOfBins_*numberOfVertices_)/threadNumber_)
#endif
#endif
        for(int i=0 ; i< (int) numberOfBins_ ; i++) {
          for(int j=0 ; j< (int) numberOfVertices_ ; j++) {
            probability_[i][j] *= normalization;
          }
        }
        return 0;
      }

      inline int setNumberOfBins(const unsigned int number) {
        numberOfBins_ = number;
        return 0;
      }

      inline int setNumberOfVertices(const unsigned int number) {
        numberOfVertices_ = number;
        return 0;
      }

      inline int setRange(const double min, const double max) {
        rangeMin_ = min;
        rangeMax_ = max;
        return 0;
      }

    protected:
      vector<double>            binValue_;
      vector<vector<double> >   probability_;
      unsigned int              numberOfBins_;
      unsigned int              numberOfInputFields_;
      unsigned int              numberOfVertices_;
      double                    rangeMin_;
      double                    rangeMax_;
      // vector<int> selection_; // TODO : selection support
  };

  class UncertainDataEstimator : public Debug{

    public:

      UncertainDataEstimator();

      ~UncertainDataEstimator();

      /// Execute the package.
      /// \return Returns 0 upon success, negative values otherwise.
      template <class dataType>
        int execute() const;

      inline int setInputScalarFields(vector<void*> *data){
        inputData_ = data;
        return 0;
      }

      inline int setOutputLowerBoundFields(vector<void*> *data){
        outputLowerBoundFields_ = data;
        return 0;
      }

      inline int setOutputUpperBoundFields(vector<void*> *data){
        outputUpperBoundFields_ = data;
        return 0;
      }

      inline int setOutputMeanFields(vector<void*> *data){
        outputMeanFields_ = data;
        return 0;
      }

      inline int setOutputProbabilityFields(vector<vector<void*>> *data){
        outputProbabilityFields_ = data;
        return 0;
      }

      inline int setComputeLowerBound(const bool &state){
        computeLowerBound_ = state;
        return 0;
      }


      inline int setComputeUpperBound(const bool &state){
        computeUpperBound_ = state;
        return 0;
      }

      inline int setBinCount(const int binCount){
        binCount_ = binCount;
        return 0;
      }

      /// Set the number of vertices in the scalar field.
      /// \param vertexNumber Number of vertices in the data-set.
      /// \return Returns 0 upon success, negative values otherwise.
      inline int setVertexNumber(const int &vertexNumber){
        vertexNumber_ = vertexNumber;
        return 0;
      }

      inline int setNumberOfInputFields(int numberOfInputFields){
        numberOfInputFields_ = numberOfInputFields;
        return 0;
      }

      inline int setNumberOfClusters(const int numberOfClusters){
        numberOfClusters_=numberOfClusters;
        return 0;
      }

      inline int setClusteringData(void* data){
        clusteringData_=data;
        return 0;
      }

    protected:

      int vertexNumber_;
      int numberOfInputFields_;
      bool computeLowerBound_;
      bool computeUpperBound_;
      int binCount_;
      vector<void*>* inputData_;
      int numberOfClusters_;
      void* clusteringData_;
      vector<void*>* outputLowerBoundFields_;
      vector<void*>* outputUpperBoundFields_;
      vector<void*>* outputMeanFields_;
      vector<vector<void*>>* outputProbabilityFields_;
  };
}

// template functions
template <class dataType> int UncertainDataEstimator::execute() const{
  Timer t;

  // Check the consistency of the variables
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfInputFields_)
    return -1;
  if(!vertexNumber_)
    return -2;
  if(!inputData_)
    return -3;
  if(!outputLowerBoundFields_)
    return -4;
  if(!outputUpperBoundFields_)
    return -5;
#endif

  // Pointers type casting
  vector<dataType*> inputData(numberOfInputFields_);
  for(int i=0; i<numberOfInputFields_; ++i)
    inputData[i]=static_cast<dataType*>((*inputData_)[i]);
  const int numberOfOutputFields=(numberOfClusters_)?numberOfClusters_:1;
  vector<double*> outputLowerBoundFields(numberOfOutputFields, nullptr);
  vector<double*> outputUpperBoundFields(numberOfOutputFields, nullptr);
  vector<double*> outputMeanFields(numberOfOutputFields, nullptr);
  vector<vector<double*>> outputProbabilityFields(numberOfOutputFields);
  for(int i=0; i<numberOfOutputFields; ++i){
    if(computeLowerBound_)
      outputLowerBoundFields[i]=static_cast<double*>((*outputLowerBoundFields_)[i]);
    if(computeUpperBound_)
      outputUpperBoundFields[i]=static_cast<double*>((*outputUpperBoundFields_)[i]);
    outputMeanFields[i]=static_cast<double*>((*outputMeanFields_)[i]);

    if(computeUpperBound_ && computeLowerBound_){
      outputProbabilityFields[i].resize(binCount_, nullptr);
      for(int j=0; j<binCount_; ++j){
        outputProbabilityFields[i][j]=static_cast<double*>((*outputProbabilityFields_)[i][j]);
        std::fill(outputProbabilityFields[i][j], outputProbabilityFields[i][j]+vertexNumber_, 0);
      }
    }
  }
  int* clusteringData=static_cast<int*>(clusteringData_);

  if(numberOfClusters_){
    for(int c = 0; c < numberOfClusters_; c++){
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(int v = 0; v < (int) vertexNumber_; v++){
        // For the lower bound scalar field
        if(computeLowerBound_){
          bool isFirst=true;
          // Loop over the inputs
          for(int inp=0 ; inp < numberOfInputFields_ ; inp++){
            if(clusteringData[inp]==c){
              if(isFirst){
                outputLowerBoundFields[c][v] = inputData[inp][v];
                isFirst=false;
              }
              else{
                if(inputData[inp][v] < outputLowerBoundFields[c][v])
                  outputLowerBoundFields[c][v] = inputData[inp][v];
              }
            }
          }
        }

        // For the upper bound scalar field
        if(computeUpperBound_){
          bool isFirst=true;
          // Loop over the inputs
          for(int inp=0 ; inp < numberOfInputFields_ ; inp++){
            if(clusteringData[inp]==c){
              if(isFirst){
                outputUpperBoundFields[c][v] = inputData[inp][v];
                isFirst=false;
              }
              else{
                if(inputData[inp][v] > outputUpperBoundFields[c][v])
                  outputUpperBoundFields[c][v] = inputData[inp][v];
              }
            }
          }
        }

        // For the mean scalar field
        int n=0;
        outputMeanFields[c][v] = 0;
        for(int inp=0; inp < numberOfInputFields_ ; inp++){
          if(clusteringData[inp]==c){
            outputMeanFields[c][v] += inputData[inp][v];
            n++;
          }
        }
        // apriori: n should be always greater or equal to 1
        outputMeanFields[c][v] /= n;

        // Histogram
        if(computeUpperBound_ && computeLowerBound_){
          const double range[2]{
            outputLowerBoundFields[c][v],
            outputUpperBoundFields[c][v]
          };

          if(range[0]==range[1]){
            for(int inp=0 ; inp<numberOfInputFields_ ; inp++){
              if(clusteringData[inp]==c){
                for(int b=0; b<binCount_; b++)
                  outputProbabilityFields[c][b][v]=-1;
              }
            }

            continue;
          }

          const double mul = double(binCount_) / (range[1]-range[0]);

          for(int inp=0 ; inp<numberOfInputFields_ ; inp++){
            if(clusteringData[inp]==c){
              int idx = floor((inputData[inp][v]-range[0])*mul);
              idx = (idx==binCount_) ? binCount_-1 : idx;
              outputProbabilityFields[c][idx][v]++;
            }
          }
        }
      }
    }
  }
  else{
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int v = 0; v < (int) vertexNumber_; v++){
      // For the lower bound scalar field
      if(computeLowerBound_){
        // Initialisation : values of the first input
        outputLowerBoundFields[0][v] = inputData[0][v];
        // Loop over the inputs
        for(int inp=1 ; inp < numberOfInputFields_ ; inp++){
          // Minimum value
          if(inputData[inp][v] < outputLowerBoundFields[0][v])
            outputLowerBoundFields[0][v] = inputData[inp][v];
        }
      }

      // For the upper bound scalar field
      if(computeUpperBound_){
        // Initialisation : values of the first input
        outputUpperBoundFields[0][v] = inputData[0][v];
        // Loop over the inputs
        for(int inp=1 ; inp < numberOfInputFields_ ; inp++){
          // Maximum value
          if(inputData[inp][v] > outputUpperBoundFields[0][v])
            outputUpperBoundFields[0][v] = inputData[inp][v];
        }
      }

      // For the mean scalar field
      outputMeanFields[0][v] = 0;
      for(int inp=0; inp < numberOfInputFields_ ; inp++)
        outputMeanFields[0][v] += inputData[inp][v];
      outputMeanFields[0][v] /= numberOfInputFields_;

      // Histogram
      if(computeUpperBound_ && computeLowerBound_){
        const double range[2]{
          outputLowerBoundFields[0][v],
          outputUpperBoundFields[0][v]
        };

        if(range[0]==range[1]){
          for(int inp=0 ; inp<numberOfInputFields_ ; inp++){
            for(int b=0; b<binCount_; b++)
              outputProbabilityFields[0][b][v]=-1;
          }

          continue;
        }

        const double mul = double(binCount_) / (range[1]-range[0]);

        for(int i=0 ; i<numberOfInputFields_ ; i++){
          int idx = floor((inputData[i][v]-range[0])*mul);
          idx = (idx==binCount_) ? binCount_-1 : idx;
          outputProbabilityFields[0][idx][v]++;
        }
      }
    }
  }

  {
    stringstream msg;
    msg << "[UncertainDataEstimator] Data-set (" << vertexNumber_
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  return 0;
}

#endif // UNCERTAINDATAESTIMATOR_H
