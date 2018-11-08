/// \ingroup base
/// \class ttk::PersistenceMap 
/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>
/// \date February 2018.
///
/// \brief TTK %persistenceMap processing package.
///
/// %PersistenceMap is a TTK processing package that computes the persistent 
/// map of a scalar field.
///
/// \sa ttk::Triangulation
/// \sa vtkPersistenceMap.cpp %for a usage example.

#pragma once

// base code includes
#include                  <FTMTreePP.h>
#include                  <Geometry.h>
#include                  <Triangulation.h>
#include                  <Wrapper.h>


namespace ttk{
  
  class PersistenceMap : public Debug{

    public:
        
      PersistenceMap();
      
      ~PersistenceMap();

      /// Execute the package.
      /// \pre If this TTK package uses ttk::Triangulation for fast mesh 
      /// traversals, the function setupTriangulation() must be called on this 
      /// object prior to this function, in a clearly distinct pre-processing 
      /// steps. An error will be returned otherwise.
      /// \note In such a case, it is recommended to exclude 
      /// setupTriangulation() from any time performance measurement.
      /// \param persistenceScale Persistence scale.
      /// \return Returns 0 upon success, negative values otherwise.
      template <class dataType>
        int execute(
          const double &persistenceThreshold,
          const double &persistenceScale, 
          const double &persistenceSpread,
          const double &defaultAmplitude,
          const double &defaultDeviation,
          const int &criticalType);
    
      /// Pass a pointer to an input array representing a scalarfield.
      /// The expected format for the array is the following:
      /// <vertex0-component0> <vertex0-component1> ... <vertex0-componentN>
      /// <vertex1-component0> <vertex1-component1> ... <vertex1-componentN>
      /// <vertexM-component0> <vertexM-component1> ... <vertexM-componentN>.
      /// The array is expected to be correctly allocated. 
      /// \param data Pointer to the data array.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa setVertexNumber() and setDimensionNumber().
      inline int setInputDataPointer(void *data){
        inputData_ = data;
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
      inline int setOutputDataPointer(void *data){
        outputData_ = (double *) data;
        return 0;
      }
     
      /// Setup a (valid) triangulation object for this TTK base object.
      ///
      /// \pre This function should be called prior to any usage of this TTK 
      /// object, in a clearly distinct pre-processing step that involves no 
      /// traversal or computation at all. An error will be returned otherwise.
      ///
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement. Therefore, it is recommended to 
      /// call this function ONLY in the pre-processing steps of your program. 
      /// Note however, that your triangulation object must be valid when 
      /// calling this function (i.e. you should have filled it at this point, 
      /// see the setInput*() functions of ttk::Triangulation). See vtkPersistenceMap 
      /// for further examples.
      ///
      /// \param triangulation Pointer to a valid triangulation.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa ttk::Triangulation
      inline int setupTriangulation(Triangulation *triangulation){
        triangulation_ = triangulation;
       
        if(triangulation_){
          persistencePairs_.setupTriangulation(triangulation_);
        }
        
        return 0;
      }
    
    protected:
    
      ftm::FTMTreePP        persistencePairs_;
      void                  *inputData_;
      double                *outputData_;
      Triangulation         *triangulation_;
      vector<int>           offsets_;
      
  };
}

// template functions
template <class dataType> int PersistenceMap::execute(
  const double &persistenceThreshold,
  const double &persistenceScale, const double &persistenceSpread,
  const double &defaultAmplitude, const double &defaultDeviation,
  const int &criticalType){
    
  {
    stringstream msg;
    msg << "[PersistenceMap] Computing persistence map..." << endl;
    msg << "[PersistenceMap]   persistenceThreshold: " << persistenceThreshold
      << endl;
    msg << "[PersistenceMap]   persistenceScale: " 
      << persistenceScale << endl;
    msg << "[PersistenceMap]   persistenceSpread: " 
      << persistenceSpread << endl;
    msg << "[PersistenceMap]   defaultAmplitude: " 
      << defaultAmplitude << endl;
    msg << "[PersistenceMap]   defaultDeviation: " 
      << defaultDeviation << endl;
    msg << "[PersistenceMap]   criticalType: " 
      << criticalType << endl;
    dMsg(cout, msg.str(), Debug::infoMsg);
  }
  
    Timer t;
    
    // check the consistency of the variables -- to adapt
    #ifndef TTK_ENABLE_KAMIKAZE
    if(!triangulation_)
      return -1;
    if(!inputData_)
      return -2;
    if(!outputData_)
      return -3;
    if((persistenceThreshold < 0)||(persistenceThreshold > 1))
      return -4;
    #endif
    
    int vertexNumber = triangulation_->getNumberOfVertices();
    offsets_.resize(vertexNumber);
    for(int i = 0; i < vertexNumber; i++){
      offsets_[i] = i;
    }
    
    // compute the persistence
    persistencePairs_.setVertexScalars(inputData_);
    persistencePairs_.setTreeType(ftm::TreeType::Join_Split);
    persistencePairs_.setVertexSoSoffsets(offsets_.data());
    persistencePairs_.setSegmentation(false);
    persistencePairs_.setThreadNumber(1);
    persistencePairs_.build<dataType>();
    
    vector<tuple<ftm::idVertex, ftm::idVertex, dataType> > JTPairs;
    vector<tuple<ftm::idVertex, ftm::idVertex, dataType> > STPairs;

    if((criticalType == 0)||(criticalType == 2))
      persistencePairs_.computePersistencePairs<dataType>(JTPairs, true);
    if(criticalType >= 1)
      persistencePairs_.computePersistencePairs<dataType>(STPairs, false);
   
    double maxPersistenceMin = -2, secondMaxPersistenceMin = -2;
    double maxPersistenceMax = -2, secondMaxPersistenceMax = -2;
    if(criticalType >= 1){
      for(int i = 0; i < (int) STPairs.size(); i++){
        if(maxPersistenceMax < get<2>(STPairs[i])){
          secondMaxPersistenceMax = maxPersistenceMax;
          maxPersistenceMax = get<2>(STPairs[i]);
        }
        else if((get<2>(STPairs[i]) != maxPersistenceMax)
          &&(secondMaxPersistenceMax < get<2>(STPairs[i]))){
          secondMaxPersistenceMax = get<2>(STPairs[i]);
        }
      }
    }
    if((criticalType == 0)||(criticalType == 2)){
      for(int i = 0; i < (int) JTPairs.size(); i++){
        if(maxPersistenceMin < get<2>(JTPairs[i])){
          secondMaxPersistenceMin = maxPersistenceMin;
          maxPersistenceMin = get<2>(JTPairs[i]);
        }
        else if((get<2>(JTPairs[i]) != maxPersistenceMin)
          &&(secondMaxPersistenceMin < get<2>(JTPairs[i]))){
          secondMaxPersistenceMin = get<2>(JTPairs[i]);
        }
      }
    }
    
    if(secondMaxPersistenceMax > maxPersistenceMax/5.0)
      maxPersistenceMax = secondMaxPersistenceMax;
    if(secondMaxPersistenceMin > maxPersistenceMin/5.0)
      maxPersistenceMin = secondMaxPersistenceMin;
    
    if(criticalType >= 1){
      for(int i = 0; i < (int) STPairs.size(); i++){
        if(get<2>(STPairs[i]) >  maxPersistenceMax)
          get<2>(STPairs[i]) = maxPersistenceMax;
      }
    }
    if((criticalType == 0)||(criticalType == 2)){
      for(int i = 0; i < (int) JTPairs.size(); i++){
        if(get<2>(JTPairs[i]) >  maxPersistenceMin)
          get<2>(JTPairs[i]) = maxPersistenceMin;
      }
    }
    
    double maxPersistence = 0;
    if(criticalType == 0) 
      maxPersistence = maxPersistenceMin;
    else if(criticalType == 1) 
      maxPersistence = maxPersistenceMax;
    else if(maxPersistenceMin < maxPersistenceMax)
      maxPersistence = maxPersistenceMin;
    else
      maxPersistence = maxPersistenceMax;

    {
      stringstream msg;
      msg << "[PersistenceMap] Maximum persistences: " << maxPersistenceMax
        << " (second: " << secondMaxPersistenceMax << ")"
        << endl;
      msg << "[PersistenceMap] Minimum persistences: " << maxPersistenceMin
        << " (second: " << secondMaxPersistenceMin << ")"
        << endl;
      msg << "[PersistenceMap] Selected maximum persistence: " 
        << maxPersistence << endl;
      dMsg(cout, msg.str(), infoMsg);
    }
    
    bool hitMin = false, hitMax = false;
    vector<double> bboxMin(3, 0), bboxMax(3, 0);
    for(int i = 0 ; i < vertexNumber; i++){
      vector<float> p(3, 0);
      triangulation_->getVertexPoint(i, p[0], p[1], p[2]);
      for(int j = 0; j < triangulation_->getDimensionality(); j++){
        if(!std::isnan(p[j])){
          if((!hitMin)||(p[j] < bboxMin[j])){
            bboxMin[j] = p[j];
            hitMin = true;
          }
          if((!hitMax)||(p[j] > bboxMax[j])){
            bboxMax[j] = p[j];
            hitMax = true;
          }
        }
      }
    }

    double maxDistance = Geometry::distance(bboxMin.data(), bboxMax.data());
    
    {
      stringstream msg;
      msg << "[PersistenceMap] Max distance: " << maxDistance << endl;
      dMsg(cout, msg.str(), infoMsg);
    }
    
    // the following open-mp processing is only relevant for embarrassingly 
    // parallel algorithms (such as smoothing) -- to adapt
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(threadNumber_) 
    #endif
    for(int i = 0; i < vertexNumber; i++){
      vector<float> p(3, 0);
      triangulation_->getVertexPoint(i, p[0], p[1], p[2]);
      
      outputData_[i] = 0;
     
      if(criticalType >= 1){
        for(int j = 0; j < (int) STPairs.size(); j++){
          
          double persistence = get<2>(STPairs[j])/maxPersistence;
          
          if(persistence > persistenceThreshold){
            
            if(persistence > 1)
              persistence = 1;
            
            vector<float> cP(3, 0);
            int cpId = get<0>(STPairs[j]);
            triangulation_->getVertexPoint(cpId, cP[0], cP[1], cP[2]);
            double distance = 
              Geometry::distance(p.data(), cP.data()) / maxDistance;
            double sigma = 
              defaultDeviation 
                + (persistenceSpread*persistence);
            
            outputData_[i] += 
              (defaultAmplitude 
                + persistenceScale*persistence)
                *exp(-(distance*distance)/(2*sigma*sigma));
          }
        }
      }
      
      if((criticalType == 0)||(criticalType == 2)){
        for(int j = 0; j < (int) JTPairs.size(); j++){
          
          double persistence = get<2>(JTPairs[j])/maxPersistence;
          if(persistence > persistenceThreshold){
          
            if(persistence > 1)
              persistence = 1;
            
            vector<float> cP(3, 0);
            int cpId = get<0>(JTPairs[j]);
            triangulation_->getVertexPoint(cpId, cP[0], cP[1], cP[2]);
            double distance = 
              Geometry::distance(p.data(), cP.data()) / maxDistance;
          
            double sigma = 
              defaultDeviation 
                + (persistenceSpread*persistence);
           
            outputData_[i] += 
              (defaultAmplitude 
                + persistenceScale*persistence)
                *exp(-(distance*distance)/(2*sigma*sigma));
          }
        }
      }
    }
    
    {
      stringstream msg;
      msg << "[PersistenceMap] Data-set (" << vertexNumber
      << " points) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
    
    return 0;
}

// if the package is a pure template class, uncomment the following line
// #include                  <PersistenceMap.cpp>
