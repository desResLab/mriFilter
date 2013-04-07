#ifndef MRIEXCEPTION_H
#define MRIEXCEPTION_H

#include <string>

// BASE EXCEPTION
class MRIException: public std::exception{
  public:
    // Constructor and Destructor
    MRIException(const char* m):msg(m){};
    virtual ~MRIException() throw(){};
	  // Member Functions
		virtual const char* what() const throw() {return msg.c_str();}
  protected:
    std::string msg;
};

// GENERIC MRI SEQUENCE EXCEPTION
class MRISequenceException: public MRIException{
  public:
    MRISequenceException(const char* m):MRIException(m){};
    virtual ~MRISequenceException() throw(){};
    // Member Functions
	  virtual const char* what() const throw() {return msg.c_str();}
};

// PRESSSURE COMPUTATION EXCEPTION
class MRIPressureComputationException: public MRIException{
  public:
    MRIPressureComputationException(const char* m):MRIException(m){};
    virtual ~MRIPressureComputationException() throw(){};
    // Member Functions
	  virtual const char* what() const throw() {return msg.c_str();}
};

// PRESSSURE COMPUTATION EXCEPTION
class MRIMeshCompatibilityException: public MRIException{
  public:
    MRIMeshCompatibilityException(const char* m):MRIException(m){};
    virtual ~MRIMeshCompatibilityException() throw(){};
    // Member Functions
	  virtual const char* what() const throw() {return msg.c_str();}
};

// STATISTIC EVALUATION EXCEPTION
class MRIStatisticsException: public MRIException{
  public:
    MRIStatisticsException(const char* m):MRIException(m){};
    virtual ~MRIStatisticsException() throw(){};
    // Member Functions
	  virtual const char* what() const throw() {return msg.c_str();}
};

#endif // MRIEXCEPTION_H
