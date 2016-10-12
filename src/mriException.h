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

#endif // MRIEXCEPTION_H
