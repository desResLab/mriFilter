#ifndef MRIEXCEPTION_H
#define MRIEXCEPTION_H

#include <string>

// BASE EXCEPTION
class mriException: public std::exception{
  public:
    // Constructor and Destructor
    mriException(const char* m):msg(m){};
    virtual ~mriException() throw(){};
	  // Member Functions
		virtual const char* what() const throw() {return msg.c_str();}
  protected:
    std::string msg;
};

#endif // MRIEXCEPTION_H
