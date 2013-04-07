#ifndef SCHMESSAGES_H
#define SCHMESSAGES_H

#include <stdio.h>
#include <string>


// Handle the possibiliy of writing messages to some Error Log or to Native Client Instances
inline void WriteSchMessage(std::string Msg){
    printf("%s",Msg.c_str());
}

inline void WriteHeader(){
	WriteSchMessage(std::string("\n"));
	WriteSchMessage(std::string("------ Flow Manipulation Toolkit - Dr. Daniele Schiavazzi, 2012 -------\n"));
	WriteSchMessage(std::string("-----------------------------------------------------------------------\n"));
	WriteSchMessage(std::string("---- Filtering and pressure estimation of acquired velocity fields ----\n"));
	WriteSchMessage(std::string("-----------------------------------------------------------------------\n"));
	WriteSchMessage(std::string("Current Release: 0.1 beta\n"));
	WriteSchMessage(std::string("\n"));
	WriteSchMessage(std::string("Usage: fmt inputFile outputFile\n"));
	WriteSchMessage(std::string("\n"));
}

#endif //SCHMESSAGES_H
