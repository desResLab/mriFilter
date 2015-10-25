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
    WriteSchMessage(std::string("-----------------------------------\n"));
    WriteSchMessage(std::string(" Flow Manipulation Toolkit\n"));
    WriteSchMessage(std::string(" 2013 - Daniele Schiavazzi,Ph.D.\n"));
    WriteSchMessage(std::string(" Release: 0.2\n"));
    WriteSchMessage(std::string("-----------------------------------\n"));
	WriteSchMessage(std::string("\n"));
}

#endif //SCHMESSAGES_H
