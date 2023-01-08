#ifndef BEMUSE_IO_H
#define BEMUSE_IO_H

#include "BEMUse_Math_Types.h"
#include <iostream>             // Debugging
#include <fstream>              // Open / write to files
#include <iomanip>              // Set output string precision

//---- Input string handling

inline std::vector<std::string> Split(const std::string &text, char sep)
{
    std::vector<std::string> tokens;
    std::size_t start = 0, end = 0;
    while ((end = text.find(sep, start)) != std::string::npos)
    {
        if (end != start) {
          tokens.push_back(text.substr(start, end - start));
        }
        start = end + 1;
    }
    if (end != start) {
       tokens.push_back(text.substr(start));
    }
    return tokens;
}

inline std::vector<std::string> Split2(const std::string &text, char sep1, char sep2)
{
    std::vector<std::string> tokens;
    std::size_t start = 0, end = 0;
    while ( ((end = text.find(sep1, start)) != std::string::npos) || ((end = text.find(sep2, start)) != std::string::npos))
    {
        if (end != start) {
          tokens.push_back(text.substr(start, end - start));
        }
        start = end + 1;
    }
    if (end != start) {
       tokens.push_back(text.substr(start));
    }
    return tokens;
}

inline std::vector<std::string> SplitSpacesTabs(const std::string &text) {return Split2(text,' ','\t');}

//--- Checking string matching

inline bool Contains(std::string S1,std::string S2)   {return (S1.find(S2) != std::string::npos);}

//--- Writing floating point values with specified width

inline void PadString(std::string &str, const size_t num, const char paddingChar = ' ')
{
    if (num > str.size())   str.insert(str.end(), num - str.size(), paddingChar);
}

inline std::string PFS(const Real &A, int NPad=15)
{
    std::ostringstream S;
    S << A;
    S << std::setprecision(6);
    std::string String_Out = S.str();
    if (A>=0)   String_Out.insert(0," ");       // Pad the string so that negatives and positives align
    PadString(String_Out,NPad);
    return String_Out;
}

inline std::string PST(std::string str, int NPad=15)
{
    std::string String_Copy = str;
    PadString(String_Copy,NPad);
    return String_Copy;
}

//--- Create directory

#include <string>       // required for std::string
#include <sys/types.h>  // required for stat.h
#include <sys/stat.h>   // no clue why required -- man pages say so

inline int CreateDirectory(std::string Path)
{
    #if defined(_WIN32)
        return mkdir(Path.c_str());
    #else
        mode_t nMode = 0733; // UNIX style permissions
        return mkdir(Path.c_str(),nMode);
    #endif
}

#endif // BEMUSE_IO_H
