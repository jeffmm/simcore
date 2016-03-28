//Class to help with testing the time of the code
#ifndef _TIMETESTER_H_ 
#define  _TIMETESTER_H_

#include <iostream>
#include <string>

class Timer{
    private:
        double init_,
               start_,
                stop_;
    public:
        Timer();
        void StartTimer();
        void StopTimer();
        void PrintDiffTime();
        void PrintDiffTime(std::string);
        
};

#endif //_TIMETESTER_H_ 
