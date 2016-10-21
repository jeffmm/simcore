#ifndef _SIMCORE_POSIT_READER_H_
#define _SIMCORE_POSIT_READER_H_
#include "auxiliary.h"

class PositReader {
  private:
    int nchar_,
        n_posit_,
        n_steps_,
        n_objs_,
        max_objs_ = -1,
        n_misc_ = 0,
        count_ = 0; // XXX REMOVE THIS
    double *position_ = nullptr,
           *scaled_position_,
           *orientation_,
           *diameter_,
           *length_,
           *misc_;
    bool print_ = false;
    std::fstream ip_;
    SID sid_posit_;

    // Converts 2d array index to 1d contiguous array index, assuming 3 dimensions
    int ix(int x, int y) {
      return 3*x+y;
    }
    int misc_ix(int x, int y) {
      return n_misc_*x+y;
    }

    void OpenPosit(std::string file_name) {
      ip_.open(file_name, std::ios::binary | std::ios::in);
      if (!ip_.is_open()) {
        std::cout << "Failed to open posit file " << file_name << "!\n";
        exit(1);
      }
    }

    void ClosePosit() {
      ip_.close();
    }

    // Reallocate arrays in case number of objects in system changes
    void ReallocateObjs() {
      max_objs_ = n_objs_;
      DeallocateObjs();
      AllocateObjs();
    }

    void DeallocateObjs() {
      // Avoid deallocation before allocation
      if (position_ != nullptr) {
        delete[] position_;
        delete[] scaled_position_;
        delete[] orientation_;
        delete[] diameter_;
        delete[] length_;
        if (n_misc_ > 0)
          delete[] misc_;
      }
    }
    void AllocateObjs() {
      if (max_objs_ < 0)
        error_exit("ERROR: max_objs_ unallocated in posit reader!\n");
      position_ = new double[max_objs_*3];
      scaled_position_ = new double[max_objs_*3];
      orientation_ = new double[max_objs_*3];
      diameter_ = new double[max_objs_];
      length_ = new double[max_objs_];
      if (n_misc_>0)
        misc_ = new double[max_objs_*n_misc_];
    }

    // Simply prints the information from the posit file for debugging purposes
    void PrintAll() {
      if (!print_) return;
      for (int i=0; i<n_objs_; ++i) {
        printf("position_ %d: {%2.2f, %2.2f, %2.2f}\n",i+1,position_[ix(i,0)],position_[ix(i,1)],position_[ix(i,2)]);
        printf("scaled_position_ %d: {%2.2f, %2.2f, %2.2f}\n",i+1,scaled_position_[ix(i,0)],scaled_position_[ix(i,1)],scaled_position_[ix(i,2)]);
        printf("orientation_ %d: {%2.2f, %2.2f, %2.2f}\n",i+1,orientation_[ix(i,0)],orientation_[ix(i,1)],orientation_[ix(i,2)]);
        printf("diameter_ %d: %2.2f\n",i+1,diameter_[i]);
        printf("length_ %d: %2.2f\n",i+1,length_[i]);
        for (int j=0; j<n_misc_; ++j)
          printf("misc_%d %d: %2.2f\n",j,i,misc_[misc_ix(i,j)]);
      }
    }

    // Read n_posit, n_steps, and species type etc;
    void Init() {
      ip_.read(reinterpret_cast<char*>(&nchar_), sizeof(int));
      std::string sid_str(nchar_, ' ');
      ip_.read(&sid_str[0], nchar_);
      ip_.read(reinterpret_cast<char*>(&(n_steps_)), sizeof(int));
      ip_.read(reinterpret_cast<char*>(&(n_posit_)), sizeof(int));
      if (print_) {
        std::cout << sid_str << "\n";
        printf("n_steps: %d\n",n_steps_);
        printf("n_posit_: %d\n",n_posit_);
      }
      sid_posit_ = StringToSID(sid_str);
      switch (sid_posit_) {
        case (SID::filament): 
          n_misc_=3;
          break;
        default:
          break;
      }
      ReadIteration();
    }

    // Reads next iteration of data from the posit file
    bool ReadIteration() {
      if (ip_.eof()) {
        return false;
      }
      ip_.read(reinterpret_cast<char*>(&n_objs_), sizeof(n_objs_));
      //printf("n_objs_: %d\n",n_objs_);
      if (n_objs_ > max_objs_) 
        ReallocateObjs();
      for (int i=0; i<n_objs_; ++i) {
        ip_.read(reinterpret_cast<char*>(&position_[ix(i,0)]), 3*sizeof(double));
        ip_.read(reinterpret_cast<char*>(&scaled_position_[ix(i,0)]), 3*sizeof(double));
        ip_.read(reinterpret_cast<char*>(&orientation_[ix(i,0)]), 3*sizeof(double));
        ip_.read(reinterpret_cast<char*>(&diameter_[i]), sizeof(double));
        ip_.read(reinterpret_cast<char*>(&length_[i]), sizeof(double));
        for (int j=0; j<n_misc_; ++j) {
          ip_.read(reinterpret_cast<char*>(&misc_[misc_ix(i,j)]), sizeof(double));
        }
      }
      return true;
    }

    void ReadAll() {
      bool still_reading = ReadIteration();
      do {
        PrintAll();
        still_reading = ReadIteration();
      } while (still_reading);
    }

    void SetPrint(bool p) {
      print_ = p;
    }

  public:

    // Reads number of objects, and positions/orientations, etc from posit file
    // and copies it to appropriate arrays
    // Note that all arrays are single-dimension, contiguous arrays
    bool GetNext(int *n_objs, double *pos, double *spos, double *u, double *d, double *l) {
      if (!ReadIteration()) {
        //printf("EOF reached\n");
        return false;
      }
      GetNObjs(n_objs);
      GetPosit(pos, spos, u, d, l);
      return true;
    }

    bool ReadNext() {
      return ReadIteration();
    }

    // Returns just the number of objects from the posit file
    // (for allocation of arrays, for example)
    void GetNObjs(int *n_objs) {
      *n_objs = n_objs_;
    }

    // Updates just the positions/orientations etc,
    void GetPosit(double *pos, double *spos, double *u, double *d, double *l) {
      std::copy(position_, position_+3*n_objs_, pos);
      std::copy(scaled_position_, scaled_position_+3*n_objs_, spos);
      std::copy(orientation_, orientation_+3*n_objs_, u);
      std::copy(diameter_, diameter_+n_objs_, d);
      std::copy(length_, length_+n_objs_, l);
    }
    void GetDiameter(double *d) {
      std::copy(diameter_, diameter_+n_objs_, d);
    }
    void GetLength(double *l) {
      std::copy(length_, length_+n_objs_, l);
    }
    void GetMisc(double *misc) {
      std::copy(misc_, misc_+n_misc_*n_objs_, misc);
    }

    // Print entire posit file to screen, for debugging
    void PrintPosit(std::string file_name) {
      SetPrint(true);
      OpenPosit(file_name);
      Init();
      ReadAll();
      ClosePosit();
      DeallocateObjs();
      SetPrint(false);
    }
    int NPosit() {
      return n_posit_;
    }
    int NSteps() {
      return n_steps_;
    }
    void LoadFile(std::string file_name) {
      OpenPosit(file_name);
      Init();
    }
    void CloseFile() {
      ClosePosit();
      DeallocateObjs();
      position_ = nullptr;
      n_misc_ = 0;
      max_objs_ = -1;
    }
    SID GetSID() {
      return sid_posit_;
    }
};
#endif // _SIMCORE_POSIT_READER_H_
