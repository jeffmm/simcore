#ifndef _SIMCORE_POSIT_READER_H_
#define _SIMCORE_POSIT_READER_H_
#include "auxiliary.h"

class PositReader {
  private:
    int nchar_,
        n_posit_,
        n_steps_,
        n_objs_,
        max_objs_ = 0;
    double *position_ = nullptr,
           *scaled_position_,
           *orientation_,
           *diameter_,
           *length_;
    bool print_ = false;
    std::fstream ip_;

    int ix(int x, int y) {
      return 3*x+y;
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

    void ReallocateObjs() {
      max_objs_ = n_objs_;
      DeallocateObjs();
      AllocateObjs();
    }

    void DeallocateObjs() {
      if (position_ != nullptr) {
        delete[] position_;
        delete[] scaled_position_;
        delete[] orientation_;
        delete[] diameter_;
        delete[] length_;
      }
    }
    void AllocateObjs() {
      position_ = new double[max_objs_*3];
      scaled_position_ = new double[max_objs_*3];
      orientation_ = new double[max_objs_*3];
      diameter_ = new double[max_objs_];
      length_ = new double[max_objs_];
    }

    void PrintAll() {
      if (!print_) return;
      for (int i=0; i<n_objs_; ++i) {
        printf("position_ %d: {%2.2f, %2.2f, %2.2f}\n",i+1,position_[ix(i,0)],position_[ix(i,1)],position_[ix(i,2)]);
        printf("scaled_position_ %d: {%2.2f, %2.2f, %2.2f}\n",i+1,scaled_position_[ix(i,0)],scaled_position_[ix(i,1)],scaled_position_[ix(i,2)]);
        printf("orientation_ %d: {%2.2f, %2.2f, %2.2f}\n",i+1,orientation_[ix(i,0)],orientation_[ix(i,1)],orientation_[ix(i,2)]);
        printf("diameter_ %d: %2.2f\n",i+1,diameter_[i]);
        printf("length_ %d: %2.2f\n",i+1,length_[i]);
      }
    }

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
    }

    bool ReadIteration() {
      if (ip_.eof()) return false;
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

    bool GetNext(int *n_objs, double *pos, double *spos, double *u, double *d, double *l) {
      if (!ReadIteration()) {
        //printf("EOF reached\n");
        return false;
      }
      *n_objs = n_objs_;
      std::copy(position_, position_+3*n_objs_, pos);
      std::copy(scaled_position_, scaled_position_+3*n_objs_, spos);
      std::copy(orientation_, orientation_+3*n_objs_, u);
      std::copy(diameter_, diameter_+n_objs_, d);
      std::copy(length_, length_+n_objs_, l);
      return true;
    }

    bool GetNObjs(int *n_objs) {
      if (ip_.eof()) return false;
      ip_.read(reinterpret_cast<char*>(&n_objs_), sizeof(n_objs_));
      if (n_objs_ > max_objs_) 
        ReallocateObjs();
      *n_objs = n_objs_;
      return true;
    }

    bool GetPosit(double *pos, double *spos, double *u, double *d, double *l) {
      if (ip_.eof()) return false;
      for (int i=0; i<n_objs_; ++i) {
        ip_.read(reinterpret_cast<char*>(&position_[ix(i,0)]), 3*sizeof(double));
        ip_.read(reinterpret_cast<char*>(&scaled_position_[ix(i,0)]), 3*sizeof(double));
        ip_.read(reinterpret_cast<char*>(&orientation_[ix(i,0)]), 3*sizeof(double));
        ip_.read(reinterpret_cast<char*>(&diameter_[i]), sizeof(double));
        ip_.read(reinterpret_cast<char*>(&length_[i]), sizeof(double));
      }
      std::copy(position_, position_+3*n_objs_, pos);
      std::copy(scaled_position_, scaled_position_+3*n_objs_, spos);
      std::copy(orientation_, orientation_+3*n_objs_, u);
      std::copy(diameter_, diameter_+n_objs_, d);
      std::copy(length_, length_+n_objs_, l);
      return true;
    }

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
    }
};
#endif // _SIMCORE_POSIT_READER_H_
