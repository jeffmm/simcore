// Implementation for nexp lookup tables

#include "lookup_table.h"

void LookupTable::Init(int pndim, std::vector<double> *x,
                       double (*func)(std::vector<double> &x ,void *params),
                       void *params) {
  ndim_ = pndim;

  x_ = new std::vector<double>[ndim_];
  for (int i = 0; i < ndim_; ++i) {
    x_[i] = x[i];
  }

  if (ndim_ == 2) {
    std::vector<double> xvec;
    xvec.resize(2);
    table_.resize(x_[0].size() * x_[1].size());
    for (int i = 0; i < x_[0].size(); ++i) {
      for (int j = 0; j < x_[1].size(); ++j) {
        xvec[0] = x_[0][i];
        xvec[1] = x_[1][j];
        table_[i * x_[1].size() + j] = func(xvec, params);
      }
    }
  }

  if (debug_trace) {
    printf("NExp lookup table ->\n");
    for (int i = 0; i < table_.size(); ++i) {
      printf("[%d] -> %2.8f\n", i, table_[i]);
    }
  }
}

void LookupTable::OutputBinary(std::string outfile) {
  FILE *f_out = std::fopen(outfile.c_str(), "w");
  
  std::fwrite(&ndim_, sizeof(int), 1, f_out);
  for (int i = 0; i < ndim_; ++i) {
    int size = x_[i].size();
    std::fwrite(&size, sizeof(int), 1, f_out);
  }
  for (int i = 0; i < ndim_; ++i)
    std::fwrite(x_[i].data(), sizeof(double), x_[i].size(), f_out);

  int n_linear = 1; for (int i = 0; i < ndim_; ++i) n_linear *= x_[i].size();

  std::fwrite(table_.data(), sizeof(double), n_linear, f_out);
  std::fclose(f_out);
}

double LookupTable::Lookup(double x, double y) {
  double xvec[2] = {x,y};
  return Lookup(xvec);
}

double LookupTable::Lookup(double *x) {
  int lowindex[2];
  double t[2];
  for (int i = 0; i < ndim_; ++i) {
    double L = x_[i].back() - x_[i].front();
    int n_bin = x_[i].size()-1;
    lowindex[i] = (int) (x[i]/L * n_bin);
    double pos = x[i];
    if (lowindex[i] >= n_bin) {
      lowindex[i] = n_bin - 1;
      pos = x_[i].back();
    }
    
    t[i] = (pos - x_[i][lowindex[i]])/(x_[i][lowindex[i]+1] - x_[i][lowindex[i]]);
  }
  double y1 = table_[lowindex[0] * x_[1].size() + lowindex[1]];
  double y2 = table_[(1 + lowindex[0]) * x_[1].size() + lowindex[1]];
  double y3 = table_[(1 + lowindex[0]) * x_[1].size() + lowindex[1] + 1];
  double y4 = table_[lowindex[0] * x_[1].size() + lowindex[1] + 1];
  double returnval = (1 - t[0]) * (1 - t[1]) * y1 + t[0] * (1 - t[1]) * y2 +
      t[0] * t[1] * y3 + (1 - t[0]) * t[1] * y4;

  return returnval;
}

double LookupTable::Invert(int dim, double u, double val[]) {
  int lowindex[ndim_];
  double t[ndim_];
  double lookupvec[ndim_];
  if (dim != 0) {
    std::cerr << "Invert only works in '0' dim for now." << std::endl;
  }
  for (int i = 0; i < ndim_; ++i) {
    if (i != dim) {
      // std::vector<double>::iterator
      //     low=std::lower_bound (x_[i].begin(), x_[i].end(), val[i]) - 1;
      // lowindex[i] = low - x_[i].begin();
      // if (lowindex[i] == - 1) {
      //     lowindex[i] = 0;
      //     low++;
      // }
      // else if (lowindex[i] == x_[i].size() - 1) {
      //     lowindex[i] = x_[i].size()-2;
      //     low--;
      // }
      double L = x_[i].back() - x_[i].front();
      int n_bin = x_[i].size()-1;
      int i_bin = static_cast<int> (val[i]/L * n_bin);
      if (i_bin >= n_bin)
          i_bin = n_bin - 1;
      lowindex[i] = i_bin;
      //printf("%g %g %g %d %d\n", x[i], i_bin * (L/n_bin), *low, i_bin, lowindex[i]);
      double *low = &x_[i][i_bin];
      
      t[i] = (val[i] - *low)/(*(low + 1) - *low);
    }
  }
  val[dim] = x_[dim][x_[dim].size()-1];
  double max = Lookup(val);

  int lowi = 0;
  int highi = x_[dim].size() - 1;

  double f_low, f_high;
  while (lowi <= highi) {
    lowindex[dim] = (highi+lowi)/2;
    val[dim] = x_[dim][lowindex[dim]];
    f_low  = Lookup(val)/max; //table_[lowindex[0] * x_[1].size() + lowindex[1]];
    val[dim] = x_[dim][lowindex[dim]+1];
    f_high = Lookup(val)/max; //table_[lowindex[0] * x_[1].size() + lowindex[1] +
    //((dim == 0) ? x_[1].size() : 1)];
    if (f_low < u) {
      if (f_high > u)
        break;
      lowi = lowindex[dim];
    }
    else if (f_high > u) {
      highi = lowindex[dim];
    }
  } 
  double y1 = table_[lowindex[0] * x_[1].size() + lowindex[1]];
  double y2 = table_[(1 + lowindex[0]) * x_[1].size() + lowindex[1]];
  double y3 = table_[(1 + lowindex[0]) * x_[1].size() + lowindex[1] + 1];
  double y4 = table_[lowindex[0] * x_[1].size() + lowindex[1] + 1];
  double a = (-(1 - t[1]) * y1 + (1 - t[1]) * y2 + t[1] * y3 -t[1] * y4)/max;
  double b = ((1-t[1]) * y1 + t[1] * y4)/max;
  val[dim] = (u - b) * (x_[dim][lowindex[dim]+1] - x_[dim][lowindex[dim]]) / a +
      x_[dim][lowindex[dim]];
  return val[dim];
}

void LookupTable::OutputInterpolatedGrid(std::vector<double> *x, std::string outfile) {
  std::vector<double> outtable;
  outtable.resize(x[0].size() * x[1].size());
  
  for (int i = 0; i < x[0].size(); ++i) {
    for (int j = 0; j < x[1].size(); ++j) {
      double myx[] = {x[0][i],x[1][j]};
      outtable[i * x[1].size() + j] = Lookup(myx);
    }
  }

  FILE *f_out = std::fopen(outfile.c_str(), "w");
  
  std::fwrite(&ndim_, sizeof(int), 1, f_out);
  int size = x[0].size();
  std::fwrite(&size, sizeof(int), 1, f_out);
  size = x[1].size();
  std::fwrite(&size, sizeof(int), 1, f_out);
  
  std::fwrite(x[0].data(), sizeof(double), x[0].size(), f_out);
  std::fwrite(x[1].data(), sizeof(double), x[1].size(), f_out);

  std::fwrite(outtable.data(), sizeof(double), outtable.size(), f_out);

  std::fclose(f_out);
}
