#include <assert.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <sys/stat.h>
#include <yaml-cpp/yaml.h>

#define NINT(x) ((x) < 0.0 ? (int)((x)-0.5) : (int)((x) + 0.5))

class VanHove {
private:
  float **posits_;
  float ***hist_;
  float **hist1d_;
  float *time_;
  int *lag_times_;
  int *n_avg_;
  int *handedness_;

  bool handedness_analysis_ = false;
  bool handedness_similar_ = true;
  float system_size_ = -1;
  float inv_size_ = -1;
  float density_ = -1;
  // TODO Initialize the next few parameters from yaml file
  float equil_percentage_ = 0.1;
  float bin_resolution_ = 0.5;
  int n_lags_ = 30;
  int sample_freq_ = 100;
  int max_lag_ = 15000;
  int n_bins_1d_ = -1;
  int n_rows_ = -1;
  int n_filaments_ = -1;
  std::string file_root_;
  std::string input_file_;
  void LoadParams(const YAML::Node &node);
  void BinDist(int t, float *r1, float *r2, bool double_count);
  void VerifyFile() const;
  void ReadLineNumber();
  void ReadFileHeader(std::ifstream &in);
  void RunSelfAnalysis();
  void RunCollectiveAnalysis();
  void AverageHistogram();
  void WriteOutputHeader();
  void OutputHistogram();
  void ClearHistogram();
  bool CheckHandedness(int i, int j);

  std::ofstream output_;
  std::ofstream output1d_;
  VanHove(const VanHove &) = delete;            // non construction-copyable
  VanHove &operator=(const VanHove &) = delete; // non copyable

public:
  VanHove(const YAML::Node &node);
  ~VanHove();
  void RunAnalysis();
};

VanHove::VanHove(const YAML::Node &node) {
  LoadParams(node);
  VerifyFile();
  ReadLineNumber();

  std::ifstream input(input_file_);
  ReadFileHeader(input);
  printf("File has %d rows of data and %d filaments\n", n_rows_, n_filaments_);
  int n_equil = static_cast<int>(equil_percentage_ * n_rows_);
  n_rows_ -= n_equil;
  if (max_lag_ >= n_rows_) {
    max_lag_ = n_rows_ - 1;
  }
  lag_times_ = new int[n_lags_];
  n_avg_ = new int[n_lags_];
  float lag_dist = log(max_lag_) / (n_lags_ - 1);
  lag_times_[0] = 0;
  for (int i = 1; i < n_lags_; ++i) {
    lag_times_[i] = (int)exp(i * lag_dist);
    if (lag_times_[i] <= lag_times_[i - 1]) {
      lag_times_[i] = lag_times_[i - 1] + 1;
    }
  }
  time_ = new float[n_rows_];
  // (N * (N-1) / 2 min distances, with 2 variables for each
  posits_ = new float *[n_rows_];
  for (int i = 0; i < n_rows_; ++i) {
    posits_[i] = new float[2 * n_filaments_];
  }
  handedness_ = new int[n_filaments_];
  hist_ = new float **[n_lags_];
  for (int i = 0; i < n_lags_; ++i) {
    hist_[i] = new float *[n_bins_1d_];
    for (int j = 0; j < n_bins_1d_; ++j) {
      hist_[i][j] = new float[n_bins_1d_];
    }
  }
  hist1d_ = new float *[n_lags_];
  for (int i = 0; i < n_lags_; ++i) {
    hist1d_[i] = new float[2 * n_bins_1d_];
  }
  // time, then x, y, z, ux, uy, uz for each filament
  std::string line;
  std::string z; // dummy variable
  int t = -n_equil;
  while (std::getline(input, line)) {
    if (t >= 0) {
      std::stringstream ss(line);
      ss >> time_[t];
      for (int i = 0; i < n_filaments_; ++i) {
        ss >> posits_[t][2 * i];
        ss >> posits_[t][2 * i + 1];
        // Ignore third dimension and orientations
        for (int j = 0; j < 4; ++j) {
          ss >> z;
        }
      }
    }
    t++;
  }
  assert(t == n_rows_);
  input.close();
  if (handedness_analysis_) {
    input.open(file_root_ + ".handedness.analysis");
    std::getline(input, line); // header
    std::getline(input, line); // n_filaments, intrinsic_curvature
    std::getline(input, line); // header
    std::getline(input, line); // handednesses
    std::stringstream ss(line);
    for (int i = 0; i < n_filaments_; ++i) {
      ss >> handedness_[i];
    }
    input.close();
  }
}

VanHove::~VanHove() {
  for (int i = 0; i < n_rows_; ++i) {
    delete[] posits_[i];
  }
  for (int i = 0; i < n_lags_; ++i) {
    for (int j = 0; j < n_bins_1d_; ++j) {
      delete[] hist_[i][j];
    }
    delete[] hist_[i];
    delete[] hist1d_[i];
  }
  delete[] hist_;
  delete[] hist1d_;
  delete[] posits_;
  delete[] time_;
  delete[] lag_times_;
  delete[] n_avg_;
  delete[] handedness_;
  printf("Analysis complete\n");
}

void VanHove::BinDist(int t, float *r1, float *r2, bool double_count) {
  float ds[2];
  float dr_mag = 0;
  for (int i = 0; i < 2; ++i) {
    float s1 = inv_size_ * r1[i];
    float s2 = inv_size_ * r2[i];
    s1 -= NINT(s1);
    s2 -= NINT(s2);
    ds[i] = s2 - s1;
    ds[i] -= NINT(ds[i]);
    dr_mag += ds[i] * ds[i];
  }
  /* Rescale dr_mag to 0-1 */
  dr_mag = sqrt(2 * dr_mag);
  /* I've encountered that ds[i] = 0.5 exactly, within float precision,
     so make sure we bin this value properly if it occurs */
  int i = static_cast<int>(n_bins_1d_ * (ds[0] + 0.5));
  if (i == n_bins_1d_) {
    i = n_bins_1d_ - 1;
  }
  int j = static_cast<int>(n_bins_1d_ * (ds[1] + 0.5));
  if (j == n_bins_1d_) {
    j = n_bins_1d_ - 1;
  }
  hist_[t][i][j] += 1.0;

  i = static_cast<int>(2 * n_bins_1d_ * dr_mag);
  if (i == 2 * n_bins_1d_) {
    i = 2 * n_bins_1d_ - 1;
  }
  hist1d_[t][i] += 1.0;
  /* For collective analysis, count the other's min distance vector */
  if (double_count) {
    hist1d_[t][i] += 1.0;
    i = static_cast<int>(n_bins_1d_ * (-ds[0] + 0.5));
    if (i == n_bins_1d_) {
      i = n_bins_1d_ - 1;
    }
    j = static_cast<int>(n_bins_1d_ * (-ds[1] + 0.5));
    if (j == n_bins_1d_) {
      j = n_bins_1d_ - 1;
    }
    hist_[t][i][j] += 1.0;
  }
}

void VanHove::LoadParams(const YAML::Node &node) {
  file_root_ = node["run_name"].as<std::string>();
  // TODO generalize to more than one filament type
  YAML::Node fnode = node["filament"][0];
  std::string spec_name = fnode["name"].as<std::string>();
  file_root_ = file_root_ + "_filament_" + spec_name;
  input_file_ = file_root_ + ".msd.analysis";
  system_size_ = 2.0 * node["system_radius"].as<float>();
  if (fnode["randomize_intrinsic_curvature_handedness"].as<bool>() &&
      (fnode["intrinsic_curvature"].as<float>() != 0 ||
       fnode["radius_of_curvature"].as<float>() > 0)) {
    handedness_analysis_ = true;
  }
  inv_size_ = 1.0 / system_size_;
  n_bins_1d_ = static_cast<int>(system_size_ / bin_resolution_);
  if (n_bins_1d_ % 2 == 0) {
    n_bins_1d_++;
  }
  if (node["n_dim"].as<int>() != 2) {
    printf("Error: VanHove2D called on 3D data\n");
    exit(1);
  } else if (node["n_periodic"].as<int>() != 2) {
    printf("Error: VanHove2D analysis must be with periodic boundary"
           " conditions\n");
    exit(1);
  }
}

void VanHove::VerifyFile() const {
  struct stat buffer;
  if (stat(input_file_.c_str(), &buffer) == 0) {
    printf("Reading from file %s\n", input_file_.c_str());
  } else {
    printf("Error: Failed to open MSD analysis file %s\n", input_file_.c_str());
    exit(1);
  }
}

void VanHove::ReadLineNumber() {
  int n_lines = 0;
  std::ifstream input(input_file_);
  std::string line;
  while (std::getline(input, line)) {
    ++n_lines;
  }
  input.close();
  n_rows_ = n_lines - 3;
  assert(n_rows_ > 0);
}

void VanHove::ReadFileHeader(std::ifstream &in) {
  std::string z; // dummy variable
  std::string line;
  // Skip first line
  std::getline(in, line);
  std::getline(in, line);
  std::stringstream ss(line);
  // Read nsteps, nspec, delta, filament number
  ss >> z >> z >> z >> n_filaments_;
  // Skip third line
  std::getline(in, line);
  density_ = n_filaments_ * inv_size_ * inv_size_;
}

void VanHove::RunAnalysis() {
  RunSelfAnalysis();
  RunCollectiveAnalysis();
  /* Run again on dissimilar handedness particles if we are distinguishing
     between similar and dissimilar handedness */
  if (handedness_analysis_) {
    handedness_similar_ = !handedness_similar_;
    RunCollectiveAnalysis();
  }
}

void VanHove::WriteOutputHeader() {
  output_ << "n_bins_1d n_frames\n";
  output_ << n_bins_1d_ << " " << n_lags_ << "\n";
  output_ << "lag_times\n";
  output_ << time_[lag_times_[0]] - time_[lag_times_[0]];
  for (int i = 1; i < n_lags_; ++i) {
    output_ << " " << time_[lag_times_[i]] - time_[lag_times_[0]];
  }
  output_ << "\n";
  output1d_ << "n_bins_1d n_frames\n";
  output1d_ << 2 * n_bins_1d_ << " " << n_lags_ << "\n";
  output1d_ << "lag_times\n";
  output1d_ << time_[lag_times_[0]] - time_[lag_times_[0]];
  for (int i = 1; i < n_lags_; ++i) {
    output1d_ << " " << time_[lag_times_[i]] - time_[lag_times_[0]];
  }
  output1d_ << "\n";
}

void VanHove::AverageHistogram() {
  for (int t = 0; t < n_lags_; ++t) {
    float factor = n_filaments_ * n_avg_[t];
    for (int i = 0; i < n_bins_1d_; ++i) {
      for (int j = 0; j < n_bins_1d_; ++j) {
        hist_[t][i][j] /= factor;
      }
    }
  }
  for (int t = 0; t < n_lags_; ++t) {
    float factor = n_filaments_ * n_avg_[t];
    for (int i = 0; i < 2 * n_bins_1d_; ++i) {
      hist1d_[t][i] /= factor;
    }
  }
}

void VanHove::OutputHistogram() {
  output_ << "n_samples_per_frame_per_filament\n";
  output_ << n_avg_[0];
  output1d_ << "n_samples_per_frame_per_filament\n";
  output1d_ << n_avg_[0];
  for (int i = 1; i < n_lags_; ++i) {
    output_ << " " << n_avg_[i];
    output1d_ << " " << n_avg_[i];
  }
  output_ << "\n";
  output1d_ << "\n";
  for (int t = 0; t < n_lags_; ++t) {
    for (int i = 0; i < n_bins_1d_; ++i) {
      for (int j = 0; j < n_bins_1d_; ++j) {
        if (j > 0)
          output_ << " ";
        output_ << hist_[t][i][j];
      }
      output_ << "\n";
    }
  }
  for (int t = 0; t < n_lags_; ++t) {
    for (int i = 0; i < 2 * n_bins_1d_; ++i) {
      if (i > 0)
        output1d_ << " ";
      output1d_ << hist1d_[t][i];
    }
    output1d_ << "\n";
  }
}

void VanHove::ClearHistogram() {
  for (int t = 0; t < n_lags_; ++t) {
    for (int i = 0; i < n_bins_1d_; ++i) {
      for (int j = 0; j < n_bins_1d_; ++j) {
        hist_[t][i][j] = 0.0;
      }
    }
  }
  for (int t = 0; t < n_lags_; ++t) {
    for (int i = 0; i < 2 * n_bins_1d_; ++i) {
      hist1d_[t][i] = 0.0;
    }
  }
}

void VanHove::RunCollectiveAnalysis() {
  fflush(stdout);
  std::string collective_output_name =
      file_root_ + ".van_hove_collective.analysis";
  std::string collective_output_name_1d =
      file_root_ + ".van_hove_collective_1d.analysis";
  if (handedness_analysis_ && handedness_similar_) {
    printf("Beginning Van Hove Collective Distribution analysis of %s for"
           " filaments with similar handedness\n",
           file_root_.c_str());
    collective_output_name =
        file_root_ + ".van_hove_collective_similar.analysis";
    collective_output_name_1d =
        file_root_ + ".van_hove_collective_similar_1d.analysis";
  } else if (handedness_analysis_) {
    printf("Beginning Van Hove Collective Distribution analysis of %s for"
           " filaments with dissimilar handedness\n",
           file_root_.c_str());
    collective_output_name =
        file_root_ + ".van_hove_collective_dissimilar.analysis";
    collective_output_name_1d =
        file_root_ + ".van_hove_collective_dissimilar_1d.analysis";
  } else {
    printf("Beginning Van Hove Collective Distribution analysis of %s\n",
           file_root_.c_str());
  }

  output_.open(collective_output_name);
  output1d_.open(collective_output_name_1d);
  WriteOutputHeader();
  std::fill(n_avg_, n_avg_ + n_lags_, 0);
  ClearHistogram();
  for (int t0 = 0; t0 < n_rows_ - sample_freq_; t0 += sample_freq_) {
    for (int t = 0; t < n_lags_; ++t) {
      int T = t0 + lag_times_[t];
      if (T > n_rows_ - 1)
        continue;
      n_avg_[t]++;
      for (int i = 0; i < n_filaments_ - 1; ++i) {
        float *r1 = &posits_[t0][2 * i];
        for (int j = i + 1; j < n_filaments_; ++j) {
          if (CheckHandedness(i, j))
            continue;
          float *r2 = &posits_[T][2 * j];
          BinDist(t, r1, r2, true);
        }
      }
    }
  }
  AverageHistogram();
  OutputHistogram();
  output_.close();
  output1d_.close();
}

void VanHove::RunSelfAnalysis() {
  printf("Beginning Van Hove Self Distribution analysis of %s\n",
         file_root_.c_str());
  fflush(stdout);
  std::string self_output_name = file_root_ + ".van_hove_self.analysis";
  std::string self_output_name_1d = file_root_ + ".van_hove_self_1d.analysis";
  output_.open(self_output_name);
  output1d_.open(self_output_name_1d);
  WriteOutputHeader();
  std::fill(n_avg_, n_avg_ + n_lags_, 0);
  ClearHistogram();
  for (int t0 = 0; t0 < n_rows_ - sample_freq_; t0 += sample_freq_) {
    for (int t = 0; t < n_lags_; ++t) {
      int T = t0 + lag_times_[t];
      if (T > n_rows_ - 1)
        continue;
      n_avg_[t]++;
      for (int i = 0; i < n_filaments_; ++i) {
        float *r1 = &posits_[t0][2 * i];
        float *r2 = &posits_[T][2 * i];
        BinDist(t, r1, r2, false);
      }
    }
  }
  AverageHistogram();
  OutputHistogram();
  output_.close();
  output1d_.close();
}

// 'continue' upon return true
bool VanHove::CheckHandedness(int i, int j) {
  if (handedness_similar_) {
    // Continue if handedness is dissimilar
    return (handedness_[i] != handedness_[j]);
  } else {
    // Continue if handedness is similar
    return (handedness_[i] == handedness_[j]);
  }
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("Missing <parameters yaml file> argument.\n");
    printf("Usage: %s <parameters yaml file>\n", argv[0]);
    return 1;
  }

  std::string fname(argv[1]);
  YAML::Node node = YAML::LoadFile(fname);
  VanHove vh(node);
  vh.RunAnalysis();
  return 0;
}
