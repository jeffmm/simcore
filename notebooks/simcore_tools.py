"""Useful functions for parsing/analyzing simcore output files

Functions include spec_reader and posit_reader for reading simcore .spec and .posit
files.
"""

__version__ = "0.1"
__author__ = "Jeffrey M. Moore"

from struct import unpack, calcsize
from mmap import ACCESS_READ, mmap
from pathlib import Path

import pandas as pd


class spec_reader:
    """Convert simcore .spec files into human readable formats.

    Example:
        sreader = spec_reader("my_spec_file.spec")
        sreader.make_human_readable()
        params, specs = sreader.get_params_specs()

    """

    def __init__(self, spec_file):
        """Class initialization

        Args:
            spec_file (string or Path): Location of spec file

        """
        self.input_file = Path(spec_file)
        if not self.input_file.is_file():
            raise FileNotFoundError("Could not find file {}", self.input_file)
        self.output_file = Path(str(spec_file) + ".hr")
        self.bpos = 0

    def write(self, data_array, fout):
        """TODO."""
        fout.write(str(data_array[0]))
        if len(data_array) > 1:
            for datum in data_array[1:]:
                fout.write(" " + str(datum))
        fout.write("\n")

    def read(self, fmt, buffer):
        """TODO"""
        old_bpos = self.bpos
        self.bpos += calcsize(fmt)
        retval = unpack(fmt, buffer[old_bpos : self.bpos])
        if len(retval) == 1:
            return retval[0]
        else:
            return retval

    def make_human_readable(self):
        """Convert spec file from binary to plaintext with extension '.hr'"""

        with open(self.input_file, "rb", 0) as f, mmap(
            f.fileno(), 0, access=ACCESS_READ
        ) as s, open(self.output_file, "w+") as fout:

            header = self.read("iid", s)
            self.write(
                [
                    "n_steps",
                    "n_spec",
                    "delta",
                    "n_filaments",
                    "diameter",
                    "length",
                    "bond_length",
                    "persistence_length",
                    "n_sites",
                ],
                fout,
            )

            n_samples = header[0] // header[1]
            num = self.read("i", s)
            details = self.read("ddd", s)
            n_sites = self.read("i", s)
            posits = list(self.read("ddd" * n_sites, s))
            perlen = self.read("d", s)
            self.bpos += calcsize("c")
            self.write(list(header) + [num] + list(details) + [perlen, n_sites], fout)
            self.write(["site_positions_x_y_z"], fout)
            for fil in range(1, num):
                self.bpos += calcsize("dddi")
                posits += list(self.read("ddd" * n_sites, s))
                self.bpos += calcsize("dc")
            self.write(posits, fout)
            for sample in range(1, n_samples):
                self.bpos += calcsize("i")
                posits = []
                for fil in range(num):
                    self.bpos += calcsize("dddi")
                    posits += list(self.read("ddd" * n_sites, s))
                    self.bpos += calcsize("dc")
                self.write(posits, fout)
        self.bpos = 0

    def get_params_specs(self):
        """Return spec file parameters and specs timeseries as dataframes

        Returns:
            (DataFrame, DataFrame): Pandas DataFrame of parameters from spec file header
                and Pandas DataFrame of specs timeseries.

        """
        if not self.output_file.is_file():
            raise FileNotFoundError(
                "Human readable file not found: {}".format(self.output_file)
            )
        params = pd.read_csv(self.output_file, delim_whitespace=True, nrows=1)
        specs = pd.read_csv(
            self.output_file, skiprows=3, delim_whitespace=True, header=None
        )
        n_sites = params.n_sites[0]
        n_fils = params.n_filaments[0]
        fil_labels = [
            i
            for sub in [
                ["fil{:03d}".format(i)] * 3 * n_sites for i in range(n_fils)
            ]
            for i in sub
        ]
        site_labels = [
            i
            for sub in [
                ["site{:03d}".format(i)] * 3 for i in range(n_sites)
            ] * n_fils
            for i in sub
        ]
        arrays = [fil_labels, site_labels, ["x", "y", "z"] * n_sites * n_fils]
        columns = pd.MultiIndex.from_arrays(
            arrays, names=["filament", "site", "coord"]
        )
        specs.columns = columns
        specs.index.name = "time"
        return params, specs

    def rm_human_readable(self):
        """TODO"""
        if self.output_file.is_file():
            self.output_file.unlink()
        else:
            raise FileNotFoundError(
                "Human readable file not found: {}".format(self.output_file)
            )


class posit_reader:
    """Convert simcore .posit files into human readable formats.

    Example:
        preader = posit_reader("my_posit_file.posit")
        preader.make_human_readable()
        params, posits = preader.get_params_posits()

    """

    def __init__(self, posit_file):
        """Class initialization

        Args:
            posit_file (string or Path): Location of .posit file

        """

        self.input_file = Path(posit_file)
        self.output_file = Path("{}.hr".format(posit_file))

    def write_to_output(self, data_array, file):
        """Write array of data to file as strings

        Args:
            data_array (array): Array of data to write to file
            file (file): Output file

        """
        file.write(str(data_array[0]))
        if len(data_array) > 1:
            for datum in data_array[1:]:
                file.write(" " + str(datum))
        file.write("\n")

    def read_from_binary(self, type_string, file):
        """Read bytes corresponding from type_string from file

        Args:
            type_string (str): String corresponding to types, e.g. "iidd" for 2 ints and
                2 doubles
            file (file): Input file

        Returns:
            (array) Buffer from file corresponding to type string

        """
        return unpack(type_string, file.read(calcsize(type_string)))

    def make_human_readable(self):
        """Convert posit file from binary into plaintext with extension '.hr'"""
        with open(self.input_file, "rb") as fin, open(self.output_file, "w+") as fout:
            # Read header: n_steps, n_posit, timestep (int int double)
            header = self.read_from_binary("iid", fin)

            # Write header
            self.write_to_output(["n_steps", "n_posit", "delta"], fout)
            self.write_to_output(header, fout)

            # Read body
            self.write_to_output(
                [
                    "n_particles",
                    "pos_x",
                    "pos_y",
                    "pos_z",
                    "spos_x",
                    "spos_y",
                    "spos_z",
                    "u_x",
                    "u_y",
                    "u_z",
                    "diameter",
                    "length",
                ],
                fout,
            )
            n_samples = header[0] // header[1]
            for sample in range(n_samples):
                # read number of particles (int), 4 bytes
                num = self.read_from_binary("i", fin)[0]
                # Read num_particles * posits, each of which are 11 doubles (88 bytes)
                # These correspond to avg_pos (x,y,z), scaled_avg_pos (x,y,z),
                # orientation (x,y,z), diameter, and length
                posits = self.read_from_binary("ddddddddddd" * num, fin)
                self.write_to_output([num] + list(posits), fout)

    def get_params_posit(self):
        """Return posit parameters and posit file

        Returns:
            (DataFrame, DataFrame) Pandas DataFrames of posit parameters and posit
                timeseries

        """
        if not self.output_file.is_file():
            raise FileNotFoundError(
                "Could not find human readable file:{}".format(self.output_file)
            )
        params = pd.read_csv(
            self.output_file, delim_whitespace=True, nrows=1
        )
        posits = pd.read_csv(
            self.output_file, skiprows=3, delim_whitespace=True, header=None
        )
        return params, posits

    def rm_human_readable(self):
        """Remove human readable file to save disk space"""
        if (self.output_file.is_file()):
            self.output_file.unlink()
        else:
            raise FileNotFoundError(
                "Could not find human readable file:{}".format(self.output_file)
            )
