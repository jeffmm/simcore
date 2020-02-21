"""Useful functions for parsing/analyzing simcore output files

Functions include SpecReader and PositReader for reading simcore .spec and .posit
files as well as FileRenamer to rename simcore output files en masse corresponding to
parameters in params.yaml files.
"""

__version__ = "0.1"
__author__ = "Jeffrey M. Moore"

from struct import unpack, pack, calcsize
from mmap import ACCESS_READ, mmap
from pathlib import Path
import re
import yaml

import numpy as np
import pandas as pd


class SpecReader:
    """Convert simcore .spec files into human readable formats.

    Example:
        sreader = SpecReader("my_spec_file.spec")
        sreader.make_human_readable()
        params, specs = sreader.get_params_specs()

    Functions:
        make_human_readable
        get_params_specs
        convert_old_spec
        rm_human_readable
    """

    def __init__(self, spec_file):
        """Class initialization

        Args:
            spec_file (string or Path): Location of spec file

        """
        self._input_file = Path(spec_file)
        if not self._input_file.is_file():
            raise FileNotFoundError("Could not find file {}", self._input_file)
        self._output_file = Path(str(spec_file) + ".hr")
        self._bpos = 0

    def _write(self, data_array, fout):
        """TODO."""
        fout.write(str(data_array[0]))
        if len(data_array) > 1:
            for datum in data_array[1:]:
                fout.write(" " + str(datum))
        fout.write("\n")

    def _read(self, fmt, buffer):
        """TODO"""
        old_bpos = self._bpos
        self._bpos += calcsize(fmt)
        retval = unpack(fmt, buffer[old_bpos : self._bpos])
        if len(retval) == 1:
            return retval[0]
        else:
            return retval

    def _writeb(self, data_array, write_as, fout):
        """Write data_array to output file fout as binary

        Args:
            data_array (list): list of data to write
            write_as (str): format to write data e.g. "iid" for int int double
            fout (file): binary output file

        """
        if len(write_as) == 1:
            fout.write(pack(write_as, data_array))
        else:
            fout.write(pack(write_as, *data_array))

    def convert_old_spec(self):
        """Convert old filament spec file from binary to binary with prefix 'new_'"""
        out_file = Path(self._input_file.parent, "new_" + self._input_file.name)
        with open(self._input_file, "rb", 0) as f, mmap(
            f.fileno(), 0, access=ACCESS_READ
        ) as s, open(out_file, "wb+") as fout:
            header = self._read("iid", s)
            self._writeb(header, "iid", fout)
            while True:
                try:
                    num = self._read("i", s)
                    self._writeb(num, "i", fout)
                    for i in range(num):
                        # diameter, length, lp, friction_par, friction_perp, l_bond
                        diameter, length, lp = self._read("ddd", s)
                        self._bpos += calcsize("dd")
                        l_bond, n_bonds = self._read("di", s)
                        # rtail[3], rhead[3]
                        r_site = np.array(self._read("ddd", s))
                        self._bpos += calcsize("ddd")
                        # u[3] for each bond
                        self._writeb([diameter, length, l_bond, n_bonds + 1],
                                     "dddi", fout)
                        self._writeb(list(r_site), "ddd", fout)
                        for bond in range(n_bonds):
                            u_bond = np.array(self._read("ddd", s))
                            r_site += l_bond * u_bond
                            self._writeb(list(r_site), "ddd", fout)
                        self._writeb([lp, 0], "dB", fout)

                except Exception as e:
                    print("EOF reached:", e)
                    break
        self._bpos = 0

    def make_human_readable(self):
        """Convert spec file from binary to plaintext with extension '.hr'"""

        with open(self._input_file, "rb", 0) as f, mmap(
            f.fileno(), 0, access=ACCESS_READ
        ) as s, open(self._output_file, "w+") as fout:

            header = self._read("iid", s)
            self._write(
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
            num = self._read("i", s)
            details = self._read("ddd", s)
            n_sites = self._read("i", s)
            posits = list(self._read("ddd" * n_sites, s))
            perlen = self._read("d", s)
            self._bpos += calcsize("c")
            self._write(list(header) + [num] + list(details) + [perlen, n_sites], fout)
            self._write(["site_positions_x_y_z"], fout)
            for fil in range(1, num):
                self._bpos += calcsize("dddi")
                posits += list(self._read("ddd" * n_sites, s))
                self._bpos += calcsize("dc")
            self._write(posits, fout)
            for sample in range(1, n_samples):
                self._bpos += calcsize("i")
                posits = []
                for fil in range(num):
                    self._bpos += calcsize("dddi")
                    posits += list(self._read("ddd" * n_sites, s))
                    self._bpos += calcsize("dc")
                self._write(posits, fout)
        self._bpos = 0

    def get_params_specs(self):
        """Return spec file parameters and specs timeseries as dataframes

        Returns:
            (DataFrame, DataFrame): Pandas DataFrame of parameters from spec file header
                and Pandas DataFrame of specs timeseries.

        """
        if not self._output_file.is_file():
            raise FileNotFoundError(
                "Human readable file not found: {}".format(self._output_file)
            )
        params = pd.read_csv(self._output_file, delim_whitespace=True, nrows=1)
        specs = pd.read_csv(
            self._output_file, skiprows=3, delim_whitespace=True, header=None
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
        if self._output_file.is_file():
            self._output_file.unlink()
        else:
            raise FileNotFoundError(
                "Human readable file not found: {}".format(self._output_file)
            )


class PositReader:
    """Convert simcore .posit files into human readable formats.

    Example:
    >>> from simcore_tools import PositReader
    >>> preader = PositReader("my_posit_file.posit")
    >>> preader.make_human_readable()
    >>> params, posits = preader.get_params_posits()

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
