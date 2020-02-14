import re
import os
import yaml


class SimcoreRenamer:
    """Rename simcore parameter and output files and run_name-type parameters

    Replaces a substring in a parameter file and associated output files with a format
    specified by the user.  Creates a new parameter file with updated run_name,
    load_name, etc parameters, and renames all associated output files (e.g. spec files,
    analysis files, etc) using the replacement substring.

    By default, the renamer will display the renaming strategy and ask for confirmation
    before proceeding.

    Example:

    # Substring of original file name as a regex (e.g. replaces v000_r000)
    substring_to_replace = "v[0-9]+_r[0-9]+"

    # Python style formatter (e.g. v000_r000 -> pf0.1_sp005_lp0005)
    replacement = "pf{}_sp{:03d}_lp{:04d}"

    # List of yaml keys to extract contents from parameter file for formatter, with ':'
    # delimitting sub nodes (e.g. yaml_node['filament']['packing_fraction'] for 'pf{}'
    # formatter, etc)
    formatter_contents = ['filament:packing_fraction', 'soft_potential_mag',
    'filament:perlen_ratio']

    renamer = SimcoreRenamer(substring_to_replace, replacement, formatter_contents)
    renamer.rename("activeff_v000_r000_reload001_params.yaml")

    """

    def __init__(self, original_substring, replacement_string, format_contents=[]):
        """Initialize SimcoreRenamer with renaming options.

        Args:
        original_substring (str): Regex of substring to replace in original filenames
        (e.g. v[0-9]+_r[0-9]+)
        replacement_string (str): Python-style formatter string (e.g. 'foo{}_bar{:03d}')
        format_contents (list of str, optional): List of yaml keys to use for '{}' in
            replacement string. The delimitter ':' can be used to denote substrings
            (e.g. ['node', 'node:subnode']).  Defaults to no special formatting.

        """

        self.original = original_substring
        self.replacement = replacement_string
        self.formatter = format_contents
        self.regex = None
        self.new_substring = None
        num_curly = self.replacement.count("{")
        if len(self.formatter) != num_curly:
            raise ValueError(
                "Formatter contents should have as many replacement values "
                "as the number of formatters in the replacement string"
            )

    def rename(self, param_file_name, confirm=True):
        """Rename parameter file and associated output files using initialized renamer rules.

        Args:
        param_file_name (str): Full path name to parameter file (colocated with files to
            rename)
        confirm (bool, optional): Show rename strategy and ask for confirmation before
            renaming files. Default to True.

        Returns:
        None

        """

        if not os.path.exists(param_file_name):
            raise ValueError("Could not locate file", param_file_name)
        if param_file_name.find("_params.yaml") == -1:
            raise ValueError(
                "Parameter file must be named in the usual way: 'run-name_params.yaml'"
            )
        rel_path, raw_param_file_name = os.path.split(param_file_name)
        if rel_path == "":
            rel_path = "."
        prefix = raw_param_file_name[: raw_param_file_name.find("_params.yaml")]

        # Get files corresponding to this parameter file
        files = [
            file
            for file in os.listdir(rel_path)
            if (file.find(prefix) == 0 and file.find("params.yaml") == -1)
        ]

        # Read parameter file and find appropriate replacement contents
        try:
            self.set_rename_params(param_file_name)
        except ValueError as err:
            raise err

        if confirm:
            print("Renaming strategy:")
            new_param_file_name = re.sub(
                self.regex, self.new_substring, raw_param_file_name
            )
            print(
                "  Create new parameter file:\n   ",
                os.path.join(rel_path, new_param_file_name),
            )
            if len(files) > 0:
                print("  Rename files:")

        # Rename all the files (or print rename strategy)
        for file in files:
            self.rename_file(os.path.join(rel_path, file), confirm)

        if confirm:
            # Request user confirmation
            user = input("Proceed with renaming? (y/N) ")
            if user == "y" or user == "Y":
                print("Renaming files")
                self.rename(param_file_name, confirm=False)
            else:
                print("Aborting rename")
        else:
            # Create new parameter file with changes filename values
            with open(param_file_name, "r") as pfile:
                lines = pfile.readlines()
            new_param_file_name = re.sub(
                self.regex, self.new_substring, raw_param_file_name
            )
            with open(os.path.join(rel_path, new_param_file_name), "w") as pfile:
                for line in lines:
                    pfile.write(re.sub(self.regex, self.new_substring, line))

    def set_rename_params(self, param_file_name):
        """Set internal rename parameters, such as formatter parameter values.

        Loads parameter yaml file and finds associated values from keys listed formatter
            contents. Currently only looks at the first subspecies in a species subnode.

        Args:
        param_file_name (str): original parameter file used to find parameter values

        Returns:
        None

        """
        rep_map = []
        # Load yaml file and find keys given by formatter contents
        with open(param_file_name) as pfile:
            yfile = yaml.safe_load(pfile)
            for item in self.formatter:
                try:
                    item = item.split(":")
                    if len(item) > 1:
                        temp = yfile[item[0]]
                        #  TODO: allow this to work for multiple filament etc types,
                        #  right now we are only picking first subspecies in list
                        if isinstance(temp, list):
                            temp = temp[0]
                        rep_map.append(temp[item[1]])
                    else:
                        rep_map.append(yfile[item[0]])
                except Exception:
                    raise ValueError(
                        "Parameter file does not have key given by formatter",
                        param_file_name,
                        item,
                    )
            self.regex = re.compile(self.original)
            self.new_substring = self.replacement.format(*rep_map)

    def rename_file(self, fname, confirm=True):
        """Rename file fname using internal replacement rules."""
        path, file = os.path.split(fname)
        file = re.sub(self.regex, self.new_substring, file)
        new_fname = os.path.join(path, file)
        # Only print the renaming rules if we are confirming
        if confirm:
            print("   ", fname, "->", new_fname)
        else:
            os.rename(fname, new_fname)
