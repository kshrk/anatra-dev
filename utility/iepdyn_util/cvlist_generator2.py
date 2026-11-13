#!/usr/bin/env python3

######################################################################
# @package cvlist_generator
# @brief   CV file list generator for IEPDYN 
# @author  Kento Kasahara (KK)
#
# (c) Copyright 2026 Kento Kasahara. All rights reserved.
######################################################################
# usage: ./cvlist_generator.py
#        (No arguments are required)
#
# You will be asked
# Enter the root directory to search
# Enter the target extension (e.g., txt)
# Enter the maximum search directory depth
# Enter directory level counted from the deepest directory for grouping (default: 1)
# Did these files come from perturbed dynamics? (y/n)
# [if yes] Enter state numbers used for Rij calculation (space-separated, empty input allowed)
#       
# Enter fraction of files to be listed from each state (default: 1.0 (corresponding to all the files))
# Enter filename to save the file list
# [if yes] Enter filename to save state_id and Rij_flag


import os
import readline
from pathlib import Path
from collections import defaultdict
import math

def collect_files(root_dir, extension, max_depth):
    """
    Search files with the specified extension
    up to the given directory depth.

    Symbolic links to directories are also followed.
    """

    root = Path(root_dir).resolve()
    matched_files = []

    # Prevent infinite recursion through symbolic links
    visited = set()

    for dirpath, dirnames, filenames in os.walk(
            root,
            followlinks=True):

        current = Path(dirpath)

        # Skip directories already visited through another symlink
        real_current = current.resolve()
        if real_current in visited:
            dirnames[:] = []
            continue
        visited.add(real_current)

        # Calculate depth from the root directory
        rel_parts = current.relative_to(root).parts
        depth = len(rel_parts)

        if depth > max_depth:
            dirnames[:] = []
            continue

        for filename in filenames:
            if filename.endswith(extension):
                matched_files.append(current / filename)

    return sorted(matched_files)


def group_files(files, group_level):
    """
    Group files by directory level counted
    from the deepest directory.
    """

    grouped_files = defaultdict(list)

    for path in files:

        parent_parts = path.parent.parts

        # Determine grouping directory
        if len(parent_parts) >= group_level:
            group_dir = parent_parts[-group_level]
        else:
            group_dir = path.parent.name

        grouped_files[group_dir].append(path)

    return grouped_files


def select_fraction(file_list, fraction):
    """
    Select the last fraction of files.

    fraction = 1.0 -> all files
    fraction = 0.9 -> last 90%
    fraction = 0.1 -> last 10%
    """

    n_files = len(file_list)

    if fraction >= 1.0:
        return file_list

    n_select = max(1, math.ceil(n_files * fraction))

    return file_list[-n_select:]

def logged_input(prompt):

    user_input = input(prompt)

    with open("cvgen.log", "a") as f:
        f.write(f"{prompt}{user_input}\n")

    return user_input


def setup_tab_completion():

    readline.parse_and_bind("tab: complete")

    def complete_path(text, state):

        if state == 0:

            # Special handling for "~"
            if text == "~":
                matches = ["~/"]
            else:

                # Expand "~" for searching
                expanded = os.path.expanduser(text)

                # Determine search directory and filename prefix
                if text.endswith("/"):
                    search_dir = expanded
                    basename = ""
                    typed_dir = text
                else:
                    dirname = os.path.dirname(expanded)
                    basename = os.path.basename(expanded)

                    search_dir = dirname if dirname else "."
                    typed_dir = (
                        os.path.dirname(text) + "/"
                        if "/" in text else ""
                    )

                matches = []

                try:
                    with os.scandir(search_dir) as entries:
                        for entry in entries:

                            if entry.name.startswith(basename):

                                candidate = typed_dir + entry.name

                                if entry.is_dir():
                                    candidate += "/"

                                matches.append(candidate)

                except OSError:
                    pass

            complete_path.matches = sorted(matches)

        try:
            return complete_path.matches[state]
        except IndexError:
            return None

    readline.set_completer(complete_path)

    # Allow spaces in paths
    readline.set_completer_delims("\t\n")

def main():

    setup_tab_completion()

    print("=== File Search Utility ===\n")

    # (1) Root directory
    root_dir = logged_input("(1) Enter the root directory to search: ").strip()

    # (2) Extension
    extension = logged_input(
        "(2) Enter the target extension (e.g., txt): "
    ).strip()

    # (3) Maximum depth
    max_depth = int(logged_input(
        "(3) Enter the maximum search directory depth: "
    ).strip())

    print("\nSearching files...\n")

    matched_files = collect_files(
        root_dir,
        extension,
        max_depth
    )

    if not matched_files:
        print("No matching files were found.")
        return

    print("=== Matched Files ===\n")

    for path in matched_files:
        print(path)

    print("\n----------------------------------------\n")

    # (4) Grouping level
    group_level = int(logged_input(
        "(4) Enter directory level counted from the deepest directory for grouping "
        "(default: 1): "
    ).strip() or "1")

    grouped_files = group_files(
        matched_files,
        group_level
    )

#    print("\n=== Grouped Files ===")

    sorted_groups = sorted(grouped_files.keys())

    # Assign state numbers
    print("\n----------------------------------------\n")

    auto_state = logged_input(
        "Assign state numbers automatically? (y/n): "
    ).strip().lower()
    
    state_map = {}
    
    if auto_state == "y":
    
        for idx, group in enumerate(sorted_groups, start=1):
            state_map[group] = idx
    
    else:
    
        print("\nDetected groups:")
    
        for idx, group in enumerate(sorted_groups, start=1):
            print(f"  Group {idx}: {group}")

        while True:
        
            user_input = logged_input(
                "\nEnter state numbers for each group "
                "(e.g., 1 2 3 4 4): "
            ).strip()
        
            state_ids = [int(x) for x in user_input.split()]
        
            if len(state_ids) < len(sorted_groups):
                print(
                    f"Please enter at least {len(sorted_groups)} integers."
                )
                continue
        
            # Ignore excess state numbers
            state_ids = state_ids[:len(sorted_groups)]
        
            break

#        while True:
#    
#            user_input = logged_input(
#                "\nEnter state numbers for each group "
#                "(e.g., 1 2 3 4 4): "
#            ).strip()
#    
#            state_ids = [int(x) for x in user_input.split()]
#    
#            if len(state_ids) != len(sorted_groups):
#                print(
#                    f"Please enter exactly {len(sorted_groups)} integers."
#                )
#                continue
#    
#            break
    
        for group, state_id in zip(sorted_groups, state_ids):
            state_map[group] = state_id
    
    print("\n=== Grouped Files ===")
    
    for idx, group in enumerate(sorted_groups, start=1):
    
        print(
            f"\n--- Group {idx}: {group} "
            f"(State {state_map[group]}) ---"
        )
    
        for filepath in sorted(grouped_files[group]):
            print(filepath)

#    # Assign state numbers
#    state_map = {}
#
#    for idx, group in enumerate(sorted_groups, start=1):
#
#        state_map[group] = idx
#
#        print(f"\n--- State {idx}: {group} ---")
#
#        for filepath in sorted(grouped_files[group]):
#            print(filepath)

    print("\n----------------------------------------\n")

    # Ask whether files are perturbed
    perturbed = logged_input(
        "Did these files come from perturbed dynamics? (y/n): "
    ).strip().lower()

    rij_state_set = set()

    if perturbed == "y":

        print("\nAvailable states:")

        for group in sorted_groups:
            print(f"  {state_map[group]} : {group}")

        user_input = logged_input(
            "\nEnter state numbers used for Rij calculation "
            "(space-separated, empty input allowed): "
        ).strip()

        if user_input:
            rij_state_set = {
                int(x) for x in user_input.split()
            }

    print("\n----------------------------------------\n")

    # Fraction selection
    fraction = float(logged_input(
        "Enter fraction of files to be listed from each state "
        "(default: 1.0 (corresponding to all the files)): "
    ).strip() or "1.0")

    # Prepare output lines
    file_lines = []
    state_lines = []
    
    # Prepare group-wise output lines
    group_file_lines = []
    group_state_lines = []

    for group in sorted_groups:
    
        state_id = state_map[group]
    
        rij_flag = 1 if state_id in rij_state_set else 0
    
        selected_files = select_fraction(
            sorted(grouped_files[group]),
            fraction
        )
    
        # Group-wise output
        current_file_lines = []
        current_state_lines = []
    
        for filepath in selected_files:
    
            current_file_lines.append(str(filepath))
    
            current_state_lines.append(
                f"{state_id} {rij_flag}"
            )
    
        # Append to total lists
        file_lines.extend(current_file_lines)
        state_lines.extend(current_state_lines)
    
        # Append to group-wise lists
        group_file_lines.append(current_file_lines)
        group_state_lines.append(current_state_lines)


    print("\n=== Output Preview ===\n")

    for filepath, stateinfo in zip(file_lines, state_lines):
        print(filepath)
        print(stateinfo)

    print("\n----------------------------------------\n")

    # Ask output filename for file list
    filelist_name = logged_input(
        "Enter filename to save the file list: "
    ).strip()
    
    # Save total file list
    with open(filelist_name, "w") as f:
        for line in file_lines:
            f.write(line + "\n")
    
    print(f"File list saved to: {filelist_name}")
    
    # Save group-wise file lists
    for idx, lines in enumerate(group_file_lines, start=1):
    
        group_filename = f"{filelist_name}.{idx:04d}"
    
        with open(group_filename, "w") as f:
            for line in lines:
                f.write(line + "\n")
    
        print(f"Group {idx} file list saved to: {group_filename}")
    
    
    # Ask output filename for state/Rij data only if perturbed
    if perturbed == "y":
    
        statefile_name = logged_input(
            "Enter filename to save state_id and Rij_flag: "
        ).strip()
    
        # Save total state/Rij data
        with open(statefile_name, "w") as f:
            for line in state_lines:
                f.write(line + "\n")
    
        print(f"State/Rij data saved to: {statefile_name}")
    
        # Save group-wise state/Rij data
        for idx, lines in enumerate(group_state_lines, start=1):
    
            group_filename = f"{statefile_name}.{idx:04d}"
    
            with open(group_filename, "w") as f:
                for line in lines:
                    f.write(line + "\n")
    
            print(
                f"Group {idx} state/Rij data saved to: {group_filename}"
            )

#    # Ask output filename for file list
#    filelist_name = logged_input(
#        "Enter filename to save the file list: "
#    ).strip()
#
#    with open(filelist_name, "w") as f:
#        for line in file_lines:
#            f.write(line + "\n")
#
#    print(f"File list saved to: {filelist_name}")
#
#    # Ask output filename for state/Rij data only if perturbed
#    if perturbed == "y":
#
#        statefile_name = logged_input(
#            "Enter filename to save state_id and Rij_flag: "
#        ).strip()
#
#        with open(statefile_name, "w") as f:
#            for line in state_lines:
#                f.write(line + "\n")
#
#        print(f"State/Rij data saved to: {statefile_name}")


if __name__ == "__main__":
    main()
