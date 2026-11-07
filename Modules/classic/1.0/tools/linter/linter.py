#!/usr/bin/env python3

'''
This tool accepts Fortran files as arguments, and modifies their source code to
fit the standards outlined for CLASSIC. When passed a list of Fortran files, it
will create a backup directory in the same location, copy the files into that
backup directory, and replaces their original location with linted, prettified
versions.

Fortran 77 files are converted to free-form Fortran 90/95. Not all errors are
corrected by the linter. For example, if multiple do-loops are terminated by
a single continue statement, the linter will insert text into the code that will
cause it to fail compilation, to prompt the developer to fix this issue.

All changes made are documented in the /backups/*.comments file, and are not
overwritten if the file is re-linted.

@author: Matthew Fortier
gitlab.com/mfortier
'''

import os
import sys
import argparse
import re
import shutil
from tools.prettifier import *
from tools.structure import *
from tools.fixer import *
from tools.whitespace import *
from tools.decapitalizer import *

def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input', help='Input file(s)', nargs='+')
    parser.add_argument('--restore', help='Removes linted files and replaces them with backups', \
                        action="store_true")
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    files = []
    for item in args.input:
        if os.path.isfile(item):
            files.append(item)
    valid_files = list(filter(lambda x: (x[-2:] in ['.f','.F']) or \
                               (x[-4:] in ['.f90','.F90']),files))
    if args.restore:
        restore_files(valid_files)
    else:
        lint_files(valid_files)
    sys.exit()

# Deletes all linted files, restores the backups into their place, and deletes
# the backup folder
def restore_files(files):
    for ifile in files:
        backup = re.sub(r'(.*\/)(.*)',r'\1backups/\2',ifile)
        if not os.path.isfile(backup):
            backup = backup[:-2]
            if not os.path.isfile(backup):
                continue
        backup_dir = ''
        directory = re.match(r'(.*)(?=\/)', ifile)
        if not directory:
            backup_dir += './backups'
        else:
            backup_dir += directory.group(1)
            backup_dir += '/backups'
        os.remove(ifile)
        os.remove(backup + '.comments')
        newlocation = re.sub(r'(.*\/)backups\/(.*)',r'\1\2',backup)
        shutil.move(backup, newlocation)
        if len(os.listdir(backup_dir)) == 0:
            os.rmdir(backup_dir)

# Main function for file linting
def lint_files(files):
    for ifile in files:
        create_backup_folder(ifile)
        fname = ifile

        # run fixer if file is Fortran 77 format
        if ifile[-2:] == '.f' or ifile[-2:] == '.F':
            fx = Fixer(ifile)
            fname = ifile[:-2] + '.f90' if ifile[-2:] == '.f' else ifile[:-2] + '.F90'
            backup_f(ifile)
        elif ifile[-4:] == '.f90' or ifile[-4:] == ".F90":
            backup_f90(fname)

        # run decapitalizer (assists prettifier to identify rule violations)
        dc = Decapitalizer(fname)

        # run prettifier (many miscellaneous rules)
        pf = Prettifier(fname)

        # run structural analyzer (implicit none, unterminated loops)
        sc = StructuralAnalyzer(fname)

        # run whitespace checker (alignment)
        ws = WhitespaceChecker(fname)

        critical = sc.comments                 # for critical comments from structural analyzer
        uncorrected = pf.uncorrected_comments  # for uncorrected line comments in the prettifier
        corrected = pf.corrected_comments      # for corrected line comments in the prettifier

        # write all activity to the .f90.comments file (if the file doesn't already exist)
        comment_name = re.sub(r'(.*\/)(.*)', r'\1backups/\2.comments', ifile)
        if os.path.isfile(comment_name) or os.path.isfile(re.sub(r'(.*\/)(.*)', r'\1backups/\2.comments', ifile[:-2])):
            pass
        else:
            with open(comment_name, "w") as f:
                f.write("=========================== CRITICAL ERRORS ===========================\n\n")
                for line in critical:
                    f.write(line[1])
                    f.write('\n')
                else:
                    f.write('\nNo critical errors found\n\n')
                f.write("====================== Uncorrected Line Comments ======================\n\n")
                for line in uncorrected:
                    f.write(line[1])
                    f.write('\n')
                f.write("\n\n======================= Corrected Line Comments =======================\n\n")
                for line in corrected:
                    f.write(line[1])
                    f.write('\n')
    print("Linting completed successfully")

# Helper function for directory creation
def create_backup_folder(ifile):
    path = ''
    directory = re.match(r'(.*)(?=\/)', ifile)
    if not directory:
        path += './backups'
    else:
        path += directory.group(1)
        path += '/backups'
    try:
        os.mkdir(path)
    except:
        pass

# Helper function for .f files
def backup_f(ifile):
    dest = re.sub(r'(.*\/)(.*)', r'\1backups/\2', ifile)
    if os.path.isfile(dest):
        pass
    else:
        try:
            shutil.move(ifile, dest)
        except:
            pass

# Helper function for .f90 files
def backup_f90(ifile):
    dest = re.sub(r'(.*\/)(.*)', r'\1backups/\2', ifile)
    if os.path.isfile(dest) or os.path.isfile(dest[:-2]):
        pass
    else:
        try:
            shutil.copy(ifile, dest)
        except:
            pass

if __name__ == '__main__':
    main()
