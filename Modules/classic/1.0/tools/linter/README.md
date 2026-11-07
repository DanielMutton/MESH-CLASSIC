# lint-fortran-90

This tool is made for linting Fortran source code files, and modifying them to fit the style guide outlined for [CLASSIC](https://cccma.gitlab.io/classic_pages/info/conventions/).

## Getting Started

Python 3 is required to use this software.

Rather than being installed as a package, this linter is meant to be a standalone tool. All libraries used are part of the Python standard ecosystem. General use is as follows:

```
python3 linter.py path/to/file1.f path/to/file2.f another/path/to/a/file.f90
```

Relative paths and absolute paths are both acceptable. If a large number of files are to be linted, it may be easier to call the program like this:

```
python3 linter.py path/to/many/files/* path/to/more/files/*
```

This will cause most shells (such as bash) to automatically unpack the wildcard in the command.

.f/F and .f90/F90 files are both supported, and .f/F files will be converted to .f90/F90 by the linter. This is to support the transition from fixed-form Fortran 70 to free-form Fortran 90/95.

After being called, a backup folder will be created in any directory where files were specified. All specified files are copied into the backup directory, and will be found beside their `.comments` file listing all changes that were made.

For example, if we were to call:

```
python3 tools/linter/linter.py src/APREP.f90
```

from the base CLASSIC repository, we would end up with the following directory structure:

```
src/APREP.f90
src/backups/APREP.f
src/backups/APREP.f.comments
```

## Ignore directive

There are some circumstances where the linter needs to be disabled for a certain block of code. Using `GOTO` statements with `continue` flags causes undefined behaviour for the linter, and sometimes the prettifier can hurt bare strings with its spacing enforcement. In cases like these, it is best to use the `!ignoreLint()` directive. On the line above the code you wish to have ignored, put `!ignoreLint(x)` on its own line. When the linter encounters this, it will ignore the next `x` lines following the directive. The only module that does not obey this rule is the `fixer.py` module.

## Restoring Original Files

In case of unintended program behaviour, or a desire to return to the original files, simply call the program as you usually would if linting, but specify the `--restore` argument before or after the input files

```
python3 tools/linter/linter.py src/APREP.f --restore
```

Continuing from our previous example, this command will delete the .f90 file, delete the .comments file, restore the .f file to its original location, and delete the backups directory (if there are no remaining backups in it).

## Further Development

Currently, the `linter.py` program calls helper programs defined in the `tools` subdirectory to do the linting. There are 5 tools in use right now:
* `decapitalizer.py` - As the name implies, this tool decapitalizes Fortran key words that it encounters in the files passed in
* `fixer.py` - This tool is responsible for converting fixed-form fortran (`.f`) to free-form fortran (`.f90`)
* `structure.py` - Ensures `implicit none` is used in every block of code, looks for loops that are terminated with `continue` (depricated syntax)
* `whitespace.py` - Whitespace corrector. This is the tool that needs modification most frequently
* `prettifier.py` - Handles all other linting tasks, and definitely the biggest module. Everything from proper spacing to line length
