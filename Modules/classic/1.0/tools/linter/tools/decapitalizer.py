'''
This tool is used to decapitalize key words in Fortran 90. The list is not
necessarily exhaustive, and can be added to / subtracted from at will.

To add a word, simply put the regular expression for that word into the 'words'
list. The '\b' on either side of the word is a regular expression entity denoting
a word boundary. For example, this ensures 'IF' is decapitalized, but 'SCIF_VAR'
is not.

If you want to modify this to, say, capitalize everything instead of decapitalize,
it is as simple as modifying a few function calls. First, in the 'decapitalization'
function, add the optional argument re.IGNORECASE to the re.match() function call.
Then, in the re.replace() call, change the lower() to upper(). The impact on
performance should be negligible.
'''

import re

class Decapitalizer(object):

    words = [r'\bPROGRAM\b', r'\bMODULE\b', r'\bSUBROUTINE\b', r'\bFUNCTION\b', r'\bCALL\b',
    r'\bIF\b', r'\bTHEN\b', r'\bELSE\b', r'\bELSEIF\b', r'\bEND\b', r'\bENDIF\b',
    r'\bENDDO\b', r'\bDO\b', r'\bWHILE\b', r'\bSELECT\b', r'\bCASE\b', r'\bSELECTCASE\b',
    r'\bCONTINUE\b', r'\bIMPLICIT NONE\b', r'\bINTEGER\b', r'\bREAL\b', r'\bCOMPLEX\b', r'\bUSE\b',
    r'\bONLY\b', r'\bLOGICAL\b', r'\bDEFAULT\b', r'\bRETURN\b', r'\bINTENT\(\s*IN\s*\)', r'\bINTENT\(\s*OUT\s*\)', r'\bINTENT\(\s*INOUT\s*\)',
    r'\.AND\.', r'\.OR\.', r'\.NOT\.', r'\.EQV\.', r'\.NEQV\.', r'\.TRUE\.', r'\.FALSE\.', r'\bFORMAT\b', r'\bWRITE\b',
    r'\bCONTAINS\b', r'\bSAVE\b', r'\bDATA\b', r'\bCHARACTER\b', ]

    def __init__(self, fname):
        self.fname = fname
        with open(self.fname, 'r') as f:
            lines = f.readlines()

        self.words = r'|'.join(self.words) # join all the 'words' into a giant regular expression
        fmt = dict(
            words=self.words
        )

        # ensure that we're not changing words found inside comments
        # note the placeholder `{words}`
        rx = r'^[^\n!]*({words})'

        # insert the joined self.words expression into `rx` and compile it
        self.rule = re.compile(rx.format(**fmt))
        self.lines = lines
        self.fixedlines = []

        self.ignoreLines = 0    # flag for the !ignoreLint() directive
        self.directive = False  # flag for preprocessor directives

        self.decapitalize()     # run the decapitalizer

        with open(self.fname, 'w') as f: # rewrite the file with the new output.
            for line in self.fixedlines:
                f.write(line)

    # main decapitalization algorithm
    def decapitalize(self):
        for line in self.lines:
            if self.skippable(line): # check if there's a reason we should be skipping this line
                continue
            else:
                newline = line
                found = True
                while found: # iterate over the line until everything is decapitalized
                    match = re.match(self.rule, newline)
                    if match:
                        newline = newline.replace(match.group(1), match.group(1).lower())
                    else:
                        found = False
                self.fixedlines.append(newline)

    # checks if a line should be skipped by the linter
    def skippable(self, line):
        if self.directive: # if we're currently in a preprocessor directive
            if re.match(r'^[ \t]*#endif', line, re.IGNORECASE): # last line of the directive
                self.directive = False
            self.fixedlines.append(line)
            return True
        else:
            if re.match(r'^[ \t]*#if', line, re.IGNORECASE): # the start of a preprocessor directive
                self.fixedlines.append(line)
                self.directive = True
                return True
            else:
                lineSkip = re.match(r'^\s*!ignoreLint\((\d+)\)', line, re.IGNORECASE) # check for !ignoreLint()
                if lineSkip:
                    self.ignoreLines = int(lineSkip.group(1))
                    self.fixedlines.append(line)
                    return True
                elif self.ignoreLines > 0:
                    self.ignoreLines -= 1
                    self.fixedlines.append(line)
                    return True
        return False
