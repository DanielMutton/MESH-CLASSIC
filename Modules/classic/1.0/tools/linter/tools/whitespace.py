'''
This tool combs through the code and keeps a tally of how many indentations
each line of code should have (and automatically corrects it). Not too much to
say about it besides the logic can get a little twisted. Each line is checked
twice; once before adjusting whitespace, and once after (to affect the next line).

Note: this is the only component to print errors to standard output. Two errors
can occur, and both involve the `self.levelNumber` variable. This tracks how many
blocks in the code should be. For example, being inside a subroutine will increment
self.levelNumber by 1, as will being inside an if statement.

At the end of the program, self.levelNumber should be 0. If not, this means there
is a code block that is not terminated properly, or perhaps a bug in the linter's
ability to detect them.

The other error is thrown if self.levelNumber ever drops below 0, which should not
happen. If this error arises, check the line number where it happens.
'''
import re

class WhitespaceChecker(object):

    def __init__(self, fname):
        self.fname = fname # name of the file being checked
        with open(self.fname, 'r') as f:
            lines = f.readlines()
        self.lines = lines # stores the lines of the pre-linted file
        self.fixedlines = [] # will store the linted lines to be put back into the file
        self.ignoreLines = 0 # tracks number of lines from the !ignoreLint(x) directive
        self.directive = False # if true, we're in a preprocessor directive (and won't touch anything)
        self.levelNumber = 0 # how many block indentations in are we?
        self.multiline_spacing = 0
        self.continuedIf = False # flag for if we're in a continued if statement (important edge case)
        self.checkWhiteSpace()
        if self.levelNumber != 0: # we haven't finished with 0 whitespace; this is a problem.
            print("Whitespace error in {}: finished with levelNumber {}".format(self.fname, self.levelNumber))
        with open(self.fname, 'w') as f:
            for line in self.fixedlines:
                if not line.endswith('\n'):
                    line += '\n'
                f.write(line)

    def checkWhiteSpace(self):
        continuationLine = False # flag for if we're continuing from a previous line with '&'
        subcall = 0 # tracks the number of spaces to indent subsequent lines of a multi-line subroutine call
        for i, line in enumerate(self.lines):
            if self.ignoreLines == 0 and self.levelNumber < 0: # big error! we've got negative whitespace.
                print("Whitepsace error in {} at line {}: levelNumber < 0".format(self.fname, i+1))
                return
            # ignore blank lines; just tack a newline into the fixedlines list, and move on
            if re.match(r'^\s*\n', line, re.IGNORECASE):
                self.fixedlines.append("\n")
                continue
            # check if there's another reason we should be ignoring this line
            if self.skippable(line):
                continue
            # check if we flagged the previous line as continuation
            elif continuationLine:
                newline = ""
                # are we continuing a multi-line subroutine call?
                if subcall > 0:
                    for x in range(subcall):
                        newline += " "
                    newline += line.strip()
                # otherwise, just use the line as-is
                else:
                    newline = line
                # end of continuation lines?
                if not re.match(r'^[^\n!]*&', line, re.IGNORECASE) and \
                   not re.match(r'^\s*!', line, re.IGNORECASE):
                    continuationLine = False
                    subcall = 0
                self.fixedlines.append(newline)
                continue
            # all other cases
            else:
                newline = ""
                # check that this line follows regular formatting
                parsed_line = re.match(r'^(\s*)(\d*)(\s*)(.*)\n', line, re.IGNORECASE)
                if not parsed_line:
                    self.fixedlines.append(line)
                    continue
                self.analyzeLine1(parsed_line.group(4))

                # Here we'll consider multiline comment descriptions of variables
                linspc = re.match(r'^([^!]*)!<', line)
                if linspc:
                    self.multiline_spacing = len(linspc.group(1))
                elif not re.match(r'^(\s+)(!!.*$)', line):
                    self.multiline_spacing = 0

                linspc = re.match(r'^(\s+)(!!.*$)', line)
                if linspc and self.multiline_spacing > 0:
                    for a in range(self.multiline_spacing):
                        newline += " "
                    newline += linspc.group(2)
                    self.fixedlines.append(newline)
                    continue

                # if we have a labelled line, we must account for that in the whitespace
                if parsed_line.group(2) != "":
                    newline += parsed_line.group(2)
                    newline += " "

                # add whitespace based on the levelNumber, and account for line labelling
                num_spaces = max(0, 2*self.levelNumber - len(newline))
                for x in range(num_spaces):
                    newline += " "
                if self.continuedIf:
                    newline += "    "
                # 'contains' should be shifted back by 2 spaces
                if re.match(r'\bcontains\b', parsed_line.group(4), re.IGNORECASE):
                    newline = newline[:-2]
                newline += parsed_line.group(4)
                self.fixedlines.append(newline)
                self.analyzeLine2(parsed_line.group(4))
                if re.match(r'^[^\n!]*&', line, re.IGNORECASE) and not self.continuedIf:
                    continuationLine = True
                    subfuncall = re.match(r'^([^\n!\"\']*)\b(call|subroutine)\b([^\n!\(]*)\(', line, re.IGNORECASE)
                    declaration = re.match(r'^([^\n!:\"\']*):: ', line, re.IGNORECASE)
                    arithmetic = re.match(r'^([^\n!=\"\']*)= ', line, re.IGNORECASE)
                    if subfuncall:
                        subcall = len(subfuncall.group(0))
                    elif declaration:
                        subcall = len(declaration.group(0))
                    elif arithmetic:
                        subcall = len(arithmetic.group(0))


    # checks if a line should be skipped by the linter
    def skippable(self, line):
        # are we still in a preprocessor directive?
        if self.directive:
            if re.match(r'^[ \t]*#endif', line, re.IGNORECASE):
                self.directive = False
            self.fixedlines.append(line)
            return True
        # other possible cases:
        else:
            # are we starting a preprocessor directive?
            if re.match(r'^[ \t]*#if', line, re.IGNORECASE):
                self.fixedlines.append(line)
                self.directive = True
                return True
            else:
                # is there an !ignoreLint(x) directive?
                lineSkip = re.match(r'^\s*!ignoreLint\((\d+)\)', line, re.IGNORECASE)
                if lineSkip:
                    self.ignoreLines = int(lineSkip.group(1))
                    self.fixedlines.append(line)
                    return True
                elif self.ignoreLines > 0:
                    self.ignoreLines -= 1
                    self.fixedlines.append(line)
                    return True
        # no reason has been found to skip the line.
        return False

    # checks if we should revert whitespace on this line (eg. 'end if')
    def analyzeLine1(self, line):
        if re.match(r'^[^!]*<<<', line, re.IGNORECASE): # ignore lines that the linter has already flagged as awful
            return
        elif re.match(r'^\s*\b(end|case|else)\b', line, re.IGNORECASE):
            self.levelNumber -= 1

    # checks if we should increase whitespace on the NEXT line (eg. 'do j=1,n')
    def analyzeLine2(self, line):
        if re.match(r'^[^!]*<<<', line, re.IGNORECASE): # ignore lines that the linter has already flagged as awful
            return
        elif re.match(r'^[^\n!\'"]*\b(module(?! procedure )|interface|function|program|type |do|subroutine|case|select|else|forall|where)\b', line, re.IGNORECASE) \
        and not re.match(r'^[^\n!]*\bend\b', line, re.IGNORECASE):
            self.levelNumber += 1
        elif re.match(r'^[^!\n\'"]*\bif\b[^!\n]*\bthen\b', line, re.IGNORECASE):
            self.levelNumber += 1
        elif re.match(r'^[^!\n\'"]*(?<!end )\bif\b[^!\n]*&', line, re.IGNORECASE):
            self.continuedIf = True
            return
        elif self.continuedIf and re.match(r'^[^!\n]*\bthen\b', line, re.IGNORECASE):
            self.levelNumber += 1
        if self.continuedIf == True and re.match(r'^[^!\n]*&', line, re.IGNORECASE):
            return
        self.continuedIf = False
