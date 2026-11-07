'''
This is the structural analyzer for the linter. It is not yet very fleshed out.
The reason this is formatted differently is because it only checks for a few
specific rules (eg. every 'do' loop has an 'end do'), and has specific actions
to take if that rule is broken. Not as simple as a search-and-replace like the
prettifier or decapitalizer.

Besides the do -> end do checking, it also counts the number of program components
and compares it to the number of 'implicit none' statements. Critical errors
are put in the comments file if these do not match up.
'''

import re

class StructuralAnalyzer(object):
    def __init__(self, fname):
        with open(fname, 'r') as f:
            lines = f.readlines()
        self.lines = lines
        self.comments = []
        self.checkStructure()
        with open(fname, 'w') as f:
            for line in self.lines:
                f.write(line)

    def checkStructure(self):
        structs = []
        doStructs = []
        implicitNone = 0
        programComponents = 0
        ignoreLines = 0
        for i, line in enumerate(self.lines):
            # if a line is '!ignoreLint(x)', then the linter skips the following x lines
            lineSkip = re.match(r'^\s*!ignoreLint\((\d+)\)', line)
            if lineSkip:
                ignoreLines = int(lineSkip.group(1))
                continue
            elif ignoreLines > 0:
                ignoreLines -= 1
                continue
            elif re.match(r'^((?!\bend\b|!|\'|\").)*\bdo\b', line, re.IGNORECASE):
                doStructs.append((i+1, 'do'))
            elif re.match(r'^[ \t]*\bend\b[ \t]+\bdo\b', line, re.IGNORECASE):
                if len(doStructs) == 0:
                    structs.append((i+1, '!do'))
                else:
                    doStructs = doStructs[:-1]
            elif re.match(r'^[ \t]*\bimplicit none\b', line, re.IGNORECASE):
                implicitNone += 1
            elif re.match(r'^((?!\bend\b|!).)*\b(subroutine|module|function)\b', \
                line, re.IGNORECASE):
                programComponents += 1

        structs.extend(doStructs)
        structs.sort(key=lambda tup: tup[0])

        # TODO: make this specify which program component has no 'implicit none',
        # and mark it in the source code to cause a compiler error.
        if implicitNone < programComponents:
            self.comments.append((0, ('>>> implicit none must be present in every program component')))
            self.comments.append((0, ('>>>   (program components: {}, implicit none statements: {})' \
                .format(programComponents, implicitNone))))

        for item in structs:
            if item[1] == 'do':
                comment = self.lines[item[0]-1].find('!')
                if comment > 0:
                    self.lines[item[0]-1] = self.lines[item[0]-1].replace('!', '<<< unterminated do loop !')
                else:
                    self.lines[item[0]-1] = self.lines[item[0]-1].replace('\n', ' <<< unterminated do loop\n')
                self.comments.append((item[0], '>>> {}: unterminated do loop'.format(item[0])))
            elif item[1] == '!do':
                comment = self.lines[item[0]-1].find('!')
                if comment > 0:
                    self.lines[item[0]-1] = self.lines[item[0]-1].replace('!', '<<< "end do" has no beginning !')
                else:
                    self.lines[item[0]-1] = self.lines[item[0]-1].replace('\n', ' <<< "end do" has no beginning\n')
                self.comments.append((item[0], '>>> {}: "end do" has no beginning'.format(item[0])))
