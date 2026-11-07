'''
This tool contains the bulk of the rules enforced in the linter. Everything
from modernized comparative operators, to whitespace around if statements, can
be found here. The general structure is similar to the decapitalizer, with rules
begin laid out ahead of time before compilation.

The main difference is that each rule is specified as a 3-tuple. The first element
is a regular expression of the code to look for. The second (optional) element is
the expression to replace it with (typically using regex groups). The third element
is a comment to add to the linter's comment files for the particular file being linted.

If the second element is set to `None`, the linter will not replace the phrase,
but instead notify the user of something that must be fixed manually due to its
inherent complexity. Examples would be replacing multi-terminated do-loops.

Similar rule tuples are coupled together in lists for organizational purposes.

Code originally adapted from cphyc, though heavily modified. Some unused portions
of the code still need to be removed.

@author: cphyc
github of original author: https://github.com/cphyc
'''
import re

class FortranRules(object):
    rules = [
        (r'(\d+)(\s+)(\.\d)', r'\2\1\3', 'Erroneous whitespace'),
        (r'(\s+)\n$', r'\n', 'Whitespace at end of line'),
        [
            # spaces around operators
            (r'^([^!\n]*)([^!\s\*=])({operators})([^\s\*])', r'\1\2 \3 \4',
             'Missing spaces around operator'),
            (r'^([^!\n]*)([^!\s\*=])({operators})', r'\1\2 \3',
             'Missing space before operator'),
            (r'^([^!\n]*)({operators})([^\s\*])', r'\1\2 \3',
             'Missing space after operator')
        ],
        [
            # " :: "
            (r'(\S)::(\S)', r'\1 :: \2',
             'Missing spaces around separator'),
            (r'(\S)::', r'\1 ::',
             'Missing space before separator'),
            (r'::(\S)', r':: \1',
             'Missing space after separator'),
            (r'(real|character\([^\)\n]*\)|type\([^\)\n]*\)|complex|integer)(\s+)(?!function)(\w+)',
             r'\1\2:: \3', 'Missing :: ')
        ],
        #(r'(real|integer|complex|logical)(\s+)([^:*])', r'\1\2:: ', ' :: added into variable declaration'),
        [
            # old -> new comparative operators
            (r'\.(eq|EQ)\.', r'==', '.eq. -> =='),
            (r'\.(neq|NEQ|ne|NE)\.', r'/=', '.ne. -> /='),
            (r'\.(gt|GT)\.', r'>', '.gt. -> >'),
            (r'\.(ge|GE)\.', r'>=', '.ge. -> >='),
            (r'\.(lt|LT)\.', r'<', '.lt. -> <'),
            (r'\.(le|LE)\.', r'<=', '.le. -> <=')
        ],
        # One should write "this, here" not "this,here"
        (r'({punctuations})(\w)', r'\1 \2',
         'Missing space after punctuation'),

        # should use lowercase for reserved words
        (r'\b({types_upper})\b', None,
         'Unnecessary use of upper case'),

        # if (foo), ...
        (r'({structs})\(', r'\1 (',
         'Missing space before parenthesis'),

        # Keep lines shorter than 80 chars
        (r'^.{linelen_re}.+$', None, 'Line length > {linelen} characters'),

        # Convert tabulation to spaces
        (r'\t', '  ', 'Should use 2 spaces instead of tabulation'),

        # Fix "real*4" to "real(4)"
        (r'({types})\*(\w+)', r'\1(\2)', 'Use new syntax TYPE(kind)'),

        # Remove dead space at end of line
        (r'(^[^\n!]*?)(\s\s+)(&|then)', r'\1 \3', 'Removed dead space before line end'),
        (r'\)then\b', r') then', ')then => ) then'),
        # Fix "foo! comment" to "foo ! comment"
        (r'(\w)\!', r'\1 !', 'At least one space before comment'),


        # Fix "!bar" to "! bar"
        (r'\!(\w)', r'! \1', 'Exactly one space after comment'),

        # Fix "!bar" to "! bar"
        (r'(![!><])(\S)', r'\1 \2', 'Exactly one space after comment'),

        # Remove trailing ";"
        (r';\s*$', r'\n', 'Useless ";" at end of line'),

        [
            # Support preprocessor instruction
            (r'(?<!#)(end|END)(if|IF|do|DO|subroutine|SUBROUTINE|function|FUNCTION|program|PROGRAM|select|SELECT)', r'\1 \2',
             'Missing space after `end\''),
            (r'\belseif\b', r'else if', 'elseif -> else if')
        ],

        # Trailing whitespace
        (r'( \t)+$', r'', 'Trailing whitespaces'),

        # Clean up arguments inside function calls)
        (r'(^[^!"\'=\n]*\()\s+(\w)', r'\1\2', 'Unnecessary space in function call'),
        (r'(\w)\s+\)', r'\1)', 'Spacing after arguments'),

        # Kind should be parametrized
        (r'\(kind\s*=\s*\d\s*\)', None, 'You should use "sp" or "dp" instead'),

        # MPI
        (r'include ["\']mpif.h[\'"]', None,
         'Should use `use mpi_f08` instead (or `use mpi` if not available)'),

        # Do not use continue to end loops, do not number loops
        [
            (r'(?<![^\s])(do|DO)([ \t]*)(\d+)([ \t]+)(.)', r'do \5', 'loops should not be number-labelled'),
            (r'^(\s*)([^\s!]*)(\s*)(continue|CONTINUE)', r'   \3end do ! loop \2', '"continue" is deprecated. Use "end do"')
        ],
        # Continuation lines
        [
            (r'([^\s])(&)', r'\1 \2', 'Missing space before &'),
            (r'(^\s*)&', r'\1', 'No & needed at beginning of continued line')
        ]
    ]

    types = [r'real', r'character', r'logical', r'integer', r'function', r'subroutine', r'end', r'module',
            r'if', r'do', r'while', r'call']
    operators = [r'==', r'/=', r'<=', r'<(?!=)', r'>=', r'(?<!=)>(?!=)', r'\.and\.',
                 r'\.or\.', r'(?<![\d\.][eE])\+', r'(?<![\d\.][eE])-', r'\*\*', r'(?<!\*)\*(?!\*)', r'(?<![=<>\/])=(?![=>])']
    structs = [r'if', r'IF', r'select', r'SELECT', r'case', r'CASE', r'while', r'WHILE']
    punctuation = ['''',',''' '\)', ';'] # commas have been removed from punctuation for now

    def __init__(self, linelen=120):
        self.linelen = linelen
        operators_re = r'|'.join(self.operators)
        types_upper = r'|'.join(self.types).upper()
        types_re = r'|'.join(self.types)
        types_all = r'|'.join([types_upper, types_re])
        struct_re = r'|'.join(self.structs)
        punctuation_re = r'|'.join(self.punctuation)
        fmt = dict(     # groups of words/operators to be formatted and compiled into the rules
            operators=operators_re,
            types_upper=types_upper,
            types=types_all,
            structs=struct_re,
            punctuations=punctuation_re,
            linelen_re="{%s}" % self.linelen,
            linelen="%s" % self.linelen)

        newRules = []
        for rule in self.rules:
            newRules.append(self.format_rule(rule, fmt)) # compile and format the rules
        self.rules = newRules

    def get(self):
        return self.rules

    def format_rule(self, rule, fmt):
        if isinstance(rule, tuple): # single rule (not grouped in a list)
            rxp, replacement, msg = rule
            msg = msg.format(**fmt) if msg is not None else None
            regexp = re.compile(rxp.format(**fmt))
            return (regexp, replacement, msg)
        elif isinstance(rule, list): # recursively format grouped lists of rules
            return [self.format_rule(r, fmt) for r in rule]
        else:
            raise NotImplementedError

class Prettifier(object):
    def __init__(self, fname, linelen=120):
        self.filename = fname

        with open(self.filename, 'r') as f:
            lines = f.readlines()
        self.lines = lines

        self.rules = FortranRules(linelen=linelen)
        self.corrected_lines = [] # these will be written to a new file
        self.corrected_comments = [] # comments about bad lines that were changed
        self.uncorrected_comments = [] # comments about bad lines that were not changed
        self.runs = 1 # tracks how many passes we've performed
        self.errcount = 0
        self.modifcount = 0
        self.errors = []
        self.correctedErrors = 0

        # Check the lines
        self.check_lines()
        while self.correctedErrors > 0:
            self.correctedErrors = 0
            self.uncorrected_comments = []
            self.errors = []
            self.lines = self.corrected_lines.copy()
            self.corrected_lines = []
            self.check_lines()
            self.runs += 1
            if self.runs > 20:
                break
        print("Completed linting {} in {} runs.".format(fname, self.runs))

        self.uncorrected_comments.sort(key=lambda tup: tup[0])
        self.corrected_comments.sort(key=lambda tup: tup[0])

        # Write corrected lines to old filename
        with open(self.filename, 'w') as f:
            for line in self.corrected_lines:
                f.write(line)

    def check_lines(self):
        ignoreLines = 0
        for i, line in enumerate(self.lines):
            # if a line is '!ignoreLint(x)', then the linter skips the following x lines
            lineSkip = re.match(r'^\s*!ignoreLint\((\d+)\)', line)
            if lineSkip:
                self.corrected_lines.append(line)
                ignoreLines = int(lineSkip.group(1))
                continue
            if ignoreLines > 0:
                self.corrected_lines.append(line)
                ignoreLines -= 1
                continue
            meta = {'line': i + 1,
                    'original_line': line.replace('\n', ''),
                    'filename': self.filename}

            line, _ = self.check_ruleset(line, line, meta, self.rules.get())
            self.corrected_lines.append(line)

    def check_ruleset(self, line, original_line, meta, ruleset, depth=0):
        if isinstance(ruleset, tuple):
            rule = ruleset
            line, hints = self.check_rule(
                line, original_line, meta, rule)
        else:
            for rule in ruleset:
                line, hints = self.check_ruleset(
                    line, original_line, meta, rule, depth+1)
                # Stop after first match
                if hints > 0 and depth >= 1:
                    break

        return line, hints

    def check_rule(self, line, original_line, meta, rule):
        regexp, correction, msg = rule
        errs = 0
        hints = 0
        newLine = line
        for res in regexp.finditer(original_line):
            meta['pos'] = res.start() + 1
            hints += 1
            if correction is not None:
                newLine = regexp.sub(correction, newLine)

            meta['correction'] = newLine
            if msg is not None:
                if correction is None:
                    self.uncorrected_comments.append((int(meta["line"]), str(meta["line"]) + ":" + msg))
                else:
                    self.correctedErrors += 1
                    self.corrected_comments.append((int(meta["line"]), str(meta["line"]) + ":" + msg))

        return newLine, hints
