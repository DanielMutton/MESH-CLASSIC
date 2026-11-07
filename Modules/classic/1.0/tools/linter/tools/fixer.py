'''
This tool is used to make the conversion from fixed-form to free-form Fortran.
You may notice the structure of this program is different from the other tools;
in part because I did not write the original code, and also because this tool
requires comparing adjacent lines to, for example, determine if one line is a
continuation from the previous line.

To that end, I have only changed this code when I've deemed it absolutely
necessary, preferring to modify the other tools to fill in the gaps.

Created on Jul 15, 2012
@author: jgoppert
github of original author: https://github.com/jgoppert
'''

import re

class Fixer(object):
    '''
    This class converts fixed format fortran code to free format.
    '''

    #misc regex's
    re_line_continuation = re.compile('^(\s\s\s\s\s)([^\s])([\s]*)')
    re_hollerith = re.compile('[\s,]((\d+)H)')
    re_number_spacing = re.compile('([^A-Za-z][\d.]+)[\s]([\d.]+[^A-Za-z])')
    re_f77_comment = re.compile('^[cC*]')
    re_f90_comment = re.compile('^!')

    #variable regex
    re_var_end = re.compile('.*[A-Za-z0-9]+$')
    re_var_begin = re.compile('(^\s\s\s\s\s)([^\s])([\s]*)([A-Za-z]+[0-9]?)+')

    #number regex
    re_number_end = re.compile('.*[\d.eE]+$')
    re_number_begin = re.compile('(^\s\s\s\s\s)([^\s])([\s]*)[\d.eE]+')

    #exponent regex
    re_exponent_old = re.compile('([\d.]E)[\s]([\d]+)')

    def __init__(self, input_filename):
        '''
        Constructor
        '''
        self.input_filename=input_filename
        self.output_filename=input_filename[:-2] + '.f90' if input_filename[-2:] == '.f' else input_filename[:-2] + '.F90'



        self.process()

    def process(self):

        #read source
        with open(self.input_filename, 'r') as f:
            source = f.readlines()

        #fix source code
        for i in range(len(source)):
            line = source[i]
            line = self.remove_new_line(line)
            line = self.fix_comment(line)
            line = self.fix_exponents(line)

            if self.is_line_continuation(line):
                prev_line_index = self.find_prev_line(source,i)
                prev_line = source[prev_line_index]
                line,prev_line = self.fix_line_continuation(line,prev_line)
                source[prev_line_index] = prev_line

            line = self.fix_number_spacing(line)
            source[i] = line

        #write new source
        with open(self.output_filename, 'w') as f:
            for line in source:
                f.write(line+'\n')

    def is_line_continuation(self, line):
        if len(line) > 5 and self.re_line_continuation.match(line):
            return True
        return False

    def find_prev_line(self, source, i):
        '''
        find first non-commented previous line
        '''
        if i>0:
            for j in range(i-1, -1, -1): # -1 since python uses [i-1,-1)
                if not self.is_f90_comment(source[j]):
                    return j
        return None

    def remove_new_line(self, line):
        return re.sub(r'[\r\n]','', line)

    def is_f77_comment(self, line):
        return self.re_f77_comment.match(line)

    def is_f90_comment(self, line):
        return self.re_f90_comment.match(line)

    def fix_comment(self, line):
        if self.is_f77_comment(line):
            line = re.sub('^[Cc*]','!', line, count=1)
        return line

    def fix_line_continuation(self, line, prev_line):
        continuation_type = self.find_continuation_type(line,prev_line)

        #if no line continuation found, return
        if not continuation_type:
            return (line,prev_line)

        #add appropriate continuation to current line
        if continuation_type == "hollerith":
            line = self.re_line_continuation.sub('\g<1>&\g<3>', line, count=1)
        elif continuation_type == "var" or continuation_type == "number" :
            line = self.re_line_continuation.sub('\g<1>\g<3>&', line, count=1)
        elif continuation_type == "generic":
            line = self.re_line_continuation.sub('\g<1>&\g<3>', line, count=1)
        else:
            raise Exception("unknown continuation type")

        #add line continuation to previous line
        prev_line = prev_line + '&'

        return (line,prev_line)

    def fix_exponents(self, line):
        return self.re_exponent_old.sub('\g<1>+\g<2>', line)

    def fix_number_spacing(self, line):
        return self.re_number_spacing.sub('\g<1>\g<2>', line)

    def find_continuation_type(self, line, prev_line):
        '''
        Finds line continuation type: hollerith, number, or generic
        if not a line continuation, return None
        '''
        if not self.is_line_continuation(line):
            return None

        # test for hollerith continuation
        for hollerith in self.re_hollerith.finditer(prev_line):
            hollerith_end = hollerith.end(1)+int(hollerith.group(2))+1
            prev_line_length = len(prev_line)-1
            #print "hollerith wrap detected"
            #print "\tline number:", i # note this is prev line number (i-1)+1 since starts at 0
            #print "\tline length:", prev_line_length
            #print "\thollerith end:", hollerith_end
            if (hollerith_end > prev_line_length):
#                print "hollerith wrap detected"
#                print "\tline number:", i # note this is prev line number (i-1)+1 since starts at 0
#                print "\tline length:", prev_line_length
#                print "\thollerith end:", hollerith_end
                return "hollerith"

        #var continuation
        if self.re_var_begin.match(line) and self.re_var_end.match(prev_line):
            return "var"

        #number continuation
        if self.re_number_begin.match(line) and self.re_number_end.match(prev_line):
            return "number"

        return "generic"
