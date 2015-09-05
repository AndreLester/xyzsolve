#!/usr/bin/env python

#Copyright 2010 Kiran Pashikanti <kpatvt@gmail.com>. All rights reserved.
#
#Redistribution and use in source and binary forms, with or without modification, are
#permitted provided that the following conditions are met:
#
#   1. Redistributions of source code must retain the above copyright notice, this list of
#      conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above copyright notice, this list
#      of conditions and the following disclaimer in the documentation and/or other materials
#      provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ``AS IS'' AND ANY EXPRESS OR IMPLIED
#WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
#FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
#CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
#ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
#ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#The views and conclusions contained in the software and documentation are those of the
#authors and should not be interpreted as representing official policies, either expressed
#or implied, of Kiran Pashikanti <kpatvt@gmail.com>.

import os
import logging
from google.appengine.ext import webapp
from google.appengine.ext.webapp.util import run_wsgi_app
from google.appengine.ext.webapp import template
from Model import *
import string
from Lexer import *
from dnsqe import vector_denorm
import urllib
import traceback





class InputHandler(webapp.RequestHandler):
    
    
    def error_out(self, ex, eqn, line_count):
        if len(eqn) > 32:
            eqn = eqn[0:31]
        return '!' + str(ex) + '\n' + str(line_count-1) + '\n' + eqn + '\n'
    
    def solve_equation_set(self, eqn_string):
        ret = ''
        m = Model()
        lex = EquationLexer()
        parse = EquationParser(m)
        eqn_part = eqn_string.split('\n');
        line_count = 1
        for eqn in eqn_part:
            eqnx = comment_chomp(eqn.strip())
            line_count += 1
            if len(eqnx) <= 0:
                continue
            
            try:
                parse.parse_equation(lex.scan(eqnx))
            except (UnknownTokenError, EquationParseError, EquationSenseError), ex:
                
                return self.error_out(ex, eqnx, line_count)            
            
            
        
        if len(m.Variables) == 0:
            return 'Nothing to solve. All is well.\n'
        
        try:
            x, f, bounds = m.generate_nle_initial_guess()
            m.solve(x, f, bounds)
        except EquationSenseError, ex:
            return self.error_out(ex, eqnx, line_count)
        except (ArithmeticError, ValueError) :
            ret =  'The solver encountered an error in your equations:\n'
            ret = ret + "- You have supplied an initial guess that is not valid.\n"
            ret = ret + "- The solver is not good enough to avoid discontinuities (divide by zeros, etc) in your problem space.\n"
            ret = ret + "Most recent solution:\n"
            ret = ret + "---------------------\n"
            ret = ret + "Solution norm:" + str(vector_denorm(len(f), f)) + "\n"
            i = 0
            keys = m.Variables.keys()
            for i in xrange(0, len(x)):
                m.Variables[keys[i]] = x[i]
            ret = ret + pretty_variables(m.Variables)
            return '\n' + ret
        
        #All is well, solutions found.
        if m.info[0] == 1:
            ret = ret + "Most recent solution:\n"
            ret = ret + "---------------------\n"
            ret = ret + "Solution norm: " + str(m.fnorm) + "\n"
            ret = ret + pretty_variables(m.Variables)
            return '\n' + ret + '\n'
        else:
            ret = ret + m.info[1]
            ret = ret + "---------------------\n"
            ret = ret + "Solution norm: " + str(m.fnorm) + "\n"
            ret = ret + pretty_variables(m.Variables)
            return '\n' + ret + '\n'
            
        raise RuntimeException('@An unexpected error occured, you should report this.')

    def get(self):
        self.response.headers['Content-Type'] = 'text/plain'
        eqn_string = ''
        try:
            #eqn_string = urllib.unquote(self.request.get('value'))
            
            eqn_string = self.request.get('value')
            self.response.out.write(self.solve_equation_set(eqn_string))
        
        except Exception, ex:
            self.response.out.write("@An unexpected error occurred, you should report this.\n")
            logging.error(traceback.format_exc())
            logging.error(str(ex))
            logging.error(eqn_string)

application = webapp.WSGIApplication([('/input.py', InputHandler)], debug=True)

def main():
    run_wsgi_app(application)

if __name__ == "__main__":
    main()
    
    