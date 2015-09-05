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

from Lexer import Lexer
from odict import OrderedDict
from dnsqe import dnsqe_nice
import math, string

def comment_chomp(s):
    index = s.find("#")
    if index >= 0:
        return s[:index]

    index = s.find("!")
    if index >= 0:
        return s[:index]

    index = s.find("//")
    if index >= 0:
        return s[:index]

    return s

def pretty_variables(src):
    return_string = ''
    all_keys = list(src.keys())
    if len(all_keys) > 0:
        all_keys.sort(key=string.lower)
        max_key_length = max([len(key) for key in all_keys])
        for key in all_keys:
            return_string = return_string + key.ljust(max_key_length) + ' = ' + str(src[key]) + '\n'

    return return_string


class EquationLexer(Lexer):
    EQUATION_RULES = {
        "IDENTIFIER": r"(([a-zA-Z]\w*)((\')*\.([a-zA-Z0-9_]\w*))*)",
        "PLUS"      : r"\+",
        "MINUS"     : r"\-",
        "POW"       : r"\*\*|\^",
        "MUL"       : r"\*",
        "DIV"       : r"\/",
        "EQUALS"    : r"\=",
        "LPAREN"    : r"\(",
        "RPAREN"    : r"\)",
        "COMMA"     : r"\,",
        "VALUE"     : r"(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?",
    }

    def __init__(self):
        Lexer.__init__(self, EquationLexer.EQUATION_RULES)



class EquationParseError(Exception):
    pass


def has_attributes(s):
    s_parts = s.split('.')
    if len(s_parts) == 1:
        return False, []
    else:
        return True, s_parts[1:]


class EquationParser(object):

    def __init__(self, model):
        self.model = model

    def parse_arglist(self, token_reader, variable_list, function_call_list, constant_list):
        self.match(token_reader, 'LPAREN')
        out = '('
        next_token = token_reader.peek()
        if next_token == None:
            raise EquationParseError, 'Unexpected end of input, function call must be closed'

        if next_token[0] != 'RPAREN':
            #Dummy to let us into the loop for the first time
            next_token = ('COMMA', None)

            while next_token[0] == 'COMMA':
                out = out + self.parse_expression(token_reader,variable_list, function_call_list, constant_list)
                next_token = token_reader.peek()
                if next_token == None:
                    raise EquationParseError, 'Unexpected end of input'
                if next_token[0] == 'COMMA':
                    self.match(token_reader, next_token[0])
                    out = out + ','

        out = out + ')'
        self.match(token_reader, 'RPAREN')
        return out

    def parse_atom(self, token_reader, variable_list, function_call_list, constant_list):
        out = ''
        next_token = token_reader.peek()
        negated = False

        #If we have a negative sign, consume it.
        if next_token is not None and next_token[0] == 'MINUS':
            out = '-('
            self.match(token_reader, next_token[0])
            next_token = token_reader.peek()
            negated = True

        if next_token == None:
            raise EquationParseError, 'Unexpected end of input, need value, identifier or function call'

        #Begin a nested expression
        if next_token[0] == 'LPAREN':
            out = out + '('
            self.match(token_reader, next_token[0])
            out = out + self.parse_expression(token_reader, variable_list, function_call_list, constant_list)
            self.match(token_reader, 'RPAREN')
            out = out + ')'

        elif next_token[0] == 'IDENTIFIER':
            #check if this variable contains derivatives
            #if does, we replace each instance of ' with _prime
            #and preface the variable with _
            #Normal variables cannot start with _
            id = next_token[1]

            if id.find("'")>=0:
                id = '_' + id.replace("'", '_prime')


            self.match(token_reader, next_token[0])

            next_token = token_reader.peek()
            if next_token != None and next_token[0] == 'LPAREN':

                #This is a function call .. add to the function call list
                function_call_list.append(id)
                clobbered_id = self.model.addFunction(id)

                out = out + clobbered_id + self.parse_arglist(token_reader, variable_list, function_call_list, constant_list)

            else:
                #Just a variable, add to the variable list
                variable_list.append(id)
                #Get the clobbered id of the variable
                clobbered_id = self.model.addIdentifier(id)
                out = out + clobbered_id

        elif next_token[0] == 'VALUE':
            out = out + next_token[1]
            self.match(token_reader, next_token[0])
            constant_list.append(next_token[1])
        else:
            raise EquationParseError, 'Unexpected token %s:%s found' % (next_token[0], next_token[1])

        if negated:
            out = out + ')'

        return out

    def parse_factor(self, token_reader, variable_list, function_call_list, constant_list):
        out = self.parse_atom(token_reader, variable_list, function_call_list, constant_list)
        next_token = token_reader.peek()
        while next_token != None and (next_token[0] == 'POW'):
            out = 'math.pow(' + out + ','
            self.match(token_reader, next_token[0])
            out1 = self.parse_atom(token_reader, variable_list, function_call_list, constant_list)
            out = out + out1 + ')'
            next_token = token_reader.peek()

        return out

    def parse_term(self, token_reader, variable_list, function_call_list, constant_list):
        out = self.parse_factor(token_reader, variable_list, function_call_list, constant_list)
        next_token = token_reader.peek()
        while next_token != None and (next_token[0] == 'MUL' or next_token[0] == 'DIV'):
            self.match(token_reader, next_token[0])
            op = next_token[1]
            out = out + op + '('
            out1 = self.parse_factor(token_reader, variable_list, function_call_list, constant_list)
            out = out + out1 + ')'
            next_token = token_reader.peek()

        return out


    def parse_expression(self, token_reader, variable_list, function_call_list, constant_list):
        out = self.parse_term(token_reader, variable_list, function_call_list, constant_list)
        next_token = token_reader.peek()
        while next_token != None and (next_token[0] == 'PLUS' or next_token[0] == 'MINUS'):
            self.match(token_reader, next_token[0])
            op = next_token[1]
            out = out + op + '('
            out1 = self.parse_term(token_reader, variable_list, function_call_list, constant_list)
            out = out + out1 + ')'
            next_token = token_reader.peek()

        return out

    def parse_equation(self, token_reader):
        lhs_variable_list = []
        lhs_function_call_list = []
        lhs_constant_list = []

        rhs_variable_list = []
        rhs_function_call_list = []
        rhs_constant_list = []

        LHS = self.parse_expression(token_reader, lhs_variable_list, lhs_function_call_list, lhs_constant_list)
        self.match(token_reader, 'EQUALS')
        RHS = self.parse_expression(token_reader, rhs_variable_list, rhs_function_call_list, rhs_constant_list)


        #This is a parameter or an initial value
        if len(lhs_variable_list) == 1 and len(rhs_variable_list) == 0 and len(lhs_constant_list) == 0 and '-' not in LHS:
            param_id = lhs_variable_list[0]
            self.model.addParameter(param_id, RHS)

        elif len(rhs_variable_list) == 1 and len(lhs_variable_list) == 0 and len(rhs_constant_list) == 0 and '-' not in RHS:
            param_id = rhs_variable_list[0]
            self.model.addParameter(param_id, LHS)

        else:
            for x in lhs_variable_list:
                if not self.model.Parameters.has_key(x):
                    self.model.addVariable(x)

            for x in rhs_variable_list:
                if not self.model.Parameters.has_key(x):
                    self.model.addVariable(x)

            self.model.addEquation(LHS, RHS)


    def match(self, token_reader, token_type, token_value=None):
        try:
            tok = token_reader.next()
            if tok[0] != token_type:
                raise EquationParseError, 'Expecting token %s, found %s instead' % (token_type, tok[0])
            if token_value!=None and tok[1]!=token_value:
                raise EquationParseError, 'Expecting %s, found %s instead' % (token_value, tok[1])
        except StopIteration:
            raise EquationParseError, 'Expecting %s, end of input found instead' % (token_type)

        #print 'ACCEPT:', tok
        return tok



class EquationSenseError(Exception):
    pass


class Model(object):

    def __init__(self):
        self.Identifiers = OrderedDict()
        self.Parameters = OrderedDict()
        self.Variables = OrderedDict()
        self.Equations = OrderedDict()
        self.DifferentialEquations = OrderedDict()
        self.Functions = OrderedDict()
        self.InitialValues = OrderedDict()
        self.Bounds = {}

        #---------Some builtin functions-------
        self.Functions['exp'] = [('math.exp')]
        self.Functions['log'] = [('math.log')]
        self.Functions['log10'] = [('math.log10')]
        self.Functions['pow'] = [('math.pow')]
        self.Functions['sqrt'] = [('math.sqrt')]


        self.Functions['sin'] = [('math.sin')]
        self.Functions['cos'] = [('math.cos')]
        self.Functions['tan'] = [('math.tan')]
        self.Functions['hypot'] = [('math.hypot')]

        self.Functions['asin'] = [('math.asin')]
        self.Functions['acos'] = [('math.acos')]
        self.Functions['atan'] = [('math.atan')]


        self.Functions['radians'] = [('math.radians')]
        self.Functions['degrees'] = [('math.degrees')]


    def addInitialValue(self, name, value):
        self.InitialValues[name] = value

    def addIdentifier(self, name):
        if name not in self.Identifiers.keys():
            self.Identifiers[name] = 'id_' + str(len(self.Identifiers)).zfill(3)

        return self.Identifiers[name]

    def addVariable(self, id):
        attrs_present, parts = has_attributes(id)
        if attrs_present:
            raise EquationSenseError, 'Cannot use attributes in equation'

        if not self.Variables.has_key(id):
            self.Variables[id] = 0.0

        return 'x[' + str(self.Variables.index(id)) + ']'

    def addEquation(self, LHS, RHS, eqn_id = None):
        if eqn_id is None:
            eqn_id = '#' + str(len(self.Equations) + 1)

        self.Equations[eqn_id] = (LHS, RHS)

    def addVariableBound(self, id, lower_bound, upper_bound):
        existing = self.Bounds.get(id, None)

        if existing:
            if lower_bound is not None:
                existing = (float(lower_bound), existing[1])
            if upper_bound is not None:
                existing = (existing[0], float(upper_bound))

            self.Bounds[id] = existing
        else:
            if lower_bound is not None:
                lower_bound = float(lower_bound)
            if upper_bound is not None:
                upper_bound = float(upper_bound)

            self.Bounds[id] = (lower_bound, upper_bound)


    def addParameter(self, id, value):

        attrs_present, parts = has_attributes(id)
        if attrs_present:
            actual_id = id.split('.')[0]
            parts = [p.lower() for p in parts]

            if len(parts) == 1 and 'guess' in parts:
                self.addInitialValue(actual_id, value)
            elif len(parts) == 1 and 'min' in parts:
                self.addVariableBound(actual_id, value, None)
            elif len(parts) == 1 and 'max' in parts:
                self.addVariableBound(actual_id, None, value)
            else:
                raise EquationSenseError, 'Only supports guess, min and max attributes'

        else:
            #Check if the id is in the variable map, if it is delete it
            #and place it in the parameters map

            #This happens when a parameters is defined after its use
            #in an equation
            if self.Variables.has_key(id):
                del self.Variables[id]

            self.Parameters[id] = value

    def addFunction(self, func_id):
        if not self.Functions.has_key(func_id):
            raise EquationSenseError, '%s is not a valid function in this context.' % func_id

        return self.Functions[func_id][0]

    def generate_nle(self):
        #Make sure number of variables == number of equations
        if len(self.Variables) != len(self.Equations):
            raise EquationSenseError, 'Number of equations does not equal number of variables (No least squares solutions yet!)'

        line_count = 1
        #print self.Identifiers

        #Debugging only
        #output = 'def f_nle(n, x, f, iflag): #%d\n\tprint x\n' % (line_count) #for dnsqe

        output = 'def f_nle(n, x, f, iflag): #%d\n' % (line_count) #for dnsqe

        line_count += 1
        #output = 'def f_nle(iflag, m, n, x, f, fjac, ldfjac):\n' #for d1nlse
        #Dump out all the parameters
        for param in self.Parameters.keys():
            id = self.Identifiers[param]
            value = self.Parameters[param]
            output = output +  '\t%s = %s #%s %d\n' % (id, value, param, line_count)
            line_count += 1

        #Dump out all the variables
        var_count = 0
        for var in self.Variables.keys():
            id = self.Identifiers[var]
            output = output + '\t%s = x[%d] #%s %d\n' % (id, var_count, var, line_count)
            var_count += 1
            line_count += 1

        #Dump out all the equations
        eqn_count = 0
        for eqn in self.Equations.keys():
            LHS, RHS = self.Equations[eqn]
            eqn_string = '(%s) - (%s)' % (LHS, RHS)
            output = output +  '\tf[%d] = %s;\n' % (eqn_count, eqn_string)

            #General partial derivatives -- Begin
            eqn_count += 1
            line_count += 1

        output = output + '\treturn True\n'
        #print output
        return output

    def generate_nle_initial_guess(self, default=1.0):
        n = len(self.Variables)
        x = [default]*n
        f = [-1.0]*n
        for initvalue in self.InitialValues.keys():
            if self.Variables.has_key(initvalue):
                index = self.Variables.index(initvalue)
                x[index] = eval(self.InitialValues[initvalue])

        bounds = None
        #Do we have any bounds?
        bounds_count = len(self.Bounds.keys())
        if bounds_count > 0:
            bounds = [None]*len(self.Variables)
            i = 0
            for v in self.Variables:
                bounds_item = self.Bounds.get(v, None)
                bounds[i] = bounds_item
                i += 1

        return x, f, bounds

    def solve(self, x, f, bounds=None):
        f_nle_string = self.generate_nle()
        exec(f_nle_string)

        x, f, fnorm, info_code, info_mesg = dnsqe_nice(f_nle, None, x, bounds=bounds)

        i = 0
        keys = self.Variables.keys()
        for i in xrange(0, len(x)):
            self.Variables[keys[i]] = x[i]

        self.fnorm = fnorm;
        #print 'Norm of function residual vector:', fnorm
        self.info = info_code, info_mesg
        #print 'Info:', info_code
        #print 'Message:', info_mesg
        ret = ''
        if self.info[0] == 1:
            ret = ret + "Most recent solution:\n"
            ret = ret + "---------------------\n"
            ret = ret + "Solution norm: " + str(self.fnorm) + "\n"
            ret = ret + pretty_variables(self.Variables)
            return '\n' + ret + '\n'
        else:
            ret = ret + self.info[1]
            ret = ret + "---------------------\n"
            ret = ret + "Solution norm: " + str(self.fnorm) + "\n"
            ret = ret + pretty_variables(self.Variables)
            return '\n' + ret + '\n'


def comment_chomp(s):
    index = s.find("#")
    if index >= 0:
        return s[:index]

    index = s.find("!")
    if index >= 0:
        return s[:index]

    index = s.find("//")
    if index >= 0:
        return s[:index]

    return s



if __name__ == '__main__':
    m = Model()
    lex = EquationLexer()
    parse = EquationParser(m)

    """
    #Colebrook test
    parse.parse_equation(lex.scan("e=0.001"))
    parse.parse_equation(lex.scan("Re=1000000"))
    parse.parse_equation(lex.scan("1/sqrt(f)=-2*log(e/3.71+2.51/(Re*sqrt(f)),10)"))
    #parse.parse_equation(lex.scan("f.guess = 1"))
    parse.parse_equation(lex.scan("f.min = 0"))


    #Generate the initial guess and residual vector
    x, f, bounds = m.generate_nle_initial_guess()
    #print bounds

    #Solve the model
    m.solve(x, f, bounds)
    #print x
    """

    """
    f = file('examples/ex21.txt','rt')
    for line in f:
        new_line = comment_chomp(line).strip()
        if len(new_line) > 0:
            parse.parse_equation(lex.scan(new_line))

    x, f, bounds = m.generate_nle_initial_guess()
    print m.solve(x, f, bounds)
    """

    """
    #Parse the model
    #Rosenbrock equations, Solution is all 1
    parse.parse_equation(lex.scan("1 - x1 = 0"))
    parse.parse_equation(lex.scan("10*(x2-x1^2) = 0"))
    parse.parse_equation(lex.scan("1 - x3 = 0"))
    parse.parse_equation(lex.scan("10*(x4-x3^2) = 0"))

    #Generate the initial guess and residual vector
    x, f, bounds = m.generate_nle_initial_guess(-1.0)

    #Solve the model
    m.solve(x, f)
    print x
    """

    #Simple flash
    parse.parse_equation(lex.scan("P_IPA = exp(88.813 - 8498.6/T - 9.0766*log(T))/1000"))
    parse.parse_equation(lex.scan("P_H2O = exp(73.649 - 7258.2/T - 7.3037*log(T))/1000"))
    parse.parse_equation(lex.scan("P = x_IPA*P_IPA + x_H2O*P_H2O"))
    parse.parse_equation(lex.scan("x_IPA + x_H2O = 1"))
    parse.parse_equation(lex.scan("T = 373"))
    parse.parse_equation(lex.scan("x_IPA = 0.5"))
    parse.parse_equation(lex.scan("P.guess = 200.0"))
    parse.parse_equation(lex.scan("P_IPA.guess = 100.0"))
    parse.parse_equation(lex.scan("P_H2O.guess = 100.0"))

    #Generate the initial guess and residual vector
    x, f, bounds = m.generate_nle_initial_guess(-1.0)

    #Solve the model
    m.solve(x, f)
    print x

    """
    parse.parse_equation(lex.scan("1 - x1 = 0"))
    #Generate the initial guess and residual vector
    x, f, bounds = m.generate_nle_initial_guess(-1.0)

    #Solve the model
    m.solve(x, f)
    print x
    """
    
    """
    parse.parse_equation(lex.scan("Keq = 3.6e13"))
    parse.parse_equation(lex.scan("Y = 0.01 - 2*X"))
    parse.parse_equation(lex.scan("4*(X^3)/(Y^3) = Keq"))

    #Generate the initial guess and residual vector
    x, f, bounds = m.generate_nle_initial_guess(0.1)

    #Solve the model
    m.solve(x, f)
    print x
    """

    """
    #Badly scaled exponential
    parse.parse_equation(lex.scan("0 = exp(21000/T)/(T**2 - 1.11e11)"))
    """


    """
    #Reaction term example
    parse.parse_equation(lex.scan("DH=-57798.0"))
    parse.parse_equation(lex.scan("GAMMA=0.283E-6"))
    parse.parse_equation(lex.scan("ALPHA=7.256"))
    parse.parse_equation(lex.scan("BETA= 2.298E-3"))
    parse.parse_equation(lex.scan("DH+TR*(ALPHA+TR*(BETA/2+TR*GAMMA/3)) -298.0*(ALPHA+298.0*(BETA/2+298.0*GAMMA/3))=0"))
    """


    """
    parse.parse_equation(lex.scan("A = 0.38969"))
    parse.parse_equation(lex.scan("B = 0.55954"))
    parse.parse_equation(lex.scan("A*B*(B*(1-x)^2-A*x^2)/(x*(A-B)+B)^2+0.14845 = 0"))
    """

    """
    #Reactor example - Hard to solve
    parse.parse_equation(lex.scan("CA0-CA-theta*k1*CA/(1+KA*CB)= 0"))
    parse.parse_equation(lex.scan("CB-CB0-(theta*k1*CA/(1+KA*CB)-theta*k2*CB+theta*k2p*CC) = 0"))
    parse.parse_equation(lex.scan("CC-CC0-theta*k2*CB+theta*k2p*CC = 0"))
    parse.parse_equation(lex.scan("85*(T-T0)+0.02*(T^2-T0^2)-((16000+3*T-0.002*T^2)*((CA0-CA)/CA0)+(30000+4*T-0.003*T^2)*CC/CA0) = 0"))
    parse.parse_equation(lex.scan("k1 = 4e6*exp(-60000/(8.314*T))"))
    parse.parse_equation(lex.scan("KA = 17*exp(-7000/(8.314*T))"))
    parse.parse_equation(lex.scan("k2 = 3e4*exp(-80000/(8.314*T))"))
    parse.parse_equation(lex.scan("k2p = 3e4*exp(-90000/(8.314*T))"))
    parse.parse_equation(lex.scan("T0 = 298"))
    parse.parse_equation(lex.scan("CA0 = 3"))
    parse.parse_equation(lex.scan("CB0 = 0"))
    parse.parse_equation(lex.scan("CC0 = 0"))
    parse.parse_equation(lex.scan("theta = 300"))
    """

    """
    #Another reactor example - Hard to solve
    parse.parse_equation(lex.scan("V*(-rA)/vo*(CAO-CA) = 1.0"))
    parse.parse_equation(lex.scan("V*(-rB)/(vo*(CBO-CB)) = 1.0"))
    parse.parse_equation(lex.scan("V*rC/vo*CC = 1.0"))
    parse.parse_equation(lex.scan("V*rD/vo*CD = 1.0"))
    parse.parse_equation(lex.scan("V*rE/vo*CE = 1.0"))
    parse.parse_equation(lex.scan("5000.0*(350.0-T) - 25.0*(20.0+40.0)*(T-300.0) + V*( -rA*20000.0 + 2.0*r2C*10000.0 + 5000.0*r3E) = 0.0"))
    parse.parse_equation(lex.scan("rA = 2.0*r1B"))
    parse.parse_equation(lex.scan("rB = r1B+2*r2C"))
    parse.parse_equation(lex.scan("rC = -3.0*r1B + r2C"))
    parse.parse_equation(lex.scan("rD = -r3E - r2C"))
    parse.parse_equation(lex.scan("rE = r3E"))
    parse.parse_equation(lex.scan("r1B = -k1B*CA*CB"))
    parse.parse_equation(lex.scan("r2C = -k2C*CC*CB^2"))
    parse.parse_equation(lex.scan("r3E = k3E*CD"))
    parse.parse_equation(lex.scan("k1B = 0.4*exp((20000.0/R)*((1/300.0)-(1/T)))"))
    parse.parse_equation(lex.scan("k2C = 10.0*exp((5000.0/R)*((1/310.0)-(1/T)))"))
    parse.parse_equation(lex.scan("k3E = 10.0*exp((10000.0/R)*((1/320.0)-(1/T)))"))
    parse.parse_equation(lex.scan("R = 1.987"))
    parse.parse_equation(lex.scan("V = 500.0"))
    parse.parse_equation(lex.scan("vo = 75/3.3"))
    parse.parse_equation(lex.scan("CAO = 1.1"))
    parse.parse_equation(lex.scan("CBO = 2.2"))
    """


    """
    #Very tough flash calculation - hard to converge even with good initial guess
    parse.parse_equation(lex.scan("0 =x1-z1/(1+alpha*(k1-1))"))
    parse.parse_equation(lex.scan("0 =x2-z2/(1+alpha*(k2-1))"))
    parse.parse_equation(lex.scan("0 =x1+x2-(y1+y2) - alpha"))
    parse.parse_equation(lex.scan("p1=10^(7.62231-1417.9/(191.15+t))"))
    parse.parse_equation(lex.scan("p2=10^(8.10765-1750.29/(235+t))"))
    parse.parse_equation(lex.scan("gamma2=10^(B*x1*x1/((x1+B*x2/A)^2))"))
    parse.parse_equation(lex.scan("gamma1=10^(A*x2*x2/((A*x1/B+x2)^2))"))
    parse.parse_equation(lex.scan("k1=gamma1*p1/760"))
    parse.parse_equation(lex.scan("k2=gamma2*p2/760"))
    parse.parse_equation(lex.scan("y1=k1*x1"))
    parse.parse_equation(lex.scan("y2=k2*x2"))
    parse.parse_equation(lex.scan("t=88.538"))
    parse.parse_equation(lex.scan("B=0.7"))
    parse.parse_equation(lex.scan("A=1.7"))
    parse.parse_equation(lex.scan("z1=0.2"))
    parse.parse_equation(lex.scan("z2=0.8"))
    """
