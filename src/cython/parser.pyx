import re
    

cdef class Token:
    def __init__(self,int kind, int data_kind,str value, int lineno, int column,**kargs):
        self.kind = kind
        self.value= value
        self.data_kind = data_kind
        self.lineno = lineno
        self.column = column

cdef class Lexer:
    def __init__(self):
        self.count = 0
        self.line_num = 0
    
    
    cdef Token get_current_token(self):
        return self.tokens[self.count]
    
    
    cdef get_next_token(self):
        self.count += 1
        if self.count <len(self.tokens):
            return self.tokens[self.count]
        else:
            raise ValueError('No end charactor for expression!')
            # return Token(END,-1,'end',-1,-1)
            # self.tokens[self.count]
            return False
           
    
    cdef void tokenize(self, str code):
        keywords = ['ABS', 'MIN', 'MAX', 'SQRT', '+', '-', '*', '/', '//', '%', '**', '(', ')', '[', ']', ':', 'DIM', ',',
                    'l', 'angle', 'k1', 'k2', 'k3', 'k4', 
                    'betax', 'alphax', 'gammax', 'betay', 'alphay', 'gammay', 'etax', 'etapx', 'nux', 'nuy', 'chromx', 'chromy',
                    'emitx', 'geo1','geo2','geo3', 'totalchromx', 'totalchromy', 'circumference']
        
        token_enum = {'ADD':ADD, 'SUB':SUB, 'MUL':MUL, 'DIV':DIV, 'POW':POW, 'MOD':MOD, 
                      "FLOOR":FLOOR, "NUMBER":NUMBER, "PROPERTY":PROPERTY, "ABS":ABS, "SQRT":SQRT, 
                      "MAX":MAX, "MIN":MIN, "LPAREN":LPAREN, "RPAREN":RPAREN, "LBRA":LBRA, "RBRA":RBRA, 
                      "NAME":NAME, "INDEX":INDEX, 'SLICE':SLICE,'DIM':DIM,'COMMA':COMMA, 'DOT':DOT}
        
        token_regex = re.compile(r"""
            (?P<PROPERTY>(?<=\.)\w+)
            |(?P<ABS>ABS)
            |(?P<SQRT>SQRT)
            |(?P<MAX>MAX)
            |(?P<MIN>MIN)
            |(?P<DIM>DIM)
            |(?P<COMMA>\,)
            |(?P<NAME>[A-Za-z][A-Za-z0-9]+)
            |(?P<NUMBER>[-+]?[\d]+\.?[\d]*(?:[Ee][-+]?[\d]+)?)
            |(?P<DOT>\.)
            |(?P<LPAREN>\()
            |(?P<RPAREN>\))
            |(?P<LBRA>\[)
            |(?P<RBRA>\])
            |(?P<ID>[A-Za-z]+)
            |(?P<FLOOR>\/\/)
            |(?P<POW>\*\*)
            |(?P<MOD>\%)
            |(?P<ADD>\+)
            |(?P<SUB>\-)
            |(?P<MUL>\*)
            |(?P<DIV>\/)
            |(?P<SLICE>\:)
            |(?P<SKIP>[ \t]+)
            |(?P<MISMATCH>.)
            """ ,re.X)
        cdef:
            int ikind, line_start = 0,data_kind=-1
            str value ,kind
        self.count = 0
        self.tokens = []
        for mo in token_regex.finditer(code):
            kind = mo.lastgroup
            # ikind = token_enum[ mo.lastgroup ]
            value = mo.group()
            column = mo.start() - line_start
            if kind in ('INDEX','NUMBER'):
                ikind = token_enum[ mo.lastgroup ]
            elif kind == 'NAME':
                ikind = token_enum[ mo.lastgroup ]
            elif kind in token_enum.keys() and value in KWD_INDEX.keys():
                ikind = token_enum[ mo.lastgroup ]
                data_kind = KWD
            elif kind in token_enum.keys() and value in TWS_INDEX.keys():
                ikind = token_enum[ mo.lastgroup ]
                data_kind = TWS
            elif kind in token_enum.keys() and value in LOC_INDEX.keys():
                ikind = token_enum[ mo.lastgroup ]
                data_kind = LOC
            elif kind in token_enum.keys() and value in GLB_INDEX.keys():
                ikind = token_enum[ mo.lastgroup ]
                data_kind = GLB
            elif kind in token_enum.keys():
                ikind = token_enum[ mo.lastgroup ]
            elif kind == 'SKIP':
                continue
            elif kind == 'MISMATCH':
                raise RuntimeError(f'{value!r} unexpected on line {self.line_num}')
            else:
                raise ValueError(f'{value!r} is not proper here!')
            self.line_num += 1
            # print(f'tokennize:{ikind:<5} {data_kind:<5} {value:<10} ')
            self.tokens.append( Token(ikind, data_kind, value, self.line_num, column))
        self.tokens.append( Token(END, data_kind, value, self.line_num, column))



cdef class Parser:  # 定义语法分析器的类
    def __cinit__(self, dict elems_index):
        self.lexer = Lexer()  # 接收词法分析器对象
        self.elems_index = elems_index
        self.isdatabase = False
    
    
    cdef void set_database(self, double** kwd_properties, double** tws_properties, double** loc_properties, double* glb_properties)nogil:
        self.kwd_properties = kwd_properties
        self.tws_properties = tws_properties
        self.loc_properties = loc_properties
        self.glb_properties = glb_properties
        self.isdatabase = True
        
        
    cdef void error(self):  # 增加语法分析器的错误方法
        raise Exception('警告：错误的输入内容！')

        
    cdef void eat(self, int kind):
        if self.current_token.kind == kind:
            self.current_token = self.lexer.get_next_token()
            # if self.lexer.get_next_token() else Token(END, -1, 'end', -1, -1)
        else:
            self.error()

            
    cdef AST* parameter(self):
        cdef: 
            AST* node
            int index,parameter,data_kind
        current_token = self.current_token
        self.eat(NAME)
        self.eat(LBRA)
        name = current_token.value
        index = int(self.current_token.value)
        self.eat(NUMBER)
        self.eat(RBRA)
        self.eat(DOT)
        
        parameter = TOTAL_INDEX[self.current_token.value]
        data_kind = self.current_token.data_kind
        index = self.elems_index[name][index]
        self.eat(PROPERTY)
        node = new Parameter(index,parameter,data_kind, self.kwd_properties, self.tws_properties, self.loc_properties,self.glb_properties)
        return node
        
        
    cdef AST* factor(self):
        cdef: 
            AST* node
            Token token,current_token
        current_token = self.current_token
        if current_token.kind == NUMBER:
            self.eat(NUMBER)
            value = float(current_token.value)
            node = new Number(current_token.kind, value)
        elif current_token.kind == NAME:
            node = self.parameter()
        elif current_token.kind == MIN:
            self.eat(MIN)
            node = new Function( current_token.kind,self.factor() )
        elif current_token.kind == MAX:
            self.eat(MAX)
            node = new Function( current_token.kind,self.factor() )
        elif current_token.kind == ABS:
            self.eat(ABS)
            node = new Function( current_token.kind,self.factor() )
        elif current_token.kind == SQRT:
            self.eat(SQRT)
            node = new Function( current_token.kind,self.factor() )
        elif current_token.kind == DIM:
            self.eat(DIM)
            self.eat(LPAREN)
            node = self.factor()
            self.eat(COMMA)
            node = new Dim(node, DIM, self.expr())
            self.eat(RPAREN)
        elif current_token.kind == LPAREN:
            self.eat(LPAREN)
            node = self.expr()  # # 创建运算符节点对象
            self.eat(RPAREN)
        while self.current_token.kind in (POW, MOD, FLOOR, SLICE):
            token = self.current_token
            if token.kind == POW:
                self.eat(POW)
            elif token.kind == MOD:
                self.eat(MOD)
            elif token.kind == FLOOR:
                self.eat(FLOOR)
            elif token.kind == SLICE:
                self.eat(SLICE)
                node = new Slice(node,token.kind, self.term())
                return node
            node = new Node(node,token.kind,self.factor() )
        return node  # 返回运算符节点对象
    
    
    cdef AST* term(self):
        cdef: 
            AST* node
            Token token
        node = self.factor()  # 左侧节点对象
        while self.current_token.kind in (MUL, DIV):
            token = self.current_token
            if token.kind == MUL:
                self.eat(MUL)
            elif token.kind == DIV:
                self.eat(DIV)
            node = new Node(node, token.kind, self.factor())  # 创建运算符节点对象
        return node  # 返回节点对象


    cdef AST* expr(self):
        cdef: 
            AST* node
            Token token
        node = self.term()  # 左侧节点对象
        while self.current_token.kind in (ADD, SUB):
            token = self.current_token
            if token.kind == ADD:
                self.eat(ADD)
            elif token.kind == SUB:
                self.eat(SUB)
            node = new Node(node, token.kind, self.term())  #self.term 返回二元运算表达式的节点
        return node  # 返回运算符节点对象


    cdef AST* parse(self, str code):
        if not self.isdatabase:
            raise RuntimeError('Parser is not linked to database!')
        self.lexer.tokenize(code)
        self.current_token = self.lexer.get_current_token()
        return self.expr()

    
    
    
