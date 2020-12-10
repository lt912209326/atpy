import re

__all__=['Token','Lexer','Parser']

token_enum = {'ADD':ADD, 'SUB':SUB, 'MUL':MUL, 'DIV':DIV, 'POW':POW, 'MOD':MOD, 
              "FLOOR":FLOOR, "NUMBER":NUMBER, "PROPERTY":PROPERTY, "ABS":ABS, "SQRT":SQRT, 
              "MAX":MAX, "MIN":MIN, "LPAREN":LPAREN, "RPAREN":RPAREN, "LBRA":LBRA, "RBRA":RBRA, 
              "POSITION":POSITION, "INDEX":INDEX, 'SLICE':SLICE,'DIM':DIM,'COMMA':COMMA, 'DOT':DOT,'END':END}
enum_token = {value:key for key,value in token_enum.items()}

cdef class Token:
    def __init__(self,int kind, int data_kind,str value, int lineno, int column,**kargs):
        self.kind = kind
        self.value= value
        self.data_kind = data_kind
        self.lineno = lineno
        self.column = column

cdef class Lexer:
    def __init__(self,dict elems_index):
        self.count = 0
        self.line_num = 0
        self.elems_index = elems_index
    
    
    cdef Token get_current_token(self):
        return self.tokens[self.count]
    
    
    cdef get_next_token(self):
        self.count += 1
        if self.count <self.token_num:
            return self.tokens[self.count]
        else:
            raise ValueError('No end charactor for expression!')
            # return Token(END,-1,'end',-1,-1)
            # self.tokens[self.count]
            return False
    
    
    cdef int check_next_token(self):
        if self.count+1<self.token_num:
            return (<Token>self.tokens[self.count+1]).kind
        else:
            raise ValueError('Token error, out of memmory when indexing tokens !')
           
    
    cdef int tokenize(self, str code)except -1:
        keywords = ['ABS', 'MIN', 'MAX', 'SQRT', '+', '-', '*', '/', '//', '%', '**', '(', ')', '[', ']', ':', 'DIM', ',',
                    'l', 'angle', 'k1', 'k2', 'k3', 'k4', 
                    'betax', 'alphax', 'gammax', 'betay', 'alphay', 'gammay', 'etax', 'etapx', 'nux', 'nuy', 'chromx', 'chromy',
                    'emitx', 'geo1','geo2','geo3', 'totalchromx', 'totalchromy', 'circumference']
        
        
        token_regex = re.compile(r"""
            (?P<NAME>[A-Za-z][A-Za-z0-9]+)\[(?P<EID>\d+)\]
            |(?P<ABS>ABS)
            |(?P<SQRT>SQRT)
            |(?P<MAX>MAX)
            |(?P<MIN>MIN)
            |(?P<DIM>DIM)
            |(?P<COMMA>\,)
            |(?P<NUMBER>[\d]+\.?[\d]*(?:[Ee][-+]?[\d]+)?)
            |(?P<DOT>\.)
            |(?P<PROPERTY>[a-zA-Z0-9_]+)
            |(?P<LPAREN>\()
            |(?P<RPAREN>\))
            |(?P<LBRA>\[)
            |(?P<RBRA>\])
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
            elif kind == 'EID':
                ikind = token_enum['POSITION']
                if 'NAME' not in mo.groupdict().keys():
                    raise ValueError('Element name is not found!')
                name = mo.group('NAME')
                if name not in self.elems_index.keys():
                    raise ValueError('Element name is not found!')
                eid = int(mo.group('EID'))
                value = str(self.elems_index[name][eid])
            elif kind == 'PROPERTY':
                if value in KWD_INDEX.keys():
                    ikind = token_enum[ mo.lastgroup ]
                    value = str(KWD_INDEX[value])
                    data_kind = KWD
                elif value in TWS_INDEX.keys():
                    ikind = token_enum[ mo.lastgroup ]
                    value = str(TWS_INDEX[value])
                    data_kind = TWS
                elif value in LOC_INDEX.keys():
                    ikind = token_enum[ mo.lastgroup ]
                    value = str(LOC_INDEX[value])
                    data_kind = LOC
                elif value in GLB_INDEX.keys():
                    ikind = token_enum[ mo.lastgroup ]
                    value = str(GLB_INDEX[value])
                    data_kind = GLB
                else:
                    # assert value in TOTAL_INDEX.keys(), f'Unrecognized property {value}!'
                    raise ValueError(f'Unrecognized property {value}!')
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
        self.token_num = len(self.tokens)



cdef class Parser:  # 定义语法分析器的类
    def __cinit__(self, dict elems_index):
        self.lexer = Lexer(elems_index)  # 接收词法分析器对象
        self.elems_index = elems_index
        self.isdatabase = False
    
    
    cdef void set_database(self, double** kwd_properties, double** tws_properties, double** loc_properties, double* glb_properties)nogil:
        self.kwd_properties = kwd_properties
        self.tws_properties = tws_properties
        self.loc_properties = loc_properties
        self.glb_properties = glb_properties
        self.isdatabase = True
        
        
    cdef void error(self,int input, int current)except *:  # 增加语法分析器的错误方法
        raise ValueError(f'Error input:token {enum_token[input]} does not match input token {enum_token[current] }！')

        
    cdef void eat(self, int kind)except *:
        if self.current_token.kind == kind:
            self.current_token = self.lexer.get_next_token()
        else:
            self.error(kind,self.current_token.kind)
        

    cdef AST* function(self, int func)except NULL:
        cdef:
            int start,end, data_kind, index
            AST* node
            Token   current_token
        node =NULL
        self.eat(LPAREN)
        try:
            if self.lexer.check_next_token()==SLICE and func in (MIN, MAX):
                start = int(self.current_token.value)
                self.eat(POSITION)
                self.eat(SLICE)
                end = int(self.current_token.value)
                self.eat(POSITION)
                if start >=end:
                    raise ValueError('The first element should be frond of the second element in slice!')
                self.eat(COMMA)
                data_kind = self.current_token.data_kind
                index = int(self.current_token.value)
                self.eat(PROPERTY)
                # print('2: ',func, start, end, data_kind, index)
                # node =new Number(NUMBER,1.0)
                node = self.slice(func, start, end, data_kind, index)
            else:
                node = self.expr()
                current_token = self.current_token
                if current_token.kind == COMMA and func in (MIN, MAX,DIM):
                    self.eat(COMMA)
                    node = new BiFunction(node, func, self.expr() )
                elif current_token.kind == RPAREN and func in (ABS, SQRT):
                    node = new MonoFunction(func, node)
                else:
                    raise ValueError('Argument of function number is not correct!')
        except:
            raise ValueError('Argument of function number is not correct!')
        self.eat(RPAREN)
        return node
    
    
    cdef AST* slice(self, int func, int start, int end, int data_kind, int index)except NULL:
        cdef:
            AST* node
        if data_kind==KWD:
            node = new Slice(func, start, end, index, self.kwd_properties)
        elif data_kind==TWS:
            node = new Slice(func, start, end, index, self.tws_properties)
        elif data_kind==LOC:
            node = new Slice(func, start, end, index, self.loc_properties)
        else:
            raise ValueError('Error in parser.Parser.slice function!')
        return node
    
    
    cdef AST* property(self)except NULL:
        cdef: 
            AST*    node
            int     position=0, index, parameter,data_kind
            bint    isglb=True
        if self.current_token.kind==POSITION:
            position = int(self.current_token.value)
            self.eat(POSITION)
            self.eat(DOT)
            isglb = False
        if self.current_token.kind==PROPERTY:
            index = int(self.current_token.value)
            data_kind = self.current_token.data_kind
            self.eat(PROPERTY)
        else:
            raise ValueError('Property need to be specified!')
        if data_kind==KWD:
            if isglb:
                raise ValueError('Element position is needed for keyword property!')
            node = new Property(PROPERTY, position, index, &self.kwd_properties[position][index])
        elif data_kind==TWS:
            if isglb:
                raise ValueError('Element position is needed for twiss property!')
            node = new Property(PROPERTY, position, index, &self.tws_properties[position][index])
        elif data_kind==LOC:
            if isglb:
                raise ValueError('Element position is needed for local property!')
            node = new Property(PROPERTY, position, index, &self.loc_properties[position][index])
        elif data_kind==GLB:
            if isglb is not True:
                raise ValueError('Element position is not needed for global property!')
            node = new Property(PROPERTY, position, index, &self.glb_properties[index])
        else:
            raise ValueError('Error in parser.Parser.slice function!')
        return node
        
        
    cdef AST* factor(self)except NULL:
        cdef: 
            AST* node
            Token token,current_token
        current_token = self.current_token
        if current_token.kind == NUMBER:
            self.eat(NUMBER)
            value = float(current_token.value)
            node = new Number(current_token.kind, value)
        elif current_token.kind == SUB:
            self.eat(SUB)
            if self.current_token.kind==NUMBER:
                value = -1*float(self.current_token.value)
                node = new Number(self.current_token.kind, value)
                self.eat(NUMBER)
            else:
                node = new Number(NUMBER, 0.0)
                node = new Node(node, SUB, self.expr() )
        elif current_token.kind in (POSITION, PROPERTY):
            node = self.property()
        elif current_token.kind in (MIN, MAX, DIM, ABS, SQRT):
            self.eat(current_token.kind)
            node = self.function(current_token.kind)
        elif current_token.kind == LPAREN:
            self.eat(LPAREN)
            node = self.expr()  # # 创建运算符节点对象
            self.eat(RPAREN)
        while self.current_token.kind in (POW, MOD, FLOOR):
            token = self.current_token
            if token.kind == POW:
                self.eat(POW)
            elif token.kind == MOD:
                self.eat(MOD)
            elif token.kind == FLOOR:
                self.eat(FLOOR)
            node = new Node(node,token.kind,self.factor() )
        return node  # 返回运算符节点对象
    
    
    cdef AST* term(self)except NULL:
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


    cdef AST* expr(self)except NULL:
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


    cdef AST* parse(self, str code)except NULL:
        cdef AST*   node
        if not self.isdatabase:
            raise RuntimeError('Parser is not linked to database!')
        self.lexer.tokenize(code)
        self.current_token = self.lexer.get_current_token()
        node = self.expr()
        if self.current_token.kind != END:
            raise ValueError(f'Error in parse function: expression should ended with END rather than {enum_token[self.current_token.kind]}!')
        return node

    
    
    
