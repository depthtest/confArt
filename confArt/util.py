from abc import ABC, abstractmethod

class VariableMgr:
    ###################################################################
    def GetVarMapper(startingvar, indexes):
        '''
        Returns a mapper from indexes to DIMACS var
        :param startingvar (int): first available DIMACS var
        :param indexes list(tuples(int,int)): list of required indices, each index represented as a tuple (starting id, number of ids)
        '''
        class VarMapper:
            def __init__(self, st_var, idxs):
                self._st_var = st_var
                self._idxs = idxs
                self._params = [1]
                for id in reversed(idxs):
                    new_par = self._params[0]*(id[1])
                    self._params.insert(0, new_par)
                self._params.pop(0)
                max_combo = tuple(map(lambda x:x[0]+x[1]-1, self._idxs))
                self._maxvar = self(*max_combo)

            def __call__(self, *args):
                assert len(args) == len(self._params), f'This variable needs {len(self._params)}. It was passed {len(args)} arguments.'
                for id, p in enumerate(args):
                    assert self._idxs[id][0] <= p < (self._idxs[id][0] + self._idxs[id][1]), f'Index {id}={p} is out of range [{self._idxs[id][0]},{self._idxs[id][0] + self._idxs[id][1]-1}]'
                return self._st_var + sum(
                    map(lambda x: ((x[0]-x[2][0]))*x[1], 
                        zip(args, self._params, self._idxs)
                    )
                )

            def inv(self, intVar):
                assert self._st_var <= intVar <= self._maxvar, f'Var {intVar} is out of range [{self._st_var}, {self._maxvar}]'
                var_idxs = []
                acc = intVar - self._st_var
                for id, p in enumerate(self._params):
                    var_idxs.append((acc // p) + self._idxs[id][0])
                    acc = acc % p
                return tuple(var_idxs)

        return VarMapper(startingvar, indexes)

    ###################################################################

    def __init__(self):
        self._vars = {}
        self._last_var = 0

    @property
    def lastVar(self):
        return self._last_var

    def addVar(self, label, index_ranges, last_used_var = -1):
        assert label not in self._vars, f'Label {label} already exists in VarManager'
        if last_used_var == -1:
            last_used_var = self._last_var
        self._vars[label] = VariableMgr.GetVarMapper(last_used_var + 1, index_ranges)
        self._last_var = self._vars[label]._maxvar
        return self._last_var

    def getVarIdxRanges(self, label):
        varMapper = self._vars[label]
        return varMapper._idxs

    def getVarRange(self, label):
        varMapper = self._vars[label]
        return (
            varMapper._st_var,  ## Min Var
            varMapper._maxvar  ## Max Var
        )

    def __call__(self, label, *args):
        varMapper = self._vars[label]
        assert len(args) == len(varMapper._params), f'Variable {label} needs {len(varMapper._params)} parameters. It was passed {len(args)} arguments.'
        return varMapper(*args)

    def inv(self, intVar):
        for varL, varM in self._vars.items():
            minV, maxV = self.getVarRange(varL)
            if minV <= intVar <= maxV:
                return (varL, varM.inv(intVar))
        return None


class FormulaMgr:
    def __init__(self, maxvar = 0):
        self._maxVar = maxvar
        self._cache = {}
    
    def _genVarFunc(self):
        self._maxVar += 1
        return self._maxVar

    def LIT(self, lit):
        return FormulaMgr._LIT(self, lit)
    def OR(self, arg_lst):
        return FormulaMgr._OR(self, arg_lst)
    def AND(self, arg_lst):
        return FormulaMgr._AND(self, arg_lst)
    def IFT(self, left, right):
        return FormulaMgr._IFT(self, left, right)
    def IFF(self, left, right):
        return FormulaMgr._IFF(self, left, right)
    def NOT(self, arg):
        return FormulaMgr._NOT(self, arg)

    class LogOp(ABC):
        def __init__(self, mgr, opType, *args, **kwargs):
            self._mgr = mgr
            self._type = opType

        def _get_lit(self):
            this_args = self._get_arg_rep()
            # with the hash requires less memory, haven't check if there are collisions though...
            #this_args = hash(self._get_arg_rep())

            self.lit = self._mgr._cache.get(this_args)
            if self.lit is None:
                self.lit = self._mgr._genVarFunc()
                self._mgr._cache[this_args] = self.lit

        def _get_arg_rep(self):
            if self._type == 'lit':
                return self.lit
            if self._type == 'not' or self._type == 'ift' or self._type == 'iff':
                return (self._type, self.expr._get_arg_rep())
            # and & or types
            arg_lst_rep = []
            for param in self.arg_lst:
                arg_lst_rep.append(param._get_arg_rep())
            return (self._type, frozenset(arg_lst_rep))

        @abstractmethod
        def __call__(self):
            pass

    class _LIT(LogOp):
        def __init__(self, mgr, lit, *args, **kwargs):
            assert isinstance(lit, int), 'lit should be an integer (DIMACS variable index)'
            super().__init__(mgr, 'lit', *args, **kwargs)
            self.lit = lit
        def __call__(self):
            return self.lit, []
        def __str__(self):
            return str(self.lit)

    class _OR(LogOp):
        def __init__(self, mgr, arg_lst, *args, **kwargs):
            for elem in arg_lst:
                assert isinstance(elem, FormulaMgr.LogOp), 'arg_lst should be a list of (subclass) FormulaMgr.LogOp - use the LIT, OR, AND, NOT, IFT and IFF functions'
            super().__init__(mgr, 'or', *args, **kwargs)
            self.arg_lst = frozenset(arg_lst)
            self._get_lit()
        def __call__(self):
            ## generate clauses from arg_lst
            clauses = []
            lits = []
            for op in self.arg_lst:
                op_lit, op_clauses = op()
                lits.append(op_lit)
                clauses.extend(op_clauses)
            # OR() <-> self.lit
            clauses.append([-self.lit, *lits])
            clauses.extend([[-lit, self.lit] for lit in lits])
            ## return ensuring no repeated clauses
            return self.lit, list(map(lambda y: list(y), set(map(lambda x: frozenset(x), clauses)))) 
        def __str__(self):
            return f'OR(' + ','.join(map(lambda x: str(x), self.arg_lst)) + ')'

    class _AND(LogOp):
        def __init__(self, mgr, arg_lst, *args, **kwargs):
            super().__init__(mgr, 'and', *args, **kwargs)
            for elem in arg_lst:
                assert isinstance(elem, FormulaMgr.LogOp), 'arg_lst should be a list of (subclass) FormulaMgr.LogOp - use the LIT, OR, AND, NOT, IFT and IFF functions'
            self.arg_lst = frozenset(arg_lst)
            self._get_lit()
        def __call__(self):
            ## generate clauses from arg_lst
            clauses = []
            lits = []
            for op in self.arg_lst:
                op_lit, op_clauses = op()
                lits.append(op_lit)
                clauses.extend(op_clauses)
            # AND() <-> self.lit
            clauses.append([self.lit, *list(map(lambda x: -x, lits))])
            clauses.extend([[-self.lit, lit] for lit in lits])
            ## return ensuring no repeated clauses
            return self.lit, list(map(lambda y: list(y), set(map(lambda x: frozenset(x), clauses))))
        def __str__(self):
            return 'AND(' + ','.join(map(lambda x: str(x), self.arg_lst)) + ')'

    class _NOT(LogOp):
        def __init__(self, mgr, arg, *args, **kwargs):
            assert isinstance(arg, FormulaMgr.LogOp), 'arg should be of type (subclass) FormulaMgr.LogOp - use the LIT, OR, AND, NOT, IFT and IFF functions'
            super().__init__(mgr, 'not', *args, **kwargs)
            if type(arg) == FormulaMgr._LIT:
                self.expr =  FormulaMgr._LIT(self._mgr, -arg.lit, *args, **kwargs)
            elif type(arg) == FormulaMgr._OR:
                self.expr =  FormulaMgr._AND(self._mgr, [FormulaMgr._NOT(self._mgr, or_arg) for or_arg in arg.arg_lst], *args, **kwargs)
            elif type(arg) == FormulaMgr._AND:
                self.expr =  FormulaMgr._OR(self._mgr, [FormulaMgr._NOT(self._mgr, or_arg) for or_arg in arg.arg_lst], *args, **kwargs)
            elif type(arg) == FormulaMgr._NOT:
                self.expr =  arg
            else: # FormulaMgr._IFT or FormulaMgr._IFF
                self.expr =  FormulaMgr._NOT(arg.expr, *args, **kwargs)
            self.lit = self.expr.lit
        def __call__(self):
            return self.expr()
        def __str__(self):
            return f'NOT({self.expr})'

    class _IFT(LogOp):
        def __init__(self, mgr, arg_left, arg_right, *args, **kwargs):
            assert isinstance(arg_left, FormulaMgr.LogOp), 'arg_left should be of type (subclass) FormulaMgr.LogOp - use the LIT, OR, AND, NOT, IFT and IFF functions'
            assert isinstance(arg_right, FormulaMgr.LogOp), 'arg_right should be of type (subclass) FormulaMgr.LogOp - use the LIT, OR, AND, NOT, IFT and IFF functions'
            super().__init__(mgr, 'ift', *args, **kwargs)
            self.expr = FormulaMgr._OR(self._mgr, [
                            FormulaMgr._NOT(self._mgr, arg_left),
                            arg_right
                        ], *args, **kwargs)
            self.lit = self.expr.lit
        def __call__(self):
            return self.expr()
        def __str__(self):
            return str(self.expr)

    class _IFF(LogOp):
        def __init__(self, mgr, arg_left, arg_right, *args, **kwargs):
            assert isinstance(arg_left, FormulaMgr.LogOp), 'arg_left should be of type (subclass) FormulaMgr.LogOp - use the LIT, OR, AND, NOT, IFT and IFF functions'
            assert isinstance(arg_right, FormulaMgr.LogOp), 'arg_right should be of type (subclass) FormulaMgr.LogOp - use the LIT, OR, AND, NOT, IFT and IFF functions'
            super().__init__(mgr, 'iff', *args, **kwargs)
            self.expr = FormulaMgr._AND(self._mgr, [
                            FormulaMgr._IFT(self._mgr, arg_left, arg_right),
                            FormulaMgr._IFT(self._mgr, arg_right, arg_left)
                        ], *args, **kwargs)
            self.lit = self.expr.lit
        def __call__(self):
            return self.expr()
        def __str__(self):
            return str(self.expr)
