from functools import cache
from pm4py.objects.petri_net.obj import PetriNet, Marking
from pm4py.objects.log.util.log import project_traces
from pm4py.objects.petri_net.utils import petri_utils

from confArt.util import VariableMgr, FormulaMgr, proxyWCNF as WCNF, proxyEncoder as Encoder

@cache
def fib(n) -> int:
    if n < 3: # n==1 or n==2
        return 1
    return fib(n-1) + fib(n-2)

class confArtefact:
    EDIT_DISTANCE = 'edit'
    __DISTANCES = [EDIT_DISTANCE,]

    EXAC_ALIGNMENT = 'exact'
    MULT_ALIGNMENT = 'multi'
    ANTI_ALIGNMENT = 'anti'
    __ALIGNMENTS = [EXAC_ALIGNMENT, MULT_ALIGNMENT, ANTI_ALIGNMENT,]

    STD_COST = 'standard'
    MSY_COST = 'max_synchro'
    LIN_COST = 'decay_lin'
    EXP_COST = 'decay_exp'
    FIB_COST = 'decay_fib'
    __COSTS = [STD_COST, MSY_COST, LIN_COST, EXP_COST, FIB_COST,]

    def __init__(self, filler_transition='OMEGA', final_place='CONF_ART_FINAL'):
        self._distance = confArtefact.EDIT_DISTANCE
        self._isfullRun = False
        self._runsize = 0
        self._prefix = 0
        self._num_traces = 0
        self._optim_sum = False
        self._filler_label = filler_transition
        self._finalplace_name = final_place
        self._costfunction = confArtefact.cost_standard

    # region Properties

    @property
    def distanceMetric(self):
        return self._distance
    @distanceMetric.setter
    def distanceMetric(self, value):
        assert value in confArtefact.__DISTANCES
        self._distance = value

    @property
    def fullRun(self):
        return self._isfullRun
    @fullRun.setter
    def fullRun(self, value):
        self._isfullRun = bool(value)

    @property
    def runsize(self):
        return self._runsize
    @runsize.setter
    def runsize(self, value):
        self._runsize = int(value)

    @property
    def prefix(self):
        return self._prefix
    @prefix.setter
    def prefix(self, value):
        self._prefix = int(value)

    @property
    def numTraces(self):
        return self._num_traces
    @numTraces.setter
    def numTraces(self, value):
        self._num_traces = int(value)

    @property
    def costFunction(self):
        return self._costfunction
    @costFunction.setter
    def costFunction(self, value):
        self._costfunction = value

    # endregion

    ##########################################################################################

    def output_alignment_result(self, params, solution, cost, outfile=None, **kwargs):
        run, traces, traces_dist = self.check(solution, **kwargs)
        if outfile:
            modelname = params['model_path']
            logname = params['log_trace']
            ntrace = params['num_traces']
            func = params['cost_func']

            # REMOVE ENDING OMEGAS
            first_filler = run.find(self._filler_label)
            run = run[:(first_filler+len(self._filler_label))]

            run = '[' + run + ']'
            trace = '[' + traces[0] + ']'
            with open(outfile, 'w') as ff:
                print(f'c M {modelname}', file=ff, end='\n')  # c L modelname
                print(f'c L {logname}', file=ff, end='\n')    # c L logname
                print(f'c NT {ntrace}', file=ff, end='\n')    # c NT val
                print(f'c F {func}', file=ff, end='\n')       # c F func
                print(f'c C {cost}', file=ff, end='\n')       # c C val
                print(f'c T {trace}', file=ff, end='\n')      # c T [sym,sym,sym]
                print(f'o {run}', file=ff, end='\n')          # o [sym,sym,sym]

    def check(self, solution, verbose=False, **kwargs):
        def getActVars(solution, varname):
            varrange = self.VM.getVarRange(varname)
            return list(filter(lambda x: x>0, solution[varrange[0]-1:varrange[1]]))
        def getActualVals(vars):
            return list(map(lambda x: self.VM.inv(x)[1], vars))

        lst_transitions = self.lst_named_transitions + self.lst_silent_transitions
        # MODEL RUN FOUND
        model_vars = getActualVals(getActVars(solution, 'nu_it'))
        run = ','.join(map(lambda x: lst_transitions[x[1]].label if lst_transitions[x[1]].label else '>>', model_vars))

        # FOR EACH TRACE
        traces = []
        trace_dists = []
        for i in range(self.numEncTraces):
            lid = f'_{i}' if (self.encoding != confArtefact.EXAC_ALIGNMENT) else ''
            log_vars = getActualVals(getActVars(solution, f'lambda_ja{lid}'))
            traces.append(','.join(map(lambda x: self.lst_named_transitions[x[1]].label, log_vars)))
            alpha_vars = getActualVals(getActVars(solution, f'alpha_i{lid}'))
            beta_vars = getActualVals(getActVars(solution, f'beta_j{lid}'))
            ## This will only be correct for (multi-)alignments
            ## Anti-alignments would require a proper edit distance computation between run and trace
            trace_dists.append(len(alpha_vars) + len(beta_vars))

        if verbose:
            print('PLACES', self.lst_places)
            print('TRANSITIONS', list(map(lambda x: x.label if x.label else f'TAU|{x.name}', lst_transitions)))
            print('*'*50)
            print('RUN:', run)
            marking_vars = getActualVals(getActVars(solution, 'm_ip'))
            print('MARKINGS:',
                list(map(lambda z:
                            list(map(lambda x: str(self.lst_places[x[1]]),
                                    filter(lambda y: y[0]==z, marking_vars))),
                        range(self.runsize+1)
                    ))
            )
            maxdst = -1
            mindst = 100000000
            avgdst = 0.0
            for i in range(self.numEncTraces):
                print('-'*50)
                print(f'LOG_{i}:', traces[i])
                print(f'EDIT DST:', trace_dists[i])
                maxdst = trace_dists[i] if trace_dists[i] > maxdst else maxdst
                mindst = trace_dists[i] if trace_dists[i] < mindst else mindst
                avgdst = (i*avgdst + trace_dists[i]) / (i+1)
            print('*'*50)
            print('MAX EDIT DST', maxdst)
            print('MIN EDIT DST', mindst)
            print('AVG EDIT DST', avgdst)

        return run, traces, trace_dists

    ##########################################################################################

    def __add_finalmarking_from_any_place(self, model, **kwargs):
        self._finalplace = PetriNet.Place(f'{self._finalplace_name}_P')
        for transition in model.transitions:
            petri_utils.add_arc_from_to(transition, self._finalplace, model)
        model.places.add(self._finalplace)

    def __addfillertransition(self, model, mf, **kwargs):
        ## add filler transition to model
        self._filler_transition = PetriNet.Transition(self._filler_label, label=self._filler_label)
        for place, _ in mf.items():
            petri_utils.add_arc_from_to(place, self._filler_transition, model)
            petri_utils.add_arc_from_to(self._filler_transition, place, model)
        model.transitions.add(self._filler_transition)

    def __encodePetriNet(self, model:PetriNet, m0:Marking, mf:Marking, **kwargs):
        lst_places = sorted(list(model.places), key=lambda x:str(x))

        lst_named_transitions = sorted(
            list(filter(lambda y: y.label is not None, model.transitions)),
            key=lambda x:str(x)
        )

        lst_silent_transitions = sorted(
            list(filter(lambda y: y.label is None, model.transitions)),
            key=lambda x:str(x)
        )

        lst_transitions = lst_named_transitions + lst_silent_transitions

        self._filler_index = lst_transitions.index(self._filler_transition)

        ## create_vars
        ### m_ip (i, p) = (num instants, num places)
        self.VM.addVar('m_ip', [(0, self.runsize+1), (0, len(lst_places))])
        ### nu_it (i, t) = (num instants, num transitions)
        self.VM.addVar('nu_it', [(1, self.runsize), (0, len(lst_transitions))])

        self.formula.extend_vars(self.VM.lastVar - self.formula.max_var())

        ### initial marking
        for id, place in enumerate(lst_places):
            if place in m0:
                self.formula.add_clause([self.VM('m_ip', 0, id)])
            else:
                self.formula.add_clause([-self.VM('m_ip', 0, id)])
        ### final marking
        if self.fullRun:
            for id, place in enumerate(lst_places):
                if place in mf:
                    self.formula.add_clause([self.VM('m_ip', self.runsize, id)])
                else:
                    self.formula.add_clause([-self.VM('m_ip', self.runsize, id)])
        # Per instant i
        for i in range(1, self.runsize+1):
            current_max_var = self.formula.max_var()
            ### EO(t)
            _, clauses = Encoder.exactly_one([self.VM('nu_it', i, tid) for tid, _ in enumerate(lst_transitions)], max_var=current_max_var)
            self.formula.add_clauses(clauses)

            for tid, t in enumerate(lst_transitions):
                tV = self.VM('nu_it', i, tid)
                t_preset = set(map(lambda x: x.source, t.in_arcs))
                t_postset = set(map(lambda x: x.target, t.out_arcs))

                for pid, p in enumerate(lst_places):
                    if p in t_preset:
                        ### nu_it -> m_im1p
                        self.formula.add_clause([-tV, self.VM('m_ip', i-1, pid)])
                        ### nu_it -> -mi_p
                        if p not in t_postset:
                            self.formula.add_clause([-tV, -self.VM('m_ip', i, pid)])

                    if p in t_postset:
                        ### nu_it -> m_ip
                        self.formula.add_clause([-tV, self.VM('m_ip', i, pid)])

                    ### nu_it -> (m_i <-> m_im1p)
                    if p not in t_preset and p not in t_postset:
                        m_ip = self.VM('m_ip', i, pid)
                        m_im1p = self.VM('m_ip', i-1, pid)
                        self.formula.add_clauses([
                            [-tV, -m_ip, m_im1p],
                            [-tV, m_ip, -m_im1p],
                        ])

                if i < self.runsize:
                    ### nu_i,omega -> nu_i+1,omega
                    if tid == self._filler_index:
                        self.formula.add_clause([
                            -self.VM('nu_it', i, tid), self.VM('nu_it', i+1, tid)
                        ])

        return lst_places, lst_named_transitions, lst_silent_transitions

    def __encodeLogTrace(self, lst_named_transitions, trace, traceIdx='', **kwargs):
        ## modify log (add filler character)
        trace.append(self._filler_label)
        ## create lambda_ja vars
        last_var_used = self.VM.addVar(f'lambda_ja{traceIdx}', [(1, len(trace)), (0, len(lst_named_transitions))], self.formula.max_var())
        self.formula.extend_vars(last_var_used - self.formula.max_var())
        ## activate only the lambda_ja such that "a" corresponds to the character in pos "j"
        for j in range(len(trace)):
            for tid, t in enumerate(lst_named_transitions):
                if trace[j] == t.label: ## str(t) returns transition label (if it has, else uses name)
                    self.formula.add_clause([self.VM(f'lambda_ja{traceIdx}', j+1, tid)])
                else:
                    self.formula.add_clause([-self.VM(f'lambda_ja{traceIdx}', j+1, tid)])

    def __encodeMatching(self, lst_named_transitions, lst_silent_transitions, trace_len, tID='', **kwargs):
        ## create tau_ji vars
        self.VM.addVar(f'tau_ji{tID}', [(1, trace_len), (1, self.runsize)])
        ## create alpha_i vars
        self.VM.addVar(f'alpha_i{tID}', [(1, self.runsize)])
        ## create beta_j vars
        self.VM.addVar(f'beta_j{tID}', [(1, trace_len)])

        self.formula.extend_vars(self.VM.lastVar - self.formula.max_var())

        # region Route through tau matrix
        ## tau_1,1 and tau_L,N
        self.formula.add_clauses([
            [self.VM(f'tau_ji{tID}', 1, 1)],
            [self.VM(f'tau_ji{tID}', trace_len, self.runsize)]
        ])

        for j in range(1, trace_len+1):
            for i in range(1, self.runsize+1):
                ## tau_ji -> tau_j,i+1 or tau_j+1,i or tau_j+1,i+1
                next_tau = []
                if j<trace_len:                     next_tau.append(self.VM(f'tau_ji{tID}', j+1, i))
                if i<self.runsize:                  next_tau.append(self.VM(f'tau_ji{tID}', j, i+1))
                if j<trace_len and i<self.runsize:  next_tau.append(self.VM(f'tau_ji{tID}', j+1, i+1))
                if next_tau:
                    self.formula.add_clause([-self.VM(f'tau_ji{tID}', j, i), *next_tau])
                ## tau_ji -> tau_j,i-1 or tau_j-1,i or tau_j-1,i-1
                prev_tau = []
                if j>1:         prev_tau.append(self.VM(f'tau_ji{tID}', j-1, i))
                if i>1:         prev_tau.append(self.VM(f'tau_ji{tID}', j, i-1))
                if j>1 and i>1: prev_tau.append(self.VM(f'tau_ji{tID}', j-1, i-1))
                if prev_tau:
                    self.formula.add_clause([-self.VM(f'tau_ji{tID}', j, i), *prev_tau])

                ## (tau_ji and tau_j+1,i) -> (not tau_j,i+1 and not tau_j-1,i+1 and not tau_j+1,i-1)
                if j<trace_len:
                    vert_taus = [-self.VM(f'tau_ji{tID}', j, i), -self.VM(f'tau_ji{tID}', j+1, i)]
                    if i<self.runsize:          self.formula.add_clause([*vert_taus, -self.VM(f'tau_ji{tID}', j, i+1)])
                    if i>1:                     self.formula.add_clause([*vert_taus, -self.VM(f'tau_ji{tID}', j+1, i-1)])
                    if i<self.runsize and j>1:  self.formula.add_clause([*vert_taus, -self.VM(f'tau_ji{tID}', j-1, i+1)])
                ## (tau_ji and tau_j,i+1) -> (not tau_j+1,i and not tau_j-1,i+1 and not tau_j+1,i-1)
                if i<self.runsize:
                    hori_taus = [-self.VM(f'tau_ji{tID}', j, i), -self.VM(f'tau_ji{tID}', j, i+1)]
                    if j<trace_len:         self.formula.add_clause([*hori_taus, -self.VM(f'tau_ji{tID}', j+1, i)])
                    if j>1:                 self.formula.add_clause([*hori_taus, -self.VM(f'tau_ji{tID}', j-1, i+1)])
                    if j<trace_len and i>1: self.formula.add_clause([*hori_taus, -self.VM(f'tau_ji{tID}', j+1, i-1)])
                ## (tau_ji and tau_j+1,i+1) -> (not tau_j-1,i+1 and not tau_j+1,i-1)
                if j<trace_len and i<self.runsize:
                    diag_taus = [-self.VM(f'tau_ji{tID}', j, i), -self.VM(f'tau_ji{tID}', j+1, i+1)]
                    if j>1: self.formula.add_clause([*diag_taus, -self.VM(f'tau_ji{tID}', j-1, i+1)])
                    if i>1: self.formula.add_clause([*diag_taus, -self.VM(f'tau_ji{tID}', j+1, i-1)])
                    if j<trace_len-1:       self.formula.add_clause([*diag_taus, -self.VM(f'tau_ji{tID}', j+2, i)])
                    if i<self.runsize-1:    self.formula.add_clause([*diag_taus, -self.VM(f'tau_ji{tID}', j, i+2)])

        if not kwargs.get('disable_amo', False):
            ## AMO(up-right diagonal)
            for i in range(1 - trace_len + 1, self.runsize+1):
                current_max_var = self.formula.max_var()
                k = i
                this_diagonal = []
                for j in range(trace_len, 0, -1):
                    if 1 <= k <= self.runsize:
                        this_diagonal.append( self.VM(f'tau_ji{tID}', j, k) )
                    k += 1
                _, clauses = Encoder.at_most_one(this_diagonal, max_var=current_max_var)
                self.formula.add_clauses(clauses)
        # endregion

        # region alpha and beta implications
        f = FormulaMgr(maxvar=self.formula.max_var())

        nu_eq_lambda = {}
        for i in range(1, self.runsize+1):
            or_tau_eq_lst = []
            ### nu_it is silent
            for st_id, _ in enumerate(lst_silent_transitions):
                or_tau_eq_lst.append(
                    f.LIT(self.VM('nu_it', i, len(lst_named_transitions)+st_id))
                )

            for j in range(1, trace_len+1):
                Ni_eq_Lj_lst = []
                #### OR_t nu_it and lambda_jt
                for tid, _ in enumerate(lst_named_transitions):
                    Ni_eq_Lj_lst.append(
                        f.AND([
                            f.LIT(self.VM('nu_it', i, tid)),
                            f.LIT(self.VM(f'lambda_ja{tID}', j, tid))
                        ])
                    )
                nu_eq_lambda[(j, i)], or_clauses = f.OR(Ni_eq_Lj_lst)()
                self.formula.add_clauses(or_clauses)

                ## (tau_ji and [Nu_i = Lambda_j]) -> (tau_j+1,i+1 and not tau_j,i+1 and not tau_j+1,i)
                if i < self.runsize:
                    next_j = j if j == trace_len else (j + 1)
                    self.formula.add_clause([
                        -self.VM(f'tau_ji{tID}', j, i),
                        -nu_eq_lambda[(j, i)],
                        self.VM(f'tau_ji{tID}', next_j, i+1)
                    ])
                    if j < trace_len:
                        self.formula.add_clauses([
                            [-self.VM(f'tau_ji{tID}', j, i), -nu_eq_lambda[(j, i)], -self.VM(f'tau_ji{tID}', j, i+1)],
                            [-self.VM(f'tau_ji{tID}', j, i), -nu_eq_lambda[(j, i)], -self.VM(f'tau_ji{tID}', j+1, i)]
                        ])

                tauji_nu_neq_lambda_impl = []
                if i < self.runsize:    tauji_nu_neq_lambda_impl.append(self.VM(f'tau_ji{tID}', j, i+1))
                if j < trace_len:       tauji_nu_neq_lambda_impl.append(self.VM(f'tau_ji{tID}', j+1, i))
                ## (tau_ji and [Nu_i != Lambda_j]) -> (tau_j+1,i or tau_j,i+1)
                if tauji_nu_neq_lambda_impl:
                    self.formula.add_clause([-self.VM(f'tau_ji{tID}', j, i), nu_eq_lambda[(j, i)], *tauji_nu_neq_lambda_impl])

                ## LST tau_ji and [Nu_i = Lambda_j]
                or_tau_eq_lst.append(
                    f.AND([
                        f.LIT(self.VM(f'tau_ji{tID}', j, i)),
                        f.LIT(nu_eq_lambda[(j, i)])
                    ])
                )

            ## (OR_j tau_ji and [Nu_i = Lambda_j]) <-> alpha_i
            or_var, clauses = f.OR(or_tau_eq_lst)()
            self.formula.add_clauses(clauses)
            self.formula.add_clauses([
                [-or_var, -self.VM(f'alpha_i{tID}', i)],
                [ or_var,  self.VM(f'alpha_i{tID}', i)]
            ])

        for j in range(1, trace_len):
            tauji_tauj1i = []
            for i in range(1, self.runsize+1):
                tauji_tauj1i.append(
                    f.AND([
                        f.LIT(self.VM(f'tau_ji{tID}', j, i)),
                        f.LIT(self.VM(f'tau_ji{tID}', j+1, i))
                    ])
                )

            ## (OR_i tau_ji and tau_j+1,i) <-> beta_j
            or_var, clauses = f.OR(tauji_tauj1i)()
            self.formula.add_clauses(clauses)
            self.formula.add_clauses([
                [-or_var,  self.VM(f'beta_j{tID}', j)],
                [ or_var, -self.VM(f'beta_j{tID}', j)]
            ])

        # endregion

    def __encodeFullLog(self, lst_named_transitions, lst_silent_transitions, projected_traces, **kwargs):
        self.numTraces = min(len(projected_traces), self.numTraces)
        for traceid, trace in enumerate(projected_traces[:self.numTraces]):
            if self.prefix != 0:
                retrace = trace[:self.prefix]
            else:
                retrace = trace
            self.__encodeLogTrace(lst_named_transitions, retrace, traceIdx=f'_{traceid}', **kwargs)
            self.__encodeMatching(lst_named_transitions, lst_silent_transitions, len(retrace), tID=f'_{traceid}', **kwargs)

    def alignment(self, model, m0, mf, traces, **kwargs):
        assert(self.numTraces <= len(traces))

        self.formula = WCNF()
        self.VM = VariableMgr()

        if not self.fullRun:
            ## if unfinished log, add final place reachable from any transition in the model
            ## changing the final marking, then addfiller will create a loop on that final place
            self.__add_finalmarking_from_any_place(model, **kwargs)
            mf = Marking()
            mf[self._finalplace] = 1

        ## modify model (add filler transition)
        self.__addfillertransition(model, mf, **kwargs)
        ## encode petrinet
        self.lst_places, self.lst_named_transitions, self.lst_silent_transitions = self.__encodePetriNet(model, m0, mf, **kwargs)
        ## encode log
        projected_trace = project_traces(traces)[self.numTraces-1]
        if self.prefix != 0:
            retrace = projected_trace[:self.prefix]
        else:
            retrace = projected_trace
        #self.__encodeLogTrace(self.lst_named_transitions, projected_trace, **kwargs)
        self.__encodeLogTrace(self.lst_named_transitions, retrace, **kwargs)
        ## encode matching
        #self.__encodeMatching(self.lst_named_transitions, self.lst_silent_transitions, len(projected_trace), **kwargs)
        self.__encodeMatching(self.lst_named_transitions, self.lst_silent_transitions, len(retrace), **kwargs)
        ##encode soft clauses
        self.costFunction(self,
            'alpha_i', self.runsize, -1,
            'beta_j', len(retrace), -1
        )

        self.numEncTraces = 1
        self.encoding = confArtefact.EXAC_ALIGNMENT
        return self.formula

    def multi_alignment(self, model, m0, mf, traces, **kwargs):
        self.formula = WCNF()
        self.VM = VariableMgr()

        if not self.fullRun:
            ## if unfinished log, add final place reachable from any transition in the model
            ## changing the final marking, then addfiller will create a loop on that final place
            self.__add_finalmarking_from_any_place(model, **kwargs)
            mf = Marking()
            mf[self._finalplace] = 1

        ## modify model (add filler transition)
        self.__addfillertransition(model, mf, **kwargs)
        ## encode petrinet
        self.lst_places, self.lst_named_transitions, self.lst_silent_transitions = self.__encodePetriNet(model, m0, mf, **kwargs)

        ## encode full log and matching
        projected_traces = project_traces(traces)
        self.__encodeFullLog(self.lst_named_transitions, self.lst_silent_transitions, projected_traces, **kwargs)
        #encode soft clauses
        for traceID in range(self.numTraces):
            lenTrace = min(self.prefix, len(projected_traces[traceID])) if self.prefix != 0 else len(projected_traces[traceID])
            self.costFunction(self,
                f'alpha_i_{traceID}', self.runsize, -1,
                f'beta_j_{traceID}', lenTrace, -1
            )

        self.numEncTraces = self.numTraces
        self.encoding = confArtefact.MULT_ALIGNMENT
        return self.formula

    def anti_alignment(self, model, m0, mf, traces, **kwargs):
        self.formula = WCNF()
        self.VM = VariableMgr()

        if not self.fullRun:
            ## if unfinished log, add final place reachable from any transition in the model
            ## changing the final marking, then addfiller will create a loop on that final place
            self.__add_finalmarking_from_any_place(model, **kwargs)
            mf = Marking()
            mf[self._finalplace] = 1

        ## modify model (add filler transition)
        self.__addfillertransition(model, mf, **kwargs)
        ## encode petrinet
        self.lst_places, self.lst_named_transitions, self.lst_silent_transitions = self.__encodePetriNet(model, m0, mf, **kwargs)
        ## encode full log and matching
        projected_traces = project_traces(traces)
        self.__encodeFullLog(self.lst_named_transitions, self.lst_silent_transitions, projected_traces, **kwargs)
        #encode soft clauses
        for traceID in range(self.numTraces):
            lenTrace = min(self.prefix, len(projected_traces[traceID])) if self.prefix != 0 else len(projected_traces[traceID])
            self.costFunction(self,
                f'alpha_i_{traceID}', self.runsize, 1,
                f'beta_j_{traceID}', lenTrace, 1
            )

        self.numEncTraces = self.numTraces
        self.encoding = confArtefact.ANTI_ALIGNMENT
        return self.formula

    #region COST FUNCTIONS

    def cost_standard(self, ivar_ID, len_i, i_sign, jvar_ID, len_j, j_sign):
        for i in range(1, len_i + 1):
            self.formula.add_clause([i_sign * self.VM(ivar_ID, i)], weight=1)
        for j in range(1, len_j + 1):
            self.formula.add_clause([j_sign * self.VM(jvar_ID, j)], weight=1)

    def cost_max_synchro(self, ivar_ID, len_i, i_sign, jvar_ID, len_j, j_sign):
        for i in range(1, len_i + 1):
            self.formula.add_clause([i_sign * self.VM(ivar_ID, i)], weight=1)
        for j in range(1, len_j + 1):
            self.formula.add_clause([j_sign * self.VM(jvar_ID, j)], weight=len_i)

    def cost_decay_lin(self, ivar_ID, len_i, i_sign, jvar_ID, len_j, j_sign):
        lin_fkM = lambda k, M: M - k + 1

        for i in range(1, len_i + 1):
            self.formula.add_clause([i_sign * self.VM(ivar_ID, i)], weight=lin_fkM(i, len_i))
        for j in range(1, len_j + 1):
            self.formula.add_clause([j_sign * self.VM(jvar_ID, j)], weight=lin_fkM(j, len_j))

    def cost_decay_exp(self, ivar_ID, len_i, i_sign, jvar_ID, len_j, j_sign):
        exp_fkM = lambda k, M: 2**(M - k)

        for i in range(1, len_i + 1):
            self.formula.add_clause([i_sign * self.VM(ivar_ID, i)], weight=exp_fkM(i, len_i))
        for j in range(1, len_j + 1):
            self.formula.add_clause([j_sign * self.VM(jvar_ID, j)], weight=exp_fkM(j, len_j))

    def cost_decay_fib(self, ivar_ID, len_i, i_sign, jvar_ID, len_j, j_sign):
        fib_fkm = lambda k, M: fib(M - k + 1)

        for i in range(1, len_i + 1):
            self.formula.add_clause([i_sign * self.VM(ivar_ID, i)], weight=fib_fkm(i, len_i))
        for j in range(1, len_j + 1):
            self.formula.add_clause([j_sign * self.VM(jvar_ID, j)], weight=fib_fkm(j, len_j))

    #endregion
