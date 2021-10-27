import pathlib
import argparse
import json
import re

from pm4py import read_pnml, read_xes

from confArt.confArt import confArtefact as CA

DST_CHOICES = [CA.EDIT_DISTANCE,]
ART_FUNCS = {
    CA.EXAC_ALIGNMENT : CA.alignment,
    CA.MULT_ALIGNMENT : CA.multi_alignment,
    CA.ANTI_ALIGNMENT : CA.anti_alignment,
}

def main_encode(args, other_args):
    encoder_argp = argparse.ArgumentParser()
    encoder_argp.add_argument('--disable_amo', action='store_true')
    parser_other = encoder_argp.parse_args(other_args)

    model, m0, mf = read_pnml(str(args.model.absolute()))
    traces = read_xes(str(args.log_trace.absolute()))

    ca = CA()
    ca.distanceMetric = args.distance
    ca.fullRun = args.reach_final
    ca.runsize = args.run_size + args.max_d
    ca.numTraces = args.num_traces

    wcnf = ART_FUNCS[args.artefact](ca, model, m0, mf, traces, **vars(parser_other))
    params = {
        'artefact' : args.artefact,
        'model' : str(args.model.absolute()),
        'log_trace' : str(args.log_trace.absolute()),
        'distance' : args.distance,
        'reach_final' : args.reach_final,
        'run_size' : args.run_size,
        'max_d' : args.max_d,
        'num_traces' : args.num_traces,
    }
    params.update(vars(parser_other))
    wcnf.add_comment('Partial MaxSAT Conformance Checking Artefacts')
    wcnf.add_comment(f'PARAMS:{json.dumps(params)}')
    wcnf.write_dimacs_file(args.out_wcnf)

def main_check(args, other_args):
    checker_argp = argparse.ArgumentParser()
    checker_argp.add_argument('--verbose', action='store_true')
    parser_other = checker_argp.parse_args(other_args)

    def parse_maxsat_solution(maxsat_file):
        solution = []
        with open(maxsat_file) as mf:
            for line in mf:
                matched = re.match(r'^v (.*)', line)
                if matched is None:
                    continue
                else:
                    solution = list(map(lambda x: int(x), matched[1].split(' ')))
        return solution

    ## Search PARAMS comment in formula
    with open(args.maxsat_formula) as f:
        for line in f:
            matched = re.match(r'^c PARAMS:(.*)$', line)
            if matched is not None:
                params = json.loads(matched.group(1))
                break

    model, m0, mf = read_pnml(params.pop('model'))
    traces = read_xes(params.pop('log_trace'))

    ca = CA()
    ca.distanceMetric = params.pop('distance')
    ca.fullRun = params.pop('reach_final')
    ca.runsize = params.pop('run_size') + params.pop('max_d')
    ca.numTraces = params.pop('num_traces')

    artefact = params.pop('artefact')

    _ = ART_FUNCS[artefact](ca, model, m0, mf, traces, **params)

    solution = parse_maxsat_solution(args.maxsat_output)
    ca.check(solution, parser_other.verbose)


if __name__ == '__main__':
    argp = argparse.ArgumentParser()
    subparsers = argp.add_subparsers(required=True, dest='command')

    ## ENCODER
    encode_argp = subparsers.add_parser('encode')
    encode_argp.add_argument('model', type=pathlib.Path, help='Path to the Model PetriNet file (PNML)')
    encode_argp.add_argument('log_trace', type=pathlib.Path, help='Path to the Log file (XES)')
    encode_argp.add_argument('artefact', choices=ART_FUNCS.keys(), help='Artefact to compute')
    encode_argp.add_argument('out_wcnf', type=pathlib.Path, help='Path to the resulting WCNF formula')
    encode_argp.add_argument('--num_traces', type=int, required=True, help='Number of traces from log to encode')
    encode_argp.add_argument('--run_size', type=int, required=True, help='Size of the run (n) to encode')
    encode_argp.add_argument('--max_d', type=int, required=True, help='Max edit distance allowed to encode')
    encode_argp.add_argument('--reach_final', action='store_true', help='Encode full model -- force final marking | def: %(default)r')
    encode_argp.add_argument('--distance', choices=DST_CHOICES, default=DST_CHOICES[0])
    encode_argp.set_defaults(func=main_encode)

    ## CHECKER
    check_argp = subparsers.add_parser('check')
    check_argp.add_argument('maxsat_formula', type=pathlib.Path, help='Path to the wcnf file describing the formula')
    check_argp.add_argument('maxsat_output', type=pathlib.Path, help='Path to the output file from the maxsat solver')
    check_argp.set_defaults(func=main_check)


    args, other_args = argp.parse_known_args()
    args.func(args, other_args)

