import sys
import argparse
import pickle
import json

def load_file(figurefile):
    fh = open(figurefile, 'rb')
    result = pickle.load(fh)
    fh.close()
    return result

def process_result(result, args):
    to_remove = []
    for scenario in result:
        for scheme in result[scenario]:
            if len(args.schemes) > 0 and scheme not in args.schemes:
                to_remove.append(scheme)
                continue
            if args.params:
                if "params" not in result[scenario][scheme]:
                    to_remove.append(scheme)
                else:
                    result[scenario][scheme] = {
                        "params": result[scenario][scheme]["params"]
                    }
        for r in to_remove:
            result[scenario].pop(r, None) 

def parse_all_args():
    parser = argparse.ArgumentParser(description='Output parameters as JSON.')
    parser.add_argument('--full', help='output full file', action='store_true')
    parser.add_argument('--params', help='only output params', action='store_true')
    parser.add_argument('--pretty', help='output pretty printed params', action='store_true')
    parser.add_argument('figurefile', metavar='figure file', type=str, help='figurefile')
    parser.add_argument('schemes', metavar='schemes', type=str, nargs="*", help='only include these schemes')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_all_args()
    
    result = load_file(args.figurefile)
    if not args.full:
        process_result(result, args)
    indent = 4 if args.pretty else None
    print(json.dumps(result, sort_keys=True, indent=indent))
