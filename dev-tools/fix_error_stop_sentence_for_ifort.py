from argparse import ArgumentParser
import re


def parse_args():
    parser = ArgumentParser()

    parser.add_argument('filename')

    args = parser.parse_args()

    return args



def main(args):
    with open(args.filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    with open(args.filename, 'w', encoding='utf-8') as f:
        for line in lines:
            m = re.match(r' *if *\((.+)\) *error stop *(.+)', line)
            if m:
                f.write(f'if ({m.group(1)}) then\n  write(0, *) {m.group(2)}\n  error stop\nend if\n')
            else:
                f.write(line)




if __name__ == '__main__':
    args = parse_args()
    main(args)
