import argparse
import sys
import re

___author__ = 'Jelle Becirspahic'
___data__ = '16-01-2023'
___version__ = 'v1.0'


class VcfFilter:
    def __init__(self, vcf_file, quality, depth, output):
        self.vcf_file = open(vcf_file)
        self.quality = quality
        self.depth = depth
        self.output = output

    def filter_by_quality(self):
        message = ""
        with open(self.output, 'w') as new_file:
            for line in self.vcf_file:
                if line.startswith("#"):
                    new_file.write(line)
                try:
                    if not line.startswith('#'):
                        variant_quality =  float(line.split('\t')[5])
                        if variant_quality >= int(self.quality):
                            new_file.write(line)
                except ValueError:
                    return 0
                except TypeError:
                    message += "fill in the -q"
            if message != "":
                new_file.close()

    def filter_by_depth(self):
        with open(self.output, 'w') as new_file:
            for line in self.vcf_file:
                if line.startswith("#"):
                    new_file.write(line)

                if not line.startswith('#'):
                    match = re.search("(DP=)(\d+)", str(line.split('\t')[7]))
                    if int(match.group(2)) >= int(self.depth):
                        new_file.write(line)
            new_file.close()

    def decomposer_filter(self):
        with open(self.output, 'w') as new_file:
            for line in self.vcf_file:
                if line.startswith("#"):
                    new_file.write(line)

                else:
                    list_container = line.split('\t')
                    if len(list_container[4].split(',')) > 1:
                        for i in list_container[4].split(','):
                            list_container[4] = i
                            new_file.write('\t'.join([str(new_string)
                                                      for new_string in list_container]))
                    else:
                        new_file.write(line)
            new_file.close()

        return 0


def main(args):
    inp_args = argparse.ArgumentParser()

    inp_args.add_argument('method', help='Filter method to be used I.E. decompose, depth or quality')
    inp_args.add_argument('-i', '--infile', help='Name of the input file, with filepath')
    inp_args.add_argument('-q', '--quality', help='The upper limit of the quality of a variant.')
    inp_args.add_argument('-d', '--depth', help='The minimum depth that a variant must have')
    inp_args.add_argument('-o', '--outfile', help='Name of the output file, with filepath')

    args = inp_args.parse_args()
    method = args.method
    vcf_file = args.infile
    quality = args.quality
    depth = args.depth
    output = args.outfile

    run = VcfFilter(vcf_file, quality, depth, output)

    if method == "decompose":
        run.decomposer_filter()
    if method == "quality":
        run.filter_by_quality()
    if method == "depth":
        run.filter_by_depth()

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))