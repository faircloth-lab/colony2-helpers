#!/usr/bin/env python
# encoding: utf-8
"""
File: cervus_to_colony2.py
Author: Brant Faircloth

Created by Brant Faircloth on 17 September 2012 16:09 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import csv
import argparse
from collections import defaultdict

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "genotypes",
            help="""The input CSV file containing genotypes in CERVUS format"""
        )
    parser.add_argument(
            "mothers",
            help="""The CSV file containing putative mother ids"""
        )
    parser.add_argument(
            "fathers",
            help="""The CSV file containing putative father ids """
        )
    parser.add_argument(
            "offspring",
            help="""The CSV file containing putative offspring ids"""
        )
    parser.add_argument(
            "--start-column",
            type=int,
            default=2,
            help="""The (0-indexed) column in which genotype values start""",
        )
    return parser.parse_args()


def get_loci_from_csv_file(input, startcol=2):
    with open(input, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        temp_geno = {}
        for lineno, line in enumerate(reader):
            if lineno == 0:
                column_positions = {pos + startcol:locus for pos, locus in enumerate(line[startcol:])}
                temp_geno = {locus:[] for locus in line[startcol:]}
            else:
                for pos, allele in enumerate(line):
                    if pos in column_positions.keys():
                        temp_geno[column_positions[pos]].append(allele)
        temp_loci = set(['_'.join(name.split('_')[:-1]) for name in column_positions.values()])
        temp_loci = {name:[] for name in temp_loci}
        for fullname, v in temp_geno.iteritems():
            name = '_'.join(fullname.split('_')[:-1])
            temp_loci[name].extend(v)
        loci = {name:{} for name in temp_loci}
        #pdb.set_trace()
        for name, v in temp_loci.iteritems():
            unique = sorted(list(set(v)))
            if '*' in unique:
                unique.remove('*')
                unique.insert(0, '*')
                for pos, u in enumerate(unique):
                    loci[name][u] = pos
            else:
                for pos, u in enumerate(unique):
                    loci[name][u] = pos + 1
    return column_positions, loci


def convert_genotypes_from_csv_file(input, columns, loci, startcol=2):
    individuals = {}
    with open(input, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        for lineno, line in enumerate(reader):
            if lineno == 0:
                pass
            else:
                allele_list = []
                iden = line[0]
                for pos, allele in enumerate(line):
                    if pos in columns.keys():
                        fullname = columns[pos]
                        name = '_'.join(fullname.split('_')[:-1])
                        newval = loci[name][allele]
                        allele_list.append(' {0}'.format(newval))
                individuals[iden] = ''.join(allele_list)
    return individuals


def get_individual_classes(females, males, offspring):
    classes = {}
    for pos, f in enumerate([females, males, offspring]):
        with open(f, 'rU') as input:
            for line in input:
                ls = line.strip().split(',')
                if pos == 0:
                    classes[ls[0]] = 'females'
                elif pos == 1:
                    classes[ls[0]] = 'males'
                elif pos == 2:
                    classes[ls[0]] = 'offspring'
    return classes


def write_output(classes, converted):
    females = open('females-colony2-converted.txt', 'w')
    males = open('males-colony2-converted.txt', 'w')
    offspring = open('offspring-colony2-converted.txt', 'w')
    for iden, genotype in converted.iteritems():
        # get class of individual
        cls = classes[iden]
        if cls == 'females':
            females.write("{0} {1}\n".format(iden, genotype.lstrip(' ')))
        elif cls == 'males':
            males.write("{0} {1}\n".format(iden, genotype.lstrip(' ')))
        elif cls == 'offspring':
            offspring.write("{0} {1}\n".format(iden, genotype.lstrip(' ')))
    for f in [females, males, offspring]:
        f.close()


def main():
    args = get_args()
    classes = get_individual_classes(args.mothers, args.fathers, args.offspring)
    columns, loci = get_loci_from_csv_file(args.csv)
    converted = convert_genotypes_from_csv_file(args.csv, columns, loci)
    write_output(classes, converted)


if __name__ == '__main__':
    main()
