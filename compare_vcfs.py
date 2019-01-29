"""
 Script to compare MassArray and WGS intersected VCF results.

 Sarah Burns & Chipo Mashayamombe-Wolfgarten 29 Jan 2019
"""

from ruffus import *
import vcf
import pandas as pd
import re
import glob
import os
from datetime import datetime
import argparse

arg_parser = argparse.ArgumentParser(description='Scripts to compare massarray and wgs vcf results.')
arg_parser.add_argument('-powers', default='/network/processed/100K_VCFs/hg19_total_af_power.txt')
arg_parser.add_argument('-dest', default='/network/processed/100K_VCFs/archive')
args = arg_parser.parse_args()

output_excel = 'IdentityCheck_%s.xlsx' % datetime.now().strftime('%Y%m%d_%H%M')


def get_alleles_from_genotype(ordered_alleles, genotype):
    """ Determine alleles from genotype and list of possible alleles.
    
    E.g. genotype 0/1 and alleles ['A', 'G', 'C'] where A is ref => AG
        
    """
    alleles = ''
    genotype_list = genotype.split('/')
    for idx in genotype_list:
        if idx != '.':
            alleles += ordered_alleles[int(idx)]
    return alleles


def get_risk(fraction):
    """
    calculating the risk: chance of 2 samples having the exact same genotype profile
    """
    my_frac = str(fraction)
    # count the number of zeros after the decimal point
    zeros = len(re.search("\.(0*)", my_frac).group(1))

    # create a string with 1 and (the number of zeros+1)
    new_str = '1' + ('0' * (zeros + 1))

    # convert new_str and my_frac to numeric and multiply them to get the numerator
    numerator = pd.to_numeric(my_frac) * pd.to_numeric(new_str)

    # divide the new string by the numerator to get a 1 in ...
    denominator = int(pd.to_numeric(new_str) / numerator)
    denominator = '{0:,}'.format(denominator)

    return '1 in {}'.format(denominator)


def calculate_power(snp_list):
    # read in the file with precalculated snp-specific probabilities
    tot_power = pd.read_table(args.powers, sep='\t')

    # take a subset of SNPs that have been called in the sample
    snp_subset = tot_power.loc[tot_power['ID'].isin(snp_list)].reset_index(drop=True)

    # calculate the probability of identity and the power of exclusion
    id_prob = pd.to_numeric(snp_subset['Probability of Identity']).prod()

    risk = get_risk(id_prob)

    power = 1 - id_prob
    return power, risk


def archive_files(sample):
    for f in glob.glob('*%s*' % sample):
        # compress
        os.rename(f, os.path.join(args.dest, f))


@collate(['*.vcf'], regex('(.+GRCh3[78]).intersected.vcf'), add_inputs(r'\1.massarray.vcf'), output_excel)
def collate_vcfs(infiles, outfile):
    """ Pair massarray and wgs vcfs with same build and parse all results to excel.
 
    """

    main_cols = ['SAMPLE', 'ARRAY CALLS', 'WGS CALLS', 'ALL MATCHES', 'PROBABILITY OF UNIQUENESS', 'ODDS RATIO',
                 'HIGH QUAL MATCHES', 'HIGH QUAL PROBABILITY OF UNIQUENESS', 'HIGH QUAL ODDS RATIO']

    mismatch_cols = ['SAMPLE', 'SNP', 'WGS GENOTYPE', 'MASSARRAY GENOTYPE', 'QUALITY OF CALL', 'VCF FILTER']

    main_df = pd.DataFrame(columns=main_cols)
    mismatch_df = pd.DataFrame(columns=mismatch_cols)

    all_samples = []

    for (wgs_vcf, array_vcf) in infiles:
        
        # Get lab number
        try:
            sample_name = re.search(r'D\d{2}.\d{5}', wgs_vcf).group(0)
        except AttributeError:
            sample_name = wgs_vcf.split('.')[0].split('_')[0]
        all_samples.append(sample_name)

        array_results = {}
        wgs_results = {}
        coords_to_snp = {}

        # Parse required array results into dict e.g. { 'rs123': { 'alleles': 'AG', 'quality': 'A', 'filter': '.' } }
        array_reader = vcf.Reader(open(array_vcf, 'r'))
        for record in array_reader:
            snp_id = record.ID
            vcf_filter = ','.join(record.FILTER)
            alleles = [str(x) for x in record.ALT]
            alleles.insert(0, str(record.REF))
            coords_to_snp[(record.CHROM, record.POS)] = snp_id
            for sample in record.samples:
                gt = sample['GT']
                quality = sample['MTQ']
                alleles_in_sample = get_alleles_from_genotype(alleles, gt)
                array_results[snp_id] = {
                    'alleles': ''.join(sorted(alleles_in_sample)), 'quality': quality, 'filter': vcf_filter
                }

        # Parse required wgs results into dict e.g. { 'rs123': 'AG' }
        wgs_reader = vcf.Reader(open(wgs_vcf, 'r'))
        for record in wgs_reader:
            key = ('chr' + record.CHROM, record.POS)
            if key in coords_to_snp:
                snp_id = coords_to_snp[key]
                alleles = [str(x) for x in record.ALT]
                alleles.insert(0, record.REF)
                for sample in record.samples:
                    gt = sample['GT']
                    alleles_in_sample = get_alleles_from_genotype(alleles, gt)
                    wgs_results[snp_id] = ''.join(sorted(alleles_in_sample))

        total_snps = 0
        array_calls = 0
        wgs_calls = []
        all_matches = []
        high_quality_matches = []

        # Compare array results to wgs
        for key, value in array_results.items():
            total_snps += 1
            if value['alleles']:
                array_calls += 1  # count of snps genotyped by array
                if key in wgs_results:
                    wgs_calls.append(key)  # list of snps called by wgs
                    wgs_genotype = wgs_results[key]
                    if wgs_genotype == value['alleles']:  # if match
                        all_matches.append(key)
                        if value['quality'] in ['A', 'B']:  # A and B are high quality calls
                            high_quality_matches.append(key)
                    else:
                        mismatch_temp_df = pd.DataFrame(
                            [[sample_name, key, wgs_genotype, value['alleles'], value['quality'], value['filter']]],
                            columns=mismatch_cols
                        )
                        mismatch_df = mismatch_df.append(mismatch_temp_df)

        # calculate probabilities
        all_prob, all_risk = calculate_power(all_matches)
        high_qual_prob, high_qual_risk = calculate_power(high_quality_matches)

        temp_df = pd.DataFrame(
            [[
                sample_name,
                '%s/%s' % (array_calls, total_snps),
                '%s/%s' % (len(wgs_calls), total_snps),
                '%s/%s' % (len(all_matches), len(wgs_calls)),
                all_prob,
                all_risk,
                '%s/%s' % (len(high_quality_matches), len(wgs_calls)),
                high_qual_prob,
                high_qual_risk
            ]],
            columns=main_cols
        )
        main_df = main_df.append(temp_df)

    writer = pd.ExcelWriter(outfile)

    workbook = writer.book
    fail_format = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})

    main_df.to_excel(writer, index=False, sheet_name='IdentityCheck')
    main_ws = writer.sheets['IdentityCheck']
    main_ws.set_column('A:A', 18)
    main_ws.set_column('B:B', 12)
    main_ws.set_column('C:C', 11)
    main_ws.set_column('D:D', 13)
    main_ws.set_column('E:E', 28)
    main_ws.set_column('F:F', 15)
    main_ws.set_column('G:G', 20)
    main_ws.set_column('H:H', 39)
    main_ws.set_column('I:I', 24)
    main_ws.conditional_format(
        'D2:D%s' % (len(infiles) + 1),
        {'type': 'formula', 'criteria': '=IF(LEFT(D2,SEARCH("/",D2)-1)/MID(D2,SEARCH("/",D2)+1,99)<1,TRUE,FALSE)',
         'format': fail_format}
    )  # highlight cells in red where number of matches < number of shared snp calls

    mismatch_df.to_excel(writer, index=False, sheet_name='Mismatches')
    mismatch_ws = writer.sheets['Mismatches']
    mismatch_ws.set_column('A:A', 18)
    mismatch_ws.set_column('B:B', 10)
    mismatch_ws.set_column('C:C', 15)
    mismatch_ws.set_column('D:D', 22)
    mismatch_ws.set_column('E:E', 16)
    mismatch_ws.set_column('F:F', 15)

    writer.save()

    # move files to archive once processed
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        for s in all_samples:
            archive_files(s)

pipeline_run()

