"""

 Script to parse XML output from MassArray ExomeQC and generate a VCF for each sample in the results.

 To run, e.g.:

    python xml_parser.py -xml WMRGL_ExomeQC_validation_Plate1_031218_1_F0008858.xml -bed hg19.bed


 Sarah Burns 28 Jan 2019

"""

from xml.etree import ElementTree
from datetime import datetime
import re
import argparse
import os


arg_parser = argparse.ArgumentParser(description='Parses XML from MassArray into VCFs for each sample in results.')
arg_parser.add_argument('-xml', action='store', help='Path to XML file', required=True)
arg_parser.add_argument('-bed', action='store', help='Bed file with full SNP info', required=True)
arg_parser.add_argument('-out', action='store', help='Output directory', default='/network/processed/100K_VCFs/limbo/')
args = arg_parser.parse_args()


class ParseMassArrayXml(object):

    def __init__(self, xml_file, snp_bedfile, output_dir):
        self.xml_tree = ElementTree.parse(xml_file)
        self.snp_file = snp_bedfile
        self.output = output_dir

    def get_snps(self):
        """ Convert SNP bed file into dictionary format.

            Output e.g. {
                            'rs123': {
                                'chrom': 'chr1',
                                'pos': 123,
                                'ref': 'A',
                                'alt': 'G',
                                'gene': 'ABC1',
                                'maf': 0.44,
                                'genome_build': 'hg19'
                            }
                        }

        """
        d = {}
        with open(self.snp_file, 'r') as infile:
            for row in infile:
                if row:
                    row_split = row.strip().split('\t')
                    chrom = row_split[0]
                    pos = row_split[1]
                    name = row_split[3].split('|')
                    snp_id = name[0]
                    gene = name[1]
                    ref_allele = name[2]
                    alt_alleles = name[3]
                    freq = name[4]
                    genome = name[5]
                    d[snp_id] = {
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref_allele,
                        'alt': alt_alleles,
                        'gene': gene,
                        'maf': freq,
                        'genome_build': genome
                    }
        return d

    def get_analyser_name(self):
        """ Get name of MassArray analyser to add to VCF header.

            Output e.g. "MA4-3.3.1.26-iPLEX"

        """
        for analyser in self.xml_tree.getroot():
            for child in analyser:
                if child.tag == 'ms_name':
                    return child.text

    def get_pcr_sequences(self):
        """ Get PCR sequences for each SNP in assay and return dictionary.

            Output e.g. { 'rs123': ['GAGA', 'GAGA'] }

        """
        d = {}
        for analyser in self.xml_tree.getroot():
            for child in analyser:
                if child.tag == 'all-assays':
                    for assay in child:
                        attributes = assay.attrib
                        assay_id = attributes['id']
                        if re.match(r'rs\d+', assay_id):
                            d[assay_id] = [attributes['pcr1'], attributes['pcr2']]
        return d

    def get_results(self):
        """ Get SNP results for each sample and return dictionary.

            Output e.g. {
                            'D00.00001': {
                                'rs123': {
                                    'genotype': '0/1',
                                    'quality': 'A'
                                }
                            }
                        }

        """
        d = {}
        for analyser in self.xml_tree.getroot():
            for child in analyser:
                if child.tag == 'all-records':
                    for record in child:
                        attributes = record.attrib
                        sample = attributes['sampleId']
                        assay_id = attributes['assayId']
                        genotype = attributes['genotypeId']
                        quality = attributes['description'].split('.')[0]
                        if re.match(r'rs\d+', assay_id):
                            if sample in d:
                                if assay_id in d[sample]:
                                    for allele in list(genotype):
                                        if allele not in d[sample][assay_id]['genotype']:
                                            d[sample][assay_id]['genotype'] += allele
                                    if quality not in d[sample][assay_id]['quality']:
                                        d[sample][assay_id]['quality'].append(quality)
                                else:
                                    d[sample][assay_id] = {'genotype': genotype, 'quality': [quality]}
                            else:
                                d[sample] = {assay_id: {'genotype': genotype, 'quality': [quality]}}
        return d

    def get_snp_call_rate(self):
        """ Calculate call rates for each SNP and return dictionary. For each SNP, count the number of samples with no
            result, returning fraction.

            Output e.g. { 'rs123': 0.1 }

        """
        d_temp = {}
        results = self.get_results()
        for sample, variants in results.items():
            for snp, info in variants.items():
                if snp in d_temp:
                    d_temp[snp].append(info['genotype'])
                else:
                    d_temp[snp] = [info['genotype']]
        d = {}
        sample_count = len(results.keys())
        for key, value in d_temp.items():
            na_count = len([x for x in value if not x])
            d[key] = float(na_count) / float(sample_count)
        return d

    def write_to_vcf(self):
        """ Run all functions in class to generate VCF file for each sample.

        """

        # 1. Generate header info
        date_for_vcf = datetime.now().strftime('%Y%m%d')
        header_info = [
            '##fileformat=VCFv4.2',
            '##fileDate=%s' % date_for_vcf,
            '##source=%s' % self.get_analyser_name(),
            '##reference=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
            '##contig=<ID=chr1,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz>',
            '##contig=<ID=chr2,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr2.fa.gz>',
            '##contig=<ID=chr3,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr3.fa.gz>',
            '##contig=<ID=chr4,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr4.fa.gz>',
            '##contig=<ID=chr5,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr5.fa.gz>',
            '##contig=<ID=chr6,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr6.fa.gz>',
            '##contig=<ID=chr7,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz>',
            '##contig=<ID=chr8,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr8.fa.gz>',
            '##contig=<ID=chr9,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr9.fa.gz>',
            '##contig=<ID=chr10,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr10.fa.gz>',
            '##contig=<ID=chr11,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr11.fa.gz>',
            '##contig=<ID=chr12,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr12.fa.gz>',
            '##contig=<ID=chr13,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr13.fa.gz>',
            '##contig=<ID=chr14,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr14.fa.gz>',
            '##contig=<ID=chr15,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr15.fa.gz>',
            '##contig=<ID=chr16,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr16.fa.gz>',
            '##contig=<ID=chr17,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz>',
            '##contig=<ID=chr18,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr18.fa.gz>',
            '##contig=<ID=chr19,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr19.fa.gz>',
            '##contig=<ID=chr20,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz>',
            '##contig=<ID=chr21,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz>',
            '##contig=<ID=chr22,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz>',
            '##contig=<ID=chrM,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz>',
            '##contig=<ID=chrX,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrX.fa.gz>',
            '##contig=<ID=chrY,URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrY.fa.gz>',
        ]
        header_parameters = [
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=MTQ,Number=1,Type=String,Description="MassArray Typer quality value for SNP call. '
            'A=Conservative, B=Moderate, C=Aggressive, D=Low Probability, i=Low Intensity. A and B are considered high '
            'quality scores.">',
            '##INFO=<ID=PCR,Number=2,Type=String,Description="PCR sequences used in assay.">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Minor allele frequency from population data.">',
            '##INFO=<ID=Gene,Number=A,Type=String,Description="HGNC Gene Name for gene containing SNP.">',
            '##INFO=<ID=Build,Number=A,Type=String,Description="Genome build used to determine SNP position for assay.">',
            '##FILTER=<ID=LowCallRate,Description="SNP not called in at least 30% of samples in assay.">',
        ]

        # 2. Extract info from XML file
        results = self.get_results()
        snps = self.get_snps()
        pcr_sequences = self.get_pcr_sequences()
        call_rates = self.get_snp_call_rate()

        # 3. For each sample, create VCF, add headers, determine genotype of each SNP and write to file.
        for sample, variants in results.items():

            with open(os.path.join(self.output, '%s.vcf' % sample), 'w+') as outfile:

                header_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', str(sample)]

                outfile.write('%s\n' % '\n'.join(header_info))
                outfile.write('%s\n' % '\n'.join(header_parameters))
                outfile.write('#%s\n' % '\t'.join(header_fields))

                # for each variant, make a line to add to the file which will
                # then be sorted
                lines_to_write = []
                for snp, info in variants.items():

                    ref_allele = snps[snp]['ref']
                    alt_alleles = snps[snp]['alt']
                    alt_list = alt_alleles.split(',')

                    # Genotype formatting matches VCF v4.0 spec where ./. is no call.
                    gt_list = []
                    called_genotype = info['genotype']
                    if not called_genotype:
                        gt_list = ['.', '.']
                    elif len(called_genotype) == 1:
                        called_genotype += called_genotype
                    for allele in list(called_genotype):
                        if allele == ref_allele:
                            gt_list.append(0)
                        else:
                            if allele in alt_list:
                                idx = alt_list.index(allele)
                                gt_list.append(idx + 1)
                            else:
                                raise ValueError(
                                    'Called genotype %s not provided as possible alt in bed file. Sample %s and SNP '
                                    '%s %s.' % (called_genotype, sample, snp, alt_alleles)
                                )
                    gt = '/'.join([str(x) for x in gt_list])

                    # Threshold currently set to 0.3 (70% results have a call).
                    snp_call_rate = call_rates[snp]
                    if snp_call_rate >= 0.3:
                        vcf_filter = 'LowCallRate'
                    else:
                        vcf_filter = 'PASS'

                    snp_pcr_seqs = pcr_sequences[snp]

                    lines_to_write.append(
                        '{chr}\t{pos}\t{id}\t{ref}\t{alt}\t.\t{filter}\tAF={af};PCR={pcr};Gene={gene};Build={build}\t'
                        'GT:MTQ\t{gt}:{qual}\n'.format(
                            chr=snps[snp]['chrom'],
                            pos=snps[snp]['pos'],
                            id=snp,
                            ref=ref_allele,
                            alt=alt_alleles,
                            filter=vcf_filter,
                            af=snps[snp]['maf'],
                            pcr=','.join(snp_pcr_seqs),
                            gene=snps[snp]['gene'],
                            build=snps[snp]['genome_build'],
                            gt=gt,
                            qual=','.join(info['quality'])
                        )
                    )

                sorted_lines_to_write = sorted(
                    lines_to_write,
                    key=lambda x: (
                        # first key for sorting is the int value of chr
                        int(x.split('\t')[0][3:]),
                        # second key for sorting is the position of the variant
                        int(x.split('\t')[1])
                    )
                )

                for line in sorted_lines_to_write:
                    outfile.write(line)

# Create directory for output files and run parser.
xml_parser = ParseMassArrayXml(xml_file=args.xml, snp_bedfile=args.bed, output_dir=args.out)
xml_parser.write_to_vcf()

