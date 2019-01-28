import interpretation_requests as interp_r
import opencga
import os

class DownloadInfo:
	"""
	Class for storing information about a case for which the VCF is to be downloaded
	"""
	def __init__(self, proband_id, ir_id, d_number):
		self.d_number = d_number
		self.proband_id = proband_id
		self.ir_id = ir_id
		self.json = self.get_json()
		self.build = self.json['assembly']
		self.sample_type = self.json['sample_type']
		self.study_id = opencga.get_study_id(self.sample_type, assembly=self.build)
		self.vcf_filename = self.get_vcf_filename()
		self.file_id = opencga.find_file_id(self.study_id, 'VCF', self.vcf_filename)

	def get_json(self):
		"""
		Retrieves are returns the case JSON from the CIP-API 
		"""
		ir_id, ir_version = self.ir_id.split('-')
		return interp_r.get_interpretation_request_json(ir_id, ir_version)

	def get_vcf_filename(self):
		"""
		Retrieves and returns the VCF filename for the proband.
		The format for this is <proband-sample-id>_<participant-sample-id>.vcf.gz
		"""
		participants = self.json['interpretation_request_data']['json_request']['pedigree']['participants']
		for participant in participants:
			if self.proband_id == participant['gelId']:
				vcf_filename = participant['samples'][0] + '_' + participant['samples'][0] + '.vcf.gz'
				return vcf_filename



def get_interp_r_ids(proband_ids):
	'''
	Gets a list of all the interpretation request IDs from given list of family IDs
	Initialises a DownloadInfo object for each case.
	Cases with a blocked status are ignored.
	'''
	ids = ''
	for k,v in proband_ids.items():
		ids += k + ','
	intep_r_list = interp_r.get_interpretation_request_list(proband_id=ids)
	cases = []
	for i in intep_r_list:
		if i['last_status'] != 'blocked':
			cases.append(i)
	results = []
	for req in cases:
		results.append(DownloadInfo(req['proband'],req['interpretation_request_id'],proband_ids[req['proband']]))
	return results

def load_config():
	'''
	Reads in the config file and returns as a dictionary
	'''
	config_dict = {}
	with open("config.txt", 'r') as config_file:
		for line in config_file:
			if not line.startswith('#'):
				line = line.strip().split('=', 1)
				if len(line) == 2:
					config_dict[line[0]] = line[1]
	return config_dict

def read_infile(storage_path):
	'''
	Reads tab delimited input file consisting of GeL participant ID and D number.
	Returns a dictionary of GeL paticipant ID keys with D number values. 
	'''
	with open(storage_path + 'input.txt') as in_file:
		rows = ( line.split('\t') for line in in_file)
		case_dict = {row[0]:row[1].rstrip('\n') for row in rows}
	return case_dict


# Runs the script
config_dict = load_config()
storage_path = config_dict['storage_path']
bcftools = config_dict['bcftools']
case_dict = read_infile(storage_path)
downloadinfo_list = get_interp_r_ids(case_dict)
for downloadinfo in downloadinfo_list:
	opencga.download_file(downloadinfo.file_id, downloadinfo.study_id, downloadinfo.vcf_filename, download_folder=storage_path + 'VCFs/')
	cmd = 'zcat ' + storage_path + 'VCFs/' + downloadinfo.vcf_filename + ' | bgzip -c > ' + storage_path + 'VCFs/' + downloadinfo.vcf_filename + '.new.vcf.gz && tabix ' + storage_path + 'VCFs/' + downloadinfo.vcf_filename + '.new.vcf.gz'
	os.system(cmd)
	if downloadinfo.build == "GRCh37":
		cmd = bcftools + ' view -O v -R ' + storage_path + 'hg19.nochr.txt ' + storage_path + 'VCFs/' + downloadinfo.vcf_filename + '.new.vcf.gz > ' + storage_path + 'Intersected-VCFs/' + downloadinfo.d_number + '_GRCh37.intersected.vcf'
		os.system(cmd)
	elif downloadinfo.build == "GRCh38":
		cmd = bcftools + ' view -O v -R ' + storage_path + 'hg38.nochr.txt ' + storage_path + 'VCFs/' + downloadinfo.vcf_filename + '.new.vcf.gz > ' + storage_path + 'Intersected-VCFs/' + downloadinfo.d_number + '_GRCh38.intersected.vcf'
		os.system(cmd)