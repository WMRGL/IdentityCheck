"""

 Script to coordinate the genearation and upload of vcf snp check files to Genomics England.

 To run, e.g.:

	python generate_and_upload_snp_vcf.py


 Ed Stone 11 March 2020

"""
import os
import re
import json
from xml_parser import ParseMassArrayXml
from vcf_uploader import UploadVcf
from load_config import LoadConfig
from datetime import datetime

class GenerateAndUploadVcf(object):
	def __init__(self):
		self.xml_path = LoadConfig().load()['xml_path']
		self.vcf_path = LoadConfig().load()['vcf_path']
		self.bed_file = LoadConfig().load()['bed_file']
		self.log_file = LoadConfig().load()['log_file']
		self.write_to_log('Starting generation and upload of SNP VCFs')

	def generate_vcfs(self):
		self.write_to_log('Attempting to generate SNP VCF files from MassARRAY XMLs')
		xml_file_found = False
		for filename in os.listdir(self.xml_path):
			path = self.xml_path + filename
			if filename.endswith(".xml"):
				xml_file_found = True
				try:
					parse_xml = ParseMassArrayXml(path, self.bed_file, self.vcf_path)
					parse_xml.write_to_vcf()
					os.rename(path, self.xml_path + "processed/" + datetime.now().strftime("%Y%m%d-%H%M%S-") + filename)
					self.write_to_log('VCF files generated for {} and XML file moved to "processed" folder'.format(filename))
				except:
					os.rename(path, self.xml_path + "failed_to_process/" + datetime.now().strftime("%Y%m%d-%H%M%S-") + filename)
					self.write_to_log('ERROR: Unable to generate VCFs for {}. XML moved to "failed_to_process" folder.'.format(filename))
		if not xml_file_found:
			self.write_to_log('No XML files found in the holding directory')

	def write_to_log(self, message):
		log = '{}\t{}\n'.format(datetime.now(), message)
		with open(self.log_file, "a") as log_file:
			log_file.write(log)

	@staticmethod
	def check_participant_id(participant_id):
		error = None
		if not re.match(r'^p\d{11}$', participant_id.lower()):
			error = 'Incorrect participant ID. Received {} which does not match the required specification.'.format(participant_id)
		return participant_id.lower(), error

	@staticmethod
	def check_referral_id(referral_id):
		error = None
		if not re.match(r'^r\d{11}$', referral_id.lower()):
			error = 'Incorrect referral ID. Received {} which does not match the required specification.'.format(referral_id)
		return referral_id.lower(), error

	@staticmethod
	def check_laboratory_sample_id(laboratory_sample_id):
		error = None
		if not re.match(r'^\d{10}$', laboratory_sample_id):
			error = 'Incorrect laboratory sample ID. Received {}. Should be 10 digits'.format(laboratory_sample_id)
		return laboratory_sample_id, error

	def upload_vcfs(self):
		self.write_to_log('Attempting to upload SNP VCF files to NGIS')
		vcf_uploader = UploadVcf() # Initialising gets authentication token ready to upload VCFs
		if vcf_uploader.token:
			self.write_to_log("Authentication successful")
			vcf_file_found = False
			for filename in os.listdir(self.vcf_path):
				path = self.vcf_path + filename
				if filename.endswith(".vcf"):
					vcf_file_found = True
					filename_split = filename[:-4].split('_')
					if len(filename_split) != 4:
						os.rename(path, self.vcf_path + 'not_uploaded/' + datetime.now().strftime("%Y%m%d-%H%M%S-") + filename)
						self.write_to_log('{} does not match the expected filename convention (d.number_referralID_participantID_sampleID.xml). File moved to "not_uploaded" directory.'.format(filename))
					else:
						errors = []
						referral_id, error = self.check_referral_id(filename_split[1])
						if error:
							errors.append(error)
						participant_id, error = self.check_participant_id(filename_split[2])
						if error:
							errors.append(error)
						laboratory_sample_id, error = self.check_laboratory_sample_id(filename_split[3])
						if error:
							errors.append(error)
						if errors:
							os.rename(path, self.vcf_path + 'not_uploaded/' + datetime.now().strftime("%Y%m%d-%H%M%S-") + filename)
							self.write_to_log('ERROR: {} does not match the expected filename convention (d.number_referralID_participantID_sampleID.xml). File moved to "not_uploaded" directory.'.format(filename))
							for error in errors:
								self.write_to_log('ERROR: {}'.format(error))
						else:
							self.write_to_log('Attempting to upload {}'.format(filename))
							response = vcf_uploader.upload_vcf(referral_id, participant_id, laboratory_sample_id, path)
							if response.status_code == 200:
								if response.json()['status'] == 'success':
									self.write_to_log('{} uploaded. {}'.format(filename, response.json()['data']['message']))
									os.rename(path, self.vcf_path + 'uploaded/' + datetime.now().strftime("%Y%m%d-%H%M%S-") + filename)
									self.write_to_log('{} moved to "uploaded" directory'.format(filename))
								else: 
									self.write_to_log('ERROR: {} not uploaded. {}'.format(filename, response.json()['data']['message']))
									os.rename(path, self.vcf_path + 'upload_failed/' + datetime.now().strftime("%Y%m%d-%H%M%S-") + filename)
									self.write_to_log('{} moved to "upload_failed" directory'.format(filename))
							elif response.status_code == 400:
								self.write_to_log('ERROR: Bad request. Response code: {}.'.format(response.status_code))
								self.write_to_log('ERROR: {} not uploaded. {}'.format(filename, response.json()['data']['message']))
								os.rename(path, self.vcf_path + 'upload_failed/' + datetime.now().strftime("%Y%m%d-%H%M%S-") + filename)
								self.write_to_log('{} moved to "upload_failed" directory'.format(filename))
							else:
								self.write_to_log('ERROR: Unsuccessful response code: {}.'.format(response.status_code))

			if not vcf_file_found:
				self.write_to_log('No VCF files found in the holding directory')
		else:
			self.write_to_log("ERROR: Authentication failed with error code:", upload_vcf.token_response.status_code)



if __name__ == "__main__":
	generate_and_upload = GenerateAndUploadVcf()
	generate_and_upload.generate_vcfs()
	generate_and_upload.upload_vcfs()
	generate_and_upload.write_to_log("Finished generating and uploading VCFs")