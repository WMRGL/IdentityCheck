"""

 Script to upload vcf snp check files to Genomics England.

 To run, e.g.:

    python vcf_uploader.py -vcf path/to/vcf -p participant_id -r referral_id -l lsid


 Ed Stone 4 March 2020

"""
import argparse
import requests
import json
from load_config import LoadConfig

class UploadVcf(object):
    def __init__(self):
        self.token, self.token_response = self.get_token()

    def get_token(self):
        # Authentication credentials
        tenant_id = LoadConfig().load()['tenant_id']
        client_id = LoadConfig().load()['client_id']
        client_secret = LoadConfig().load()['client_secret']
        url = "https://login.microsoftonline.com/{tenant_id}/oauth2/token".format(tenant_id=tenant_id)
        payload = "grant_type=client_credentials"
        headers = {
            'Content-Type': "application/x-www-form-urlencoded",
            }
        response = requests.request("POST", url, data=payload, headers=headers, auth=(client_id, client_secret))
        token = None
        if response.status_code == 200:
            data = json.loads(response.text)
            token = data['access_token']
        return token, response

    def upload_vcf(self, r_id, p_id, ls_id, vcf_path):
        host = LoadConfig().load()["snp_host"]
        with open(vcf_path) as f:
            vcf = f.read()
        url = '{host}/sample_vcf/{referral_hrid}/{patient_hrid}/{lims_sample_id}'.format(host=host,
                                                                                         referral_hrid=r_id,
                                                                                         patient_hrid=p_id,
                                                                                         lims_sample_id=ls_id)
        auth_header = {'Authorization': 'Bearer {}'.format(self.token), 'Content-Type': 'text/plain'}
        tv = requests.put(url, data=vcf, headers=auth_header)
        return tv


if __name__ == "__main__":
    # parse arguments
    arg_parser = argparse.ArgumentParser(description='Uploads SNP identity check VCF to Genomics England.')
    arg_parser.add_argument('-vcf', action='store', help='Path to VCF file', required=True)
    arg_parser.add_argument('-p', action='store', help='Participant ID', required=True)
    arg_parser.add_argument('-r', action='store', help='Referral ID', required=True)
    arg_parser.add_argument('-l', action='store', help='LSID', required=True)
    args = arg_parser.parse_args()
    upload_vcf = UploadVcf()
    if upload_vcf.token:
        result = upload_vcf.upload_vcf(args.r, args.p, args.l, args.vcf)
        print("Upload response code:", result.status_code)
        print(result.json())
    else:
        print("Authentication failed with error code:", upload_vcf.token_response.status_code)