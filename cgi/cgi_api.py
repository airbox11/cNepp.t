# https://www.cancergenomeinterpreter.org/rest_api
# https://www.cancergenomeinterpreter.org/formats

import requests
import pdb
import os
from zipfile import ZipFile
import time
import sys
import os
import pdb


headers = {'Authorization': 'lvyhwind@hotmail.com 828adf7a9f548e3d6dc4'}
headers = {'Authorization': 'lvyanhongkao@gmail.com 1b707b6766a037c8a299'}

headers = {'Authorization': 'lvyhwind@hotmail.com 828adf7a9f548e3d6dc4'}


def job_launch (reference):
	payload = {'cancer_type': f'{cgiTumorType}', 'title': f'{reference}_{cgiTumorType}', 'reference': f'{reference}'}
	r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
			headers=headers,
			files={
				'mutations': open(f'{file_mutation}', 'rb')
				},
			data=payload,
			verify=False
			)
	# breakpoint()
	job_id = r.json()
	print(f'===[ cgi launch: {job_id}')
	return(r.json())

def job_access_info(job_id):
	payload={'action':'logs'}
	r = requests.get(f'https://www.cancergenomeinterpreter.org/api/v1/{job_id}', headers=headers, params=payload)
	print(r.json())
	status = r.json()['status']
	return(status)


def job_download(job_id):
	zip1 = f'cig_{job_id}.zip'

	payload={'action':'download'}
	r = requests.get(f'https://www.cancergenomeinterpreter.org/api/v1/{job_id}', headers=headers, params=payload)
	with open(zip1, 'wb') as fd:
		fd.write(r._content)

	with ZipFile(zip1, 'r') as zObject:
		zObject.extractall('./cgi')
	os.remove(zip1)

def job_log (job_id):
	payload={'action':'logs'}
	r = requests.get(f'https://www.cancergenomeinterpreter.org/api/v1/{job_id}', headers=headers, params=payload)
	dic_log = r.json()
	with open(f'cgi_{cgiTumorType}_{job_id}.log', 'w') as f: 
		for i in dic_log['logs']:
			i.strip('\n')
			f.write(f'{i}\n')

def job_delete (job_id):
	r = requests.delete(f'https://www.cancergenomeinterpreter.org/api/v1/{job_id}', headers=headers)
	print(r.json)

def get_identifier():
	r = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
	print(r.json())
	return(r.json())

def jobs_remove_all ():
	print('=== === remove older cgi tasks:')
	jobs = get_identifier()
	len1 = len(jobs)
	if len1 > 10:
		for job_id in jobs[:(len1-10)]:
			job_delete(job_id)
	get_identifier()

def job_wait_download (job_id):
	status = job_access_info(job_id) 
	if status == 'Done':
		job_download(job_id)
		job_log(job_id)
	elif status == 'Running' or status == 'Waiting':
		time.sleep(10)
		job_wait_download(job_id)
	elif status == 'Error':
		job_log(job_id)
		return(-1)



## main
workDir = sys.argv[1]
	# workDir = '/omics/groups/OE0422/internal/yanhong/all_in_one_pipeline_collection/mhc4.1/CA59/1_hla_type'
cgiTumorType = sys.argv[2].upper()
os.chdir(workDir)
vcfOnly = os.environ['vcfOnly']


file_mutation = 'mutation.tsv'
if vcfOnly=='promise':
	hg = 'hg38'
else:
	hg = 'hg19'
job_id = job_launch(hg)
job_wait_download(job_id)
jobs_remove_all()
print('===] cgi done')