import requests
headers = {'Authorization': 'lvyhwind@hotmail.com 828adf7a9f548e3d6dc4'}
payload={'action':'logs'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/job_id', headers=headers, params=payload)
r.json()