
import gspread
from oauth2client.service_account import ServiceAccountCredentials
from pdb import set_trace as stop
import pprint
import os
from vison.support import utils

scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']

credsf = os.path.join(utils.credentials_path,
                      'EuclidVis_calcampreporter.json')
creds = ServiceAccountCredentials.from_json_keyfile_name(credsf, scope)

client = gspread.authorize(creds)

sheet = client.open('TestB_5APR18').sheet1

pp = pprint.PrettyPrinter()
result = sheet.row_values(1)
#print result

pp.pprint(result)

# sheet.update_cell(1,1,"Ni3hao3!")
#result2 = sheet.row_values(1)
# pp.pprint(result2)

#row = ["I'm","Updating","a","spreadhseet","from","python"]
#index = 2
# sheet.insert_row(row,index)
