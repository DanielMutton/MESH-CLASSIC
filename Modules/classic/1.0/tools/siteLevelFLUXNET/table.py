#!/home/crp102/.conda/envs/python_dev/bin/python3.7
import os
import sys
import yaml

fluxnet_files_path = sys.argv[1]
output = sys.argv[2]



def html_table(lol):
  table = ""
  table += '<table class="classicTable">\n'
  table += '<thead><tr><th>site</th><th>lon</th><th>lat</th><th>start</th><th>end</th><th>notes</th><th>principal references</th></tr></thead>\n'
  for sublist in lol:
    table += '<tr><td>'
    table += '</td>\n<td>'.join(sublist)
    table += '</td></tr>\n'
  table +='</table>\n'
  return table

data = []
for item in os.listdir(fluxnet_files_path):
  info = []
  yamlfile = fluxnet_files_path + "/" + item + "/siteinfo.yaml"
  with open(yamlfile, "r") as f:
    contents = yaml.safe_load(f)
  info.append(contents["name"] + " (" + contents["biome"] + ")")
  info.append(str(contents["lon"]))
  info.append(str(contents["lat"]))
  info.append(str(contents["start"]))
  info.append(str(contents["end"]))
  info.append(contents["notes"].replace("\n", "<br>"))
  info.append(contents["principal_refs"])
  data.append(info)

data.sort(key=lambda x: x[0])

with open(output + "/site_info.html", "w") as f:
  f.write('''<!DOCTYPE html>\n<head><style>table.classicTable {
  border: 1px solid #1C6EA4;
  background-color: #EEEEEE;
  width: 100%;
  text-align: left;
  border-collapse: collapse;
}
table.classicTable td, table.classicTable th {
  border: 1px solid #AAAAAA;
  padding: 5px 5px;
}
table.classicTable tbody td {
  font-size: 13px;
}
table.classicTable tr:nth-child(even) {
  background: #D0E4F5;
}
table.classicTable thead {
  background: #1C6EA4;
  background: -moz-linear-gradient(top, #5592bb 0%, #327cad 66%, #1C6EA4 100%);
  background: -webkit-linear-gradient(top, #5592bb 0%, #327cad 66%, #1C6EA4 100%);
  background: linear-gradient(to bottom, #5592bb 0%, #327cad 66%, #1C6EA4 100%);
  border-bottom: 2px solid #444444;
}
table.classicTable thead th {
  font-size: 15px;
  font-weight: bold;
  color: #FFFFFF;
  text-align: left;
  border-left: 2px solid #D0E4F5;
}
table.classicTable thead th:first-child {
  border-left: none;
}

table.classicTable tfoot {
  font-size: 14px;
  font-weight: bold;
  color: #FFFFFF;
  background: #D0E4F5;
  background: -moz-linear-gradient(top, #dcebf7 0%, #d4e6f6 66%, #D0E4F5 100%);
  background: -webkit-linear-gradient(top, #dcebf7 0%, #d4e6f6 66%, #D0E4F5 100%);
  background: linear-gradient(to bottom, #dcebf7 0%, #d4e6f6 66%, #D0E4F5 100%);
  border-top: 2px solid #444444;
}
table.classicTable tfoot td {
  font-size: 14px;
}
table.classicTable tfoot .links {
  text-align: right;
}
table.classicTable tfoot .links a{
  display: inline-block;
  background: #1C6EA4;
  color: #FFFFFF;
  padding: 2px 8px;
  border-radius: 5px;
}</style></head>\n<body>''')
  f.write(html_table(data))
  f.write('</body>\n</html>')


