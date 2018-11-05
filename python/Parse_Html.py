from bs4 import BeautifulSoup
import re

def Parse_Html(file):
    name = []
    all_tax = []
    acc = []
    length = []
    with open(file, "r") as f:
        output = f.read()
    soup = BeautifulSoup(output, features="html.parser")
    for entry in soup.find_all('docsum'):
        try:
            name.append(re.search('\[(.*)\]',entry.find('item', {'name': 'Title'}).text).group(1))
        except:
            continue
        try:
            all_tax.append(entry.find('item', {'name': 'TaxId'}).text)
        except:
            continue
        try:
            acc.append(entry.find('item', {'name': 'AccessionVersion'}).text)
        except:
            continue
        try:
            length.append(entry.find('item', {'name': 'Length'}).text)
        except:
            continue

    return name, all_tax, acc, length

def Write_Html_To_File (file):
    a,b,c,d = Parse_Html(file)
    with open(file[:-4] + '_parsed.csv', 'w') as f:
        f.write('Name,Taxon_ID,Accession_Number,Length\n')
        for i in range(len(a)):
            line = a[i].strip() + ',' + b[i].strip() + ',' + c[i].strip() + ',' + d[i].strip() + '\n'
            f.write(line)
