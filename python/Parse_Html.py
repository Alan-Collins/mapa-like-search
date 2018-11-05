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
        name.append(re.search('\[(.*)\]',entry.find('item', {'name': 'Title'}).text).group(1))
        all_tax.append(entry.find('item', {'name': 'TaxId'}).text)
        acc.append(entry.find('item', {'name': 'AccessionVersion'}).text)
        length.append(entry.find('item', {'name': 'Length'}).text)

    return name, all_tax, acc, length

def Write_Html_To_File (file):
    a,b,c,d = Parse_Html(file)
    with open(file[:-4] + '_parsed.csv', 'w') as f:
        for i in range(len(a)):
            line = a[i].strip() + ',' + b[i].strip() + ',' + c[i].strip() + ',' + d[i].strip() + '\n'
            f.write(line)

Write_Html_To_File('/Users/alancollins/GitHub/File_Management/test_lapG_ACC_html_summaries/combined.txt')
