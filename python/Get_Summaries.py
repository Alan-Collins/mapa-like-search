import urllib.request
import os
import errno

def get_summaries(directory, name):

    input_dir = directory
    output_dir = str(os.path.abspath(os.path.join(input_dir ,"..") + "/" + name))
    try:
        os.mkdir(output_dir) # Create folder to contain split files
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    total = len(os.listdir(input_dir))-1
    i = 1
    for file in os.listdir(input_dir):
        if file.endswith('.txt'):
            listACC = ""
            with open(os.path.join(input_dir, file)) as f:
                for line in f.readlines():
                    listACC += line.strip() + ","
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=docsum&id=' + listACC
            response = urllib.request.urlopen(url)
            text = response.read().decode('utf-8')

            print("working on file %s/%s" %(str(i), str(total)))
            i+=1

            with open(output_dir + "/" + file[:-4] + "_html_sum.txt", 'w') as out:
                out.writelines(text)
                out.close()

get_summaries("/Users/alancollins/GitHub/File_Management/test", "test_lapG_ACC_html_summaries")

