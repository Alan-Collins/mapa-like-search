import Split_Fuse_txt
import Parse_Html
import Get_Summaries
import os

def run_all (input_name):

    input_file = input_name + ".txt"

    Split_Fuse_txt.Split_txt(input_file, 50)

    split_dir = str(os.path.dirname(os.path.realpath(input_file)) + "/" + input_name + "_split_files")

    Get_Summaries.get_summaries(split_dir, input_name + "_html_summaries")

    summaries_dir = str(os.path.abspath(os.path.join(split_dir ,"..")) + "/" + input_name + "_html_summaries")

    Split_Fuse_txt.Fuse_txt(summaries_dir, "", input_name)

    Parse_Html.Write_Html_To_File(summaries_dir + "/" + input_name + "_combined.txt")

# Retrieved list of accession numbers for LapD-like and LapG-like proteins from these two URLs respectively:
# https://www.ncbi.nlm.nih.gov/sites/entrez?LinkName=cdd_protein&db=cdd&cmd=Link&from_uid=318615 on 11/5/2018
# https://www.ncbi.nlm.nih.gov/sites/entrez?LinkName=cdd_protein&db=cdd&cmd=Link&from_uid=310549 on 11/5/2018

run_all("lapG_ACC")
run_all("lapD_ACC")







