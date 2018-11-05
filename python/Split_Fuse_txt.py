from math import ceil
import os
import errno


def Split_txt(file, piece_size):

    """
    Given a text file and a number this will create a directory and populate it with pieces of the text file that are
    "piece size" lines long.
    :param file: "filename.txt"
    :param piece_size: number of lines per output file
    :return: doesn't return anyrhing, but produces a folder full of files.
    """

    dir_path = os.path.dirname(os.path.realpath(file))
    new_fol = str(dir_path + "/" + file[:-4] + "_split_files")
    try:
        os.mkdir(new_fol) # Create folder to contain split files
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass
    with open(file) as f:
        content = f.readlines()
    length = len(content)
    pieces = ceil(length/piece_size)

    slices = []
    for i in range(pieces):
        if i == pieces - 1:
            slices.append((i*piece_size, length))
        else:
            slices.append((i*piece_size, (i+1)*piece_size))

    for i in range(len(slices)):
        start = slices[i][0]
        end = slices[i][1]
        with open(new_fol + "/" + file[:-4] + str(i) +".txt", 'w') as out:
            out.writelines(content[start:end])
            out.close()






def fuse_txt (dir, split):

    """
    Given a directory full of files, fuses them together and outputs a single file called "combined.txt"
    Adds whatever arg is given for split at the points where the files were fused.
    :param dir: directory containing text files
    :param split: something to add at the fusion points e.g. '/n'
    :return: returns nothing. Makes a big text file.
    """

    bigfile = ''
    for file in os.listdir(dir):
        if file.endswith('.txt') and file != "combined.txt":
            with open(os.path.join(dir, file)) as f:
                for line in f.readlines():
                    bigfile += line
                f.close()
            bigfile += split

    with open(os.path.join(dir, "combined.txt"), 'w') as f:
        f.writelines(bigfile)
        f.close()


fuse_txt("/Users/alancollins/GitHub/File_Management/lapG_ACC_html_summaries", "")