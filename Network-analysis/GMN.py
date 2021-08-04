'''
This script generate Graph Neural Network Embeddings in 2D format

@ Haoyang Mi
@ Johns Hopkins University, School of Medicine
'''
import os, errno


for core in range(1,38):

    # community JSON file path
    PATH = 'D:/DP/Projects/HCC/Data/Networks/' + str(core)
    PATH.replace(os.sep, '/')
    print(PATH)

    # number of files

    fileNum = len(os.listdir(PATH))

    # define output path
    OT_PATH = 'D:/DP/Projects/HCC/Data/Features/' + str(core)

    # make feature folder
    try:
        os.makedirs(OT_PATH)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    # define output file
    OT_File = OT_PATH + '/nci.csv'
    OT_File.replace(os.sep, '/')


    os.system('python graph2vec.py --input-path' + ' ' + PATH + '/' + ' ' +'--output-path' + ' ' + OT_File)

