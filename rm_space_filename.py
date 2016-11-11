"""
Name: rm_space_filename.py
Created by: Tovah Markowitz
Purpose: to replace all spaces in filenames with underscores
Checks all files within current working directory
"""
##########################
# MODULES
import os
##########################

# list files in folder
files = os.listdir(".")

# find files with spaces in name
filesSpace = [file for file in files if file.find(" ") is not -1]

# replace spaces with underscores
# and rename file
for fileName in filesSpace:
    name = fileName.replace(" ","_")
    os.rename(fileName,name)
