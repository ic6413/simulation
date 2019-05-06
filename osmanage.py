import os

def create_directory(post_process_folder_path):
    # define the name of the directory to be created
    if not os.path.isdir(post_process_folder_path):
        try:  
            os.mkdir(post_process_folder_path)
        except OSError:  
            print ("Creation of the directory %s failed" % post_process_folder_path)
        else:  
            print ("Successfully created the directory %s " % post_process_folder_path)