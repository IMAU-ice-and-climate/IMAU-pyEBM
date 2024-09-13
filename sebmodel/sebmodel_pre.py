import shutil
import sys, getopt
import os


def main(argv):
    """ SEB model prepare function, which makes the necessary files before running the model.
    Params 
    ======
    chstation: name of the station (command line -aws)
    version: 1D or 2D (command line -v)
    convert: swith to convert info file from fortran SEB model (command line -c)
    Returns
    ======
    Copy of info file in source directory and ready-to-use forcing
    """
    version = ''
    convert = False

    # Get command line arguments
    opts, args = getopt.getopt(argv,"aws:v:c:",["aws=","version="])
    for opt, arg in opts:
        if opt in ("-v", "--version"):
            version = arg
        elif opt in ("-aws", "--aws"):
            chstation = arg
        elif opt in ("-c", "--c"):
            convert = True
    print ('chstation is ', chstation)
    print ('version is ', version)
    print ('convert info file is ', convert)
    
    # Prepare output folder
    ifolder = '../input/' + chstation + '/'
    ofolder = '../output/' + chstation + '/'
    output_netcdf = chstation + '_out.nc'
    data_path = '../output/' + chstation + '/'

    if not os.path.isdir(data_path):
        os.mkdir(data_path)
    elif os.path.exists(output_netcdf):
        os.remove(output_netcdf)
    
    # Copy parameters file
    if (convert == 1):
        print('Conversion of INFO file not implemented')
    infofile = ifolder + chstation +  '_ebm.py'
    tmpinfo = 'info.py'
    shutil.copyfile(infofile, tmpinfo)

    # Add version and chstation to info file
    with open(tmpinfo, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        line = 'chstation =  "' + str(chstation) + '"' + '\n' + 'version =  "' + str(version) + '"'
        f.write(line.rstrip('\r\n') + '\n' + content)
        
    # Prepare version dependant input
    if version == '2D':
        # Downscale forcing
        exit()
    else:  
        # do nothing
        exit()
        

''' MODEL EXECUTION '''

if __name__ == "__main__":
    
   main(sys.argv[1:])