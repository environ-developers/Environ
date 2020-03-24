def load_files(file):
    loaded_file = []
    with open(file, 'r') as f:
        while True:
            try:
                line = f.readline()
                loaded_file.append(line.rstrip())
                if line[0] == 'hello':
                    print('wtf')
            except IndexError:
                break
    return loaded_file

def add_environ_dependencies(file_array):
    for i,line in enumerate(file_array):
        try:
            if line.strip().split()[0] == 'MODFLAGS=':
                file_array[i+1] = file_array[i+1]+' $(MOD_FLAG)../../Environ/Modules_Files \\'
                file_array.insert(i+2,'         $(MOD_FLAG)../../Environ/UtilXlib $(MOD_FLAG)../../Environ/FFTXlib \\')
                file_array.insert(i+3,'         $(MOD_FLAG)../../Environ/PW_files')
            elif line.strip().split()[0] == 'QEMODS=../../Modules/libqemod.a':
                file_array[i+2] = file_array[i+2]+' \\'
                file_array.insert(i+3,'''         ../../Environ/src/libenviron.a ../../Environ/UtilXlib/libenvutil.a \\
         ../../Environ/FFTXlib/libenvfft.a ../../Environ/Modules_Files/libenvmod.a ../../Environ/PW_files/libenvpw.a''')
            elif line.strip().split()[0] == 'TLDEPS=':
                file_array[i] = file_array[i]+' libenviron'
        except IndexError:
            continue
    return file_array

if __name__ == "__main__":
    import sys
    import numpy as np
    loaded_file = load_files(sys.argv[1])
    loaded_file = add_environ_dependencies(loaded_file)
    for line in loaded_file:
       print(line)