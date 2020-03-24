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
            if line.strip().split()[0] == 'MODFLAGS=$(BASEMOD_FLAGS)':
                file_array[i+1] = file_array[i+1]+' $(MOD_FLAG)../Environ/Modules_Files \\'
                file_array.insert(i+2,'         $(MOD_FLAG)../Environ/UtilXlib $(MOD_FLAG)../Environ/FFTXlib')
            elif line.strip().split()[0] == 'TLDEPS=libiotk':
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