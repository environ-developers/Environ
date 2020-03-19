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
            if line.strip().split()[0] == 'pwlibs:':
                file_array[i] = file_array[i] + ' envlibs'
            elif line.strip().split()[0] == 'bindir':
                file_array.insert(i+2,'	( cd Environ ; $(MAKE) TLDEPS= all || exit 1 )')
                file_array.insert(i+2,'envlibs :')
                file_array.insert(i+2,'')
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