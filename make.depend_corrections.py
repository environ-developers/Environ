Modules_Files = ['cell_base.o', 'constants.o', 'control_flags.o', 'fft_base.o', 
            'invmat.o', 'io_files.o', 'io_global.o', 'ions_base.o', 'kind.o',
            'mp_bands.o', 'mp_images.o', 'parameters.o', 'parser.o', 'paw_variables.o',
            'pseudo_types.o', 'radial_grids.o', 'random_numbers.o', 'recvec.o',
            'solvent_tddfpt.o', 'uspp.o', 'wrappers.o', 'ws_base.o']

PW_files = ['esm.o', 'extfield.o', 'lsda_mod.o', 'martyna_tuckerman.o',
            'scf_mod.o', 'vlocal.o', 'v_of_rho.o']

UtilXlib_files = ['clocks_handler.o', 'cuda_util.o', 'data_buffer.o', 'divide.o',
            'error_handler.o', 'find_free_unit.o', 'fletcher32_mod.o', 'hash.o',
            'mem_counter.o', 'mp.o', 'mp_bands_util.o', 'mp_base.o', 'mp_base_gpu.o',
            'parallel_include.o', 'thread_util.o', 'util_param.o']

FFTXlib_files = ['fft_error.o', 'fft_fwinv.o', 'fft_ggen.o', 'fft_helper_subroutine.o',
            'fft_interfaces.o', 'fft_interpolate.o', 'fft_parallel.o',
            'fft_param.o', 'fft_scalar.ARM_LIB.o', 'fft_scalar.DFTI.o',
            'fft_scalar.ESSL.o', 'fft_scalar.o', 'fft_scalar.FFTW.o',
            'fft_scalar.FFTW3.o', 'fft_scalar.SX6.o', 'fft_smallbox.o', 'fft_smallbox_type.o', 
            'fft_support.o', 'fft_types.o', 'fftw_interfaces.o', 'scatter_mod.o',
            'stick_base.o', 'test.o', 'test0.o', 'tg_gather.o']

def load_files(file):
    loaded_file = []
    with open(file, 'r') as f:
        while True:
            try:
                line = f.readline().strip().split()
                loaded_file.append([line[0],line[2]])
            except IndexError:
                break
    return loaded_file

def load_Environ_inc(file):
    loaded_file = []
    with open(file, 'r') as f:
        line = f.readline().strip().split()
        while True:
            try:
                line = f.readline().strip().split()
                if line[1] == '\\' :
                    loaded_file.append(line[0])
            except IndexError:
                loaded_file.append(line[0])
                break
    return loaded_file

if __name__ == "__main__":
    import sys
    import numpy as np
    loaded_file = load_files(sys.argv[1])
    Environ_src_files = load_Environ_inc(sys.argv[2])
    #initial removal
    for line in loaded_file:
        if line[0] in Environ_src_files and not line[1] in Environ_src_files:
            line[1] = line[1].replace("@env_","")
            line[1] = line[1].replace("@",".o")
            #print(line[1])
            if line[1] in Modules_Files:
                line[1] = '../Environ/Modules_Files/'+line[1]
            elif line[1] in PW_files:
                line[1] = '../Environ/PW_Files/'+line[1]
            elif line[1] in UtilXlib_files:
                line[1] = '../Environ/UtilXlib/'+line[1]
            elif line[1] in FFTXlib_files:
                line[1] = '../Environ/FFTXlib/'+line[1]
            elif line[1] == 'mp.o':
                line[1] = '../Environ/UtilXlib/'+line[1]
            elif line[1] == 'kinds.o':
                line[1] = '../Environ/Modules_Files/kind.o'
            else:
                loaded_file.remove(line)
    #double check
    for line in loaded_file:
        if '@' in line[1]:
            line[1] = line[1].replace("@env_","")
            line[1] = line[1].replace("@",".o")
            if line[1] in Modules_Files:
                line[1] = '../Environ/Modules_Files/'+line[1]
            elif line[1] in PW_files:
                line[1] = '../Environ/PW_Files/'+line[1]
            elif line[1] in UtilXlib_files:
                line[1] = '../Environ/UtilXlib/'+line[1]
            elif line[1] in FFTXlib_files:
                line[1] = '../Environ/FFTXlib/'+line[1]
            elif line[1] == 'mp.o':
                line[1] = '../Environ/UtilXlib/'+line[1]
            elif line[1] == 'kinds.o':
                line[1] = '../Environ/Modules_Files/kind.o'
    for line in loaded_file:
      print(line[0]+' : '+line[1])