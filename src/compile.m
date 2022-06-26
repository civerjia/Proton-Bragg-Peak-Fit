use_openmp = 1;
use_avx = 1;
is_avx512 = 0;
output_filename = 'pcf';
src_path = './src/*.cpp';
header_flag = ['-I','./includes/'];
output_dir = '..';
compile_mex(output_filename,output_dir,src_path,header_flag,use_openmp,use_avx,is_avx512);
%%
function compile_mex(output_filename,output_dir,src_path,header_flag,use_openmp,use_avx,is_avx512)
% output_filename : filename, (.mexw64) extension is not needed,
% output_dir : '.' or '..'
% src_path = './src/*.cpp';
% header_flag = ['-I','./includes/'];

% use_openmp = 1;
% use_avx = 1;
avx512 = is_avx512 & use_avx;
if use_openmp
    if ismac
        % Code to run on Mac platformmex 
        cpp_openmp_flag = "CXXFLAGS='$CXXFLAGS -fopenmp'";
        LinkerOptimizationFlags = 'LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp -O3"';
        mex('-output',output_filename,'-outdir',output_dir,cpp_openmp_flag,LinkerOptimizationFlags,src_path,header_flag);
    elseif isunix
        % Code to run on Linux platform
        cpp_openmp_flag = "CXXFLAGS='$CXXFLAGS -fopenmp'";
        LinkerOptimizationFlags = 'LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp -O3"';
        if use_avx
            if avx512
                avx_flag = 'CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -mavx512f"';
            else
                avx_flag = 'CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -mavx2"';
            end
        else
            avx_flag = '';
        end
        mex('-output',output_filename,'-outdir',output_dir,cpp_openmp_flag,LinkerOptimizationFlags,avx_flag,src_path,header_flag);
    elseif ispc
        % Code to run on Windows platform
        warning('Compiler should be Visual Studio on Windows!');
        cpp_openmp_flag = 'COMPFLAGS="$COMPFLAGS /openmp"';
        if use_avx
            % visual studio is using /arch:AVX2 or /arch:AVX512
            if avx512
                warning('Be careful! It can be compiled even your CPU do not support AVX512');
                cpp_openmp_flag = 'COMPFLAGS="$COMPFLAGS /openmp /arch:AVX512"';
            else
                cpp_openmp_flag = 'COMPFLAGS="$COMPFLAGS /openmp /arch:AVX2"';
            end
        end
        mex('-output',output_filename,'-outdir',output_dir,cpp_openmp_flag,src_path,header_flag);
    end
else
    % without any feature
    % add -v, display details
    mex('-v','-output',output_filename,'-outdir',output_dir,src_path,header_flag);
end
end
