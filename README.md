Multi-block-LBM-GPU

 Please refer the paper: 
 Zhang, Ya, Guang Pan, and Qiaogao Huang. "ICCM2016: The implementation of two-dimensional multi-block lattice Boltzmann method on GPU." International Journal of Computational Methods 16, no. 05 (2019): 1840002.
 DOI: 10.1142/s0219876218400029
  
 How to compile (assume cuda is installed properly, 
 https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/index.html
 https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html)
    nvcc -lcusparse -lcurand -lcublas Multi-block-LBM.cu -o test
  OR Visual Studio on Windows
 	 1. Add the above libraries (cusparse.lib, curand.lib, cublas.lib) at 
 		Project Property Pages -> Linker -> Input -> Additional Dependencies;
 	 2. Build the project;
 	 3. Please use Nsight to debug the code 
		(https://developer.nvidia.com/nsight-visual-studio-edition).
 