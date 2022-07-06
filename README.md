# Proton Bragg Peak Fit
[![View Fit Proton Bragg Peaks on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/100516-fit-proton-bragg-peaks)

- Provide 3 different IDD measurement data, data acquired from Zebra and other MLIC which can be used in Proton Radiographics
- Bortfeld function implemented in C++, provide IDD, mean gradient and jacobian outputs.
- Input Integrated Depth Dose(IDD) is suggested to rescale to [0,10] or nomalize to [0,1]
- Compiled with Visual Studio + Intel OneAPI, faster than mex in MATLAB. Highly recommanded to compile the src with VS and Intel OneAPI

How to use:  
- If your can run the first statement in Windows, do nothing. If you can't, run ./src/compile.m
- `output = bf_mex((1:64)*0.3,[15,0.3,1e-3,0.4, 12,0.4,1e-3,0.4],'idd')`
- `[x,idd_o] = precise_fit(z,idd_i,num_bp,strict);`
- `x = fast_fit(z,idd_i,num_bp);`

Details can be found in demo.m, column 1D array is prefered such as, zeros(n,1)

Zebra data:

<img src="./Zebra_fit.png" width="500" height="500">

Multi-Bragg Peaks Fit:

<img src="./IDD_fit.png" width="500" height="500">

It take 160s to fit 11057 IDD curves with 2 bragg peak model @i9-9900k

It take 120s to fit 10498 IDD curves with 2 bragg peak model @i9-9900k

Brief introduction:

Bortfeld function is an analytical approximation of the Bragg curve for therapeutic proton beams, given by

$$
\begin{align}
D(z) \approx 
\begin{cases}
\hat{D}(z) \;\;z < R_0 - 10\sigma\\
D(z) \;\; R_0 - 10\sigma\le z\le R_0+5\sigma \\
0 \;\; otherwise
\end{cases}
\end{align}
$$

z denotes the depth in cm. there are 4 parameters in bortfeld funtion $R_0, \sigma, \epsilon, \Phi_0$ and we can guess a good initial points from the table provided by bortfeld.
```
z = (1:64)*0.291; % depth in cm
[vmax,idx] = maxk(idd,k);
R0 = z(idx);% Range
alpha = 0.0022;
p = 1.77;
E0 = (R0./alpha).^(1/p);% estimated proton energy
sigma = sqrt((0.012.*R0.^0.935).^2 + (0.01.*E0).^2.*(alpha.*p.*E0.^(p-1)).^2);
epsilon = 1e-3;
Phi = zv.*zr.*epsilon;
```
or use the simple version
```
R0 = z(idx);% Range
sigma = 0.07*zr;
epsilon = 1e-3;
Phi = zv.*zr.*1e-2;
```
The depth-dose distribution in water is given by $\hat D_{H_2O}(z)$ and $D_{H_2O}(z)$:

$$
\begin{align}
\hat D_{H_2O}(z) &= \frac{\Phi_0}{1+0.012R_0}\left[17.93(R_0-z)^{-0.435}+\left(0.444+31.7\frac{\epsilon}{R_0}\right)(R_0-z)^{0.565} \right]\\
D_{H_2O}(z) &= \Phi_0\frac{e^{-\frac{(R_0-z)^2}{4\sigma^2}}\sigma^{0.565}}{1+0.012R_0} \left[11.26\frac{\mathfrak{D}(-0.565,-\frac{R_0-z}{\sigma})}{\sigma}+ \left(0.157+11.26\frac{\epsilon}{R_0}\right) \mathfrak{D}(-1.565,-\frac{R_0-z}{\sigma})\right]\\
\end{align}
$$

$\mathfrak{D}(a,x)$, is a [parabolic cylinder function](https://mathworld.wolfram.com/ParabolicCylinderFunction.html) defined in Eq.33 We can get $\mathfrak{D}(a,x)=\mathit{U}(-a-0.5,x)$, this function is defined in https://github.com/civerjia/Parabolic-Cylinder-Functions-C-

Parameters of a single bortfeld function is a 4-element 1d array, $[R_0,\sigma,\epsilon,\Phi_0]$,  n-bortfeld function is a $4n$ 1d array,$[R_1,\sigma_1,\epsilon_1,\Phi_1,R_2,\sigma_2,\epsilon_2,\Phi_2,\cdots,R_n,\sigma_n,\epsilon_n,\Phi_n]$.


Reference :
- An analytical approximation of the Bragg curve for therapeutic proton beams
- E. Cojocaru. Parabolic Cylinder Functions (https://www.mathworks.com/matlabcentral/fileexchange/22620-parabolic-cylinder-functions)

