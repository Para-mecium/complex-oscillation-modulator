# complex-oscillation-modulator

## Dependency
Symbolic Toolbox (MATLAB)

## Syntax
[solution,Error] = FMAM_ODE(
    system, observables, initialTS, initialParam, 
    paramControlled, propertyPerturb, target, 
    maxstepsizes, errorBar, truncationOrder, ptNum, PV)

## Description
This project is built to solve the complex oscillation modulation problem raised in article *Modulatability of complex oscillators*. The function 'FMAM_ODE' computes the modulation path along which the properties of the periodic solution move towards the target value. 

## Input arguments
**system**: Functions of the ordinary differential equation, specified as a $1\times N$ cell where $N$ is the dimension of the ODE. Each components of 'system' are expected to be function handles.

**observables**: Functions that maps the state variables of ODE to the observables, specified as a $1\times n$ cell where $n$ is the number of observables. Each components of 'observables' are expected to be function handles.

**initialTS**: Time series of the initial periodic solution, specified as a matrix with $N+1$ columns. The first column saves the discrete time points from $0$ to $T$ where $T$ is the initial period and the rest columns save the values of each state variables at corresponding time points, respectively.

**initialParam**: Initial parameters, specified as a $1\times m$ vector.

**paramControlled**: Modulation parameters, specified as a $2\times n_1$ cell. For each column, enter sequentially 'parameter' and the index of the modulation parameter. 

**propertyPerturb**: Properties that are modulated, specified as a cell with $n_1$ columns. For a single column, the first row specifies the class of modulated property. We use'state' for state variable, 'observable' for observable. The second row specifies the index of the component. The thrid row specifies the property: 1 for amplitude, 2 for maximum, 3 for minimum and 4 for phase. 

To modulate the period, add the column with 'parameter' and 0.

**target**: The values of target properties, specified as a $n_1\times 1$ vector. The order should be consistent with the one in **propertyPerturb**.

**maxstepsize**: The max stepsize in a single iteration, specified as a $n_1\times 1$ vector. The order should be consistent with the one in **propertyPerturb**.

**errorBar**: The values of error, specified as a $n_1\times 1$ vector. The function while keep iterating until the errors are smaller than the values in errorBar. The order should be consistent with the one in **propertyPerturb**.

**truncationOrder**: The truncation order of Fourier series, specified as an positive interger. Increasing this number can help improve the accuracy but may also slow down the computation process.

**ptNum**: The number of discrete time points in the pseudo time domain, specified as a positive integer.

**PV**: The primary variable which is reformed into the normal cosine form through the nonlinear time transformation, specified as a $2\times 1$ cell. The first column specifies the class of the variable: 'state' for state variable, 'observable' for observable. The second row specifies the variable index.