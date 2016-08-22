# tdo-demod
### Frequency demodulation for TDO in pulsed fields
-------------------------
[`demodulator.py`](./demodulator.py) is a collection of Python functions for performing frequency demodulation and other initial processing steps for TDO data taken in pulsed magnetic fields. Only two data formats are supported at the moment (see [ExampleData](./ExampleData/)), but the demodulation and spline fitting routines are general-purpose, meaning that one only needs to load data into the correct dictionary form in order to process it using `demod()`, `spline_plot()`, and `data_save()`. That form is:
`osc_dict = {'file': ..., 'dt_fast': ..., 'tdo': ..., 'dt_slow': ..., 'Field': ...}`, where `file` is a string, `dt_fast` and `dt_slow` are doubles, and `tdo` and `Field` are numpy arrays.

I have included an [example IPython notebook](./Demod_Example.ipynb) and some [example data](./ExampleData/) to demonstrate usage.

I have also included [`demod.m`](./demod.m), the original MATLAB script I wrote to demodulate data having the format of the [included LANL data](./ExampleData/LANL_example.txt). It's worse than `demodulator` in pretty much every way, but may be useful to someone.

Author: Logan Bishop-Van Horn, logan.bvh@gmail.com
