# Running massive tests: automatic testing (local)

Follow these steps in order to run automatic testing.

1. In the terminal, inside the container, enter into the `src` folder.
2. `make clean`
3. In the Makefile, comment and uncomment all the necessary sections (they are indicated)\
considering whether the serial or parallel version needs to be run.
4. In the `main_<version>.cpp`, make sure that the following lines are present:\
`OutputObject result_<function_name> = pso_<version>(<function_name>, dim, stop, n_points);`\
`result_<function_name>.output_to_file();`\
where `<function_name>` takes the name of all the functions we want to test onto.\
Moreover, make sure `result_<function_name>.terminal_info()` is commented out, otherwise a cascade\
of comments will appear while running tests.
5. `make`
6. Run (still inside the container) `./run_tests.sh --n_tests <int1> <int2> <...> --dim <int3> <int4> <...> --n_points <int5> <int6> <...> --max_iter <int7> <int8> <...> --error_tol <float1> <float2> <...> --n_cores <int8> <int9> <...>`\
For example, `./run_tests.sh --n_tests 5 --dim 3 4 --n_points 50 100 --max_iter 100 --error_tol 0.001 --n_cores 2 4`. The command is one for both serial and parallel versions.
7. In the terminal, outside the container, in folder src, run `python3 run_plots.py`. This uses some packages, which can be installed by using the package manage poetry. The dependencies are just numpy and matplotlib, you might have them in some other environment.
8. You can find all the output and plotting results in:
- `src/tests/<function_name>/<dim>/<n_points>/<n_cores>/test_<test_number>.txt`
- `src/plots/<function_name>/agg_type_<index>_<fixed_parameters>.png`\
Plot types are four and are identified in the `index` of the plot file name:
 0. $\Delta_x$ over it_n
 1. final_t over n_cores
 2. final_t over dim
 3. it_n over dim\
All the other relevant configuration parameters are fixed and specified in the plot file name in `<fixed_parameters>`. The results are averaged across `n_tests` (which is nothing but the number of identical tests to be run) and across the irrelevant configuration parameters as well.


