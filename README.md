# ROMSPY: Regional Ocean Modelling Systems preprocessing and tools
This package is used to interpolate and make adjustments to input data so that it is suitable as input to ROMS. 
##Quick Start:
1. ```pip install romspy``` to install romspy into your environment. Note that this requires you to have a C compiler which accepts CMake. Alternative approaches can be found under **Limitations**.

2. Create a script file and save this in whichever directory you desire. Guidelines found under **Options file**, examples can be found [here](https://www.github.com/saixos/romspy) under the tests folder.

3. Run the script file either from your IDE or by typing ```python options.py``` into your terminal.

4. Head out for a coffee - This could take a few minutes!

5. Navigate to your destination directory to find one or more output files.

## Script file:
 The script file determines what inputs and settings should be used. The most basic example would be:
 ```python
from romspy import PreProcessor
from romspy.settings import Settings

target_grid = "target_grid.nc"
outfile = "outfile.nc"
sources = [
    {   
        'variables': [
            {'out': 'my_variable', 'in': 'my_variable'},
        ],
        'files': ['my_input_file.nc'],
        'interpolation_method': 'bil',
    }
]

my_preprocessor = PreProcessor(target_grid, outfile, sources, Settings([]))
my_preprocessor.make()
```
To explain this line by line:

```PreProcessor``` performs all interpolation, renaming, etc. automatically based on the inputs provided.

```Settings``` runs adjustments made to variables. To give an example of an adjustment, say that the desired variable is variable a which is equal to the sum of variables a1 and a2. To produce variable a, variables a1 and a2 are first interpolated onto the desired grid, then a function in Settings is called which produces a by summing a1 and a2.

```target_grid``` is a ROMS standard grid file which is the grid that will be interpolated onto. Ensure the filepath is included if necessary.

```outfile``` is the name given to the output files. If ```outfile``` is ```"my_dir/outfile.nc"``` as in the example, then the output files will be put in the directory 'my_dir' and be named 'outfile_#_#.nc' with the #'s replaced by integers.

```sources``` contains information on which variables to interpolate and how they should be interpolated. ```sources``` is a list of dictionaries, where each dictionary is called a *group*. Each group must have three keys: ```'variables'```, ```'files'```, ```'interpolation_method'```. 
 * ```'variables'``` is a list of dictionaries, where each dictionary corresponds to a single variable. The dictionary of a single variable must have the keys ```'in'``` which is the name of the variable in the input files, and ```'out'``` which is the desired name of the variable in the output files. If a variable is to be edited by ```Settings```, then the name in ```'out'``` should correspond to the input accepted by ```Settings```. The dictionary of a single variable can also include the key ```'vertical'``` with the value ```True``` if the variable has a depth dimension and should be vertically interpolated onto an s_rho grid.
 * ```'files''``` is a list of input files where the variables can be found. Only the variables in ```'variables'``` will be extracted from the input files. The list can have multiple files if the files are separated by timesteps. If a variable is in ```'variables''``` which is not present in one of the files in ```'files'```, an error will be thrown.
 * ```'interpolation_method'``` is the horizontal interpolation method used. Vertical interpolation is always done linearly. Options are: ```'bil'```,```'bic'```,```'nn'```,```'dis'```,```'con'```,```'con2'```,```'laf'```.
 
```PreProcessor``` accepts a number of arguments detailed in **Class Overview**, but the 4 mandatory inputs are the grid to use, the output, the sources, and which Settings to use. The settings input in this example are empty, so no changes would be made. Some default settings exist and can be used by importing them ```from romspy import forcing_settings``` or ```from romspy import clim_settings```. 

```my_preprocessor``` is the ```PreProcessor``` instance with the settings as input. To run the interpolation and adjustments, call the function ```make()```.

## Settings:

Settings are created by calling ```Settings(adjustments_list, flags_dict)```.

####Adjustments:
An adjustment is a single change which should be made to the file. ```adjustments_list``` is a list of all the changes which will be made. An element of ```adjustments_list``` is a dictionary with the following keys:
 * ```'out_var_names' : set(var_name_1, var_name_2)``` Out var names points to a set of variables which this adjustment produces as output. This is used to check if the adjustment is necessary. For example, if the adjustment produces the variable ```'dQdSST'``` but you've already listed ```'dQdSST'``` as one of your ```'out'``` variables in ```sources```, then the program will not attempt to 
        'out_var_names': set(), 'in_var_names': {'sustr', 'svstr'},
        'func': str_adjustment
    },
## Customisation:

## Limitations:

###### Author information:


This is a simple example package. You can use
[Github-flavored Markdown](https://guides.github.com/features/mastering-markdown/)
to write your content.



