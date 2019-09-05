# ROMSPY: Regional Ocean Modelling Systems preprocessing and tools
This package is used to interpolate and make adjustments to input data so that it is suitable as input to ROMS. 
##Quick Start:
1. ```pip install romspy``` to install romspy into your environment. Note that this requires you to have a C compiler which accepts CMake. Alternative approaches can be found under **Limitations**.

2. Create an options file and save this in whichever directory you desire. Guidelines found under **Options file**, examples can be found [here](https://www.github.com/saixos/romspy) under the tests folder.

3. Run the options file either from your IDE or by typing ```python options.py``` into your terminal.

4. Head out for a coffee - This could take a few minutes!

5. Navigate to your destination directory to find one or more output files.

## Options file:
 The options file determines what inputs and settings should be used. The most basic example would be:
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

my_preprocessor = PreProcessor(target_grid, outfile, sources, Settings({}))
my_preprocessor.make()
```
To explain this line by line:

```PreProcessor``` is the main class used by romspy. It performs all interpolation, renaming, etc. automatically based on the inputs provided.

```Settings``` is a class which runs adjustments made to variables. To give an example of an adjustment, say that the desired variable is variable a which is equal to the sum of variables a1 and a2. To produce variable a, variables a1 and a2 are first interpolated onto the desired grid, then a function in Settings is called which produces a by summing a1 and a2.

```target_grid``` is a ROMS standard grid file which is the grid that will be interpolated onto. Ensure the filepath is included if necessary.

```outfile``` is the name given to the output files. If ```outfile``` is ```"my_dir/outfile.nc"``` as in the example, then the output files will be put in the directory 'my_dir' and be named 'outfile_#_#.nc', with the first # replaced by the group number and the second # replaced by the file number.

```sources``` contains information on which variables to interpolate and how they should be interpolated. ```sources``` is a list of dictionaries, where each dictionary is called a *group*. Each group must have three keys: ```'variables'```, ```'files'```, ```'interpolation_method'```. 
 * ```'variables'``` is a list of dictionaries, where each dictionary corresponds to a single variable. The dictionary of a single variable must have the keys ```'in'``` which is the name of the variable in the input files, and ```'out'``` which is the desired name of the variable in the output files. If a variable is to be edited by ```Settings```, then the name in ```'out'``` should correspond to the input accepted by ```Settings```. The dictionary of a single variable can also include the key ```'vertical'``` with the value ```True``` if the variable has a depth dimension and should be vertically interpolated onto an s_rho grid.
 * ```'files''``` is a list of input files where the variables can be found. Only the variables in ```'variables'``` will be extracted from the input files. The list can have multiple files if the files are separated by timesteps. If a variable is in ```'variables''``` which is not present in one of the files in ```'files'```, an error will be thrown.
 * ```'interpolation_method'``` is the horizontal interpolation method used. Vertical interpolation is always done linearly. Options are: ```''```,


## Class Overview:

## Customisation:

## Limitations:

###### Author information:


This is a simple example package. You can use
[Github-flavored Markdown](https://guides.github.com/features/mastering-markdown/)
to write your content.



