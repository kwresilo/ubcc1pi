# Standalone analysis code

## Installation
In this directory the `CMakeLists.txt` file is used when compiling the package using LArSoft. To run the analysis code outside of LArSoft,
just start up root with the compile script.

```
# Check you are in the ubcc1pi_standalone directory
pwd

# Start root and compile the code
root -l compile.C
```


## Running a macro
```
# Start root and load the compiled code
root -l compile.C

# Run the macro you wish
Test()
```


## Writing a new macro
If you want to write a new macro, just add it to the compile.C script to make sure it's picked up!
