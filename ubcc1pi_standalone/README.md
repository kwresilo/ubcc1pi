# Standalone analysis code

## Installation
In this directory the [CMakeLists.txt](ubcc1pi_standalone/CMakeLists.txt) file is used when compiling the package using LArSoft. To run the
analysis code outside of LArSoft, just start up root and compile the required code using the instruction below.

Run the set-up script for the standalone environment
```
source setup_standalone.sh
```

This script sets the environment variable `UBCC1PI_STANDALONE_DIR` and makes an alias (`ubcc1pi = "root -l -b compile.C"`) that starts up
root and compiles the required standalone code. If the code is already compiled then `ubcc1pi` will just start root and load it.
```
ubcc1pi
```


## Running a macro
As above, start root and load the compiled code using the alias:
```
ubcc1pi
```

From the root command-line, run the macro you want with the default configuration
```
ubcc1pi_macros::RunFullAnalysis();
```

To see the full list of available macros, you can type `ubcc1pi_macros::` and press TAB.


## Modifying the configuration
Start root and load the compiled code:
```
ubcc1pi
```

The configuration options are handled by the `ubcc1pi::Config` structure. Make an instance of it using the root command line:
```
ubcc1pi::Config config;
```

Print the configuration to see the defaults
```
ubcc1pi_macros::PrintConfig(config);
```

You can modify a parameter as follows - for a full list of parameters, see
[ubcc1pi_standalone/Objects/Config.h](ubcc1pi_standalone/Objects/Config.h). You can also press TAB to see the options.
```
config.files.overlaysFileName = "/path/to/my/file.root";
```

To run a macro with your custom configuration just pass it as the first argument
```
ubcc1pi_macros::RunFullAnalysis(config);
```

If you want to change the default configuration, then just modify [ubcc1pi_standalone/Objects/Config.h](ubcc1pi_standalone/Objects/Config.h).


## Writing a new macro
First declare the macro in [ubcc1pi_standalone/Macros/Macros.h](ubcc1pi_standalone/Macros/Macros.h). For example:
```
    /**
     *  @brief  Description of my new macro
     *
     *  @param  config the input config
     */
    void MyNewMacro(const Config &config = Config());
```

By convention, the macro should only take one argument which is the `ubcc1pi::Config` object. Use `const Config &config = Config()` so the
default configuration can be easily used. You can add any configuration options required into
[ubcc1pi_standalone/Objects/Config.h](ubcc1pi_standalone/Objects/Config.h).

Now add the actual macro in the [ubcc1pi_standalone/Macros/](ubcc1pi_standalone/Macros/) directory, take a look at the existing macros as an
example - you should use the `ubcc1pi_macros` namespace.

Finally, add the macro to the compilation script [ubcc1pi_standalone/compile.C](ubcc1pi_standalone/compile.C) - you just need to add the
path to your macro with the rest.
