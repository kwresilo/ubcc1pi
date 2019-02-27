# Documentation

This repository is documented using doxygen. To view the documentation for the current master branch, please go [here](https://www.hep.phy.cam.ac.uk/~asmith/ubcc1pi-docs/).
If you need to compile the documentation for another branch / commit, please `git clone` this repository and then follow these steps:

```
# Let doxygen create the documentation
cd ubcc1pi/docs
doxygen

# View the documentation in your web browser
#     Please note, if you are using an ssh session then this will open up firefox on the remote machine.
#     You may prefer to copy the docs to your local machine, or if you just need the master branch - view the docs hosted above.
firefox html/index.html &
```
