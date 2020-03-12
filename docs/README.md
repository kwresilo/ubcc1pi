# Documentation

This repository is documented using doxygen. To view the documentation for the current master branch, please go <https://a-d-smith.github.io/ubcc1pi>.
If you need to compile the documentation for another branch / commit, please `git clone` this repository, `git checkout` your desired branch
or commit and then follow these steps:

```
# Make the documentation using doxygen
cd ubcc1pi/docs
./makedox

# If you are developing the code... if there is any undocumented code, this script will tell you what and where you need to document!
# Please do that, and then re-run makedox 

# View the documentation in your web browser
#     Please note, if you are using an ssh session then this will open up firefox on the remote machine.
#     You may prefer to use a VNC or copy the docs to your local machine. If you just need the master branch - just view the docs hosted above.
firefox index.html &
```
