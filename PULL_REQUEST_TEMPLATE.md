**Sanity checklist**
- [ ] Is the feature rebased on top the tip of the branch with which you wish to merge?
- [ ] Does `mrbsetenv; mrb t -j4` come back clean?
- [ ] Does `ubcc1pi` come back clean? (you may need to run `ubcc1pi_standalone/setup_standalone.sh` first)
- [ ] Have you sourced `cleanWhitespace.sh`? Feel free to make a commit before and after you run this one if you are worried!
- [ ] Does `cd docs; ./makedox` come back clean?

# Details of changes

Details your changes here.
