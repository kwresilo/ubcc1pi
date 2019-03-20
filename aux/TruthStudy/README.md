# Truth Study

## Overview
In this MC truth study, we count the fraction of neutrino events that have the desired CC1Pi topology, and look at the kinematics of the
interactions. In order to determine if an event contiains a CC1Pi topology, we do the following...

- Find the MCTruth object that corresponds to the most energetic beam neutrino in the event (NB. no significant pileup is expected).
- Get the primary MCParticles from the neutrino interaction
- Navigate through the MCParticle hierarchy of each primary until we find a "visible" final state particle
    - Here visible = one that can leave a track or shower in the detector
    - E.g. for a pi0 primary, we consider the two decay photons to be the visible final states
- Apply thresholds on the momentum of these particles to avoid considering particles which have too little energy to leave a signal
- Insist that of these "reconstructable" final states, there is exactly:
    - 1 Mu-
    - 1 Pi+
    - N p
    - Nothing else (NB. neutrons are not deemed visible)
- Further to this, we also insist that the neutrino interaction vertex is within the fiducial volume

## Instructions

### Instructions for running locally on a single file
```
# Run the TruthStudy module
lar -c truth_study_driver.fcl -s <my_input_art_root_file>

# Now you should have an output root file
root -l truthStudy.root
```

## Output file content

The output file contains two trees.
- The first contains details of all interactions, and has branches as documented [here](https://a-d-smith.github.io/ubcc1pi/html/structubcc1pi_1_1TruthStudy_1_1InteractionOutput.html).
- The second contains details of only CC1Pi signal interactions, and has branches as documented [here](https://a-d-smith.github.io/ubcc1pi/html/structubcc1pi_1_1TruthStudy_1_1SignalOutput.html).
