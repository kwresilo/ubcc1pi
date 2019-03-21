import numpy as np
import uproot
import matplotlib.pyplot as plt
from collections import Counter

def GetMag(x, y, z):
    return (x**2 + y**2 + z**2) ** 0.5

def GetAngleToZ(x, y, z):
    return np.arctan2((x**2 + y**2) ** 0.5, z)

def GetOpeningAngle(x1, y1, z1, x2, y2, z2):
    return np.arccos(((x1 * x2) + (y1 * y2) + (z1 * z2)) / (GetMag(x1, y1, z1) * GetMag(x2, y2, z2)))

# Open the file
file_name = 'truthStudy.root'
file = uproot.open(file_name)

# Get the trees
signal_tree = file['truthStudy/signalInteractions']
interactions_tree = file['truthStudy/interactions']

# Proton multiplicity
data = signal_tree.array("nProton")
nBins = 11;
bins = np.arange(nBins + 1) - 0.5
plt.hist(data, histtype='bar', bins=bins, edgecolor='C0', facecolor='None')
plt.xticks(range(nBins))
plt.xlabel('Number of protons')
plt.ylabel('Number of signal events')
plt.savefig("truthStudy_nProtons.pdf")
plt.close()

# Neutrino energy
data = signal_tree.array("nuE")
plt.hist(data, histtype='step')
plt.xlabel('Neutrino energy / GeV')
plt.ylabel('Number of signal events')
plt.savefig("truthStudy_nuEnergy.pdf")
plt.close()

# Muon momentum
muMomX = signal_tree.array("muMomX")
muMomY = signal_tree.array("muMomY")
muMomZ = signal_tree.array("muMomZ")

data = GetMag(muMomX, muMomY, muMomZ)
plt.hist(data, histtype='step')
plt.xlabel(r'Muon momentum / GeV $c^{-1}$')
plt.ylabel('Number of signal events')
plt.savefig("truthStudy_muMom.pdf")
plt.close()

# Muon angle to Z
data = GetAngleToZ(muMomX, muMomY, muMomZ)
plt.hist(data, histtype='step')
plt.xlabel('Muon angle to Z-axis / rad')
plt.ylabel('Number of signal events')
plt.savefig("truthStudy_muAngle.pdf")
plt.close()

# Pion momentum
piMomX = signal_tree.array("piMomX")
piMomY = signal_tree.array("piMomY")
piMomZ = signal_tree.array("piMomZ")

data = GetMag(piMomX, piMomY, piMomZ)
plt.hist(data, histtype='step')
plt.xlabel(r'Pion momentum / GeV $c^{-1}$')
plt.ylabel('Number of signal events')
plt.savefig("truthStudy_piMom.pdf")
plt.close()

# Pion angle to Z
data = GetAngleToZ(piMomX, piMomY, piMomZ)
plt.hist(data, histtype='step')
plt.xlabel('Pion angle to Z-axis / rad')
plt.ylabel('Number of signal events')
plt.savefig("truthStudy_piAngle.pdf")
plt.close()

# Muon/Pion opening angle
data = GetOpeningAngle(muMomX, muMomY, muMomZ, piMomX, piMomY, piMomZ)
plt.hist(data, histtype='step')
plt.xlabel('Muon-Pion opening angle / rad')
plt.ylabel('Number of signal events')
plt.savefig("truthStudy_muPiOpeningAngle.pdf")
plt.close()

# Interaction types table
print("Interaction : Count")
interactions = signal_tree.array("interaction")
frequencies = Counter(interactions)
for interaction in frequencies:
    print("{} : {}".format(interaction, frequencies[interaction]))
