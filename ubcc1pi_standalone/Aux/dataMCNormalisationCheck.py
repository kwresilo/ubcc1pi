import math

# Beam on data
nEvents_dataBNB = 28678
POT_dataBNB = 1.467e+20 # tor875_wcut
triggers_dataBNB = 32603210.0 # E1DCNT_wcut

# Overlays
nEvents_overlays = 127407
POT_overlays = 1.089591e+21

# Beam off data
nEvents_dataEXT = 10384
triggers_dataEXT = 62920813.0 # EXT

# ------------------------------------------------------------------------------------------------------------------------------------------

# Get the errors on the values assuming Poisson stats
nEvents_dataBNB_err = math.sqrt(nEvents_dataBNB)
POT_dataBNB_err = math.sqrt(POT_dataBNB)
triggers_dataBNB_err = math.sqrt(triggers_dataBNB)
nEvents_overlays_err = math.sqrt(nEvents_overlays)
POT_overlays_err = math.sqrt(POT_overlays)
nEvents_dataEXT_err = math.sqrt(nEvents_dataEXT)
triggers_dataEXT_err = math.sqrt(triggers_dataEXT)

# Get the scaling factors
scale_overlays = POT_dataBNB / POT_overlays
scale_overlays_err = scale_overlays * math.sqrt(math.pow(POT_dataBNB_err / POT_dataBNB, 2) + math.pow(POT_overlays_err / POT_overlays, 2))

scale_dataEXT = triggers_dataBNB / triggers_dataEXT
scale_dataEXT_err = scale_dataEXT * math.sqrt(math.pow(triggers_dataBNB_err / triggers_dataBNB, 2) + math.pow(triggers_dataEXT_err / triggers_dataEXT, 2))

# Scale the number of events
scaledEvents_overlays = nEvents_overlays * scale_overlays
scaledEvents_overlays_err = scaledEvents_overlays * math.sqrt(math.pow(nEvents_overlays_err / nEvents_overlays, 2) + math.pow(scale_overlays_err / scale_overlays, 2))

scaledEvents_dataEXT = nEvents_dataEXT * scale_dataEXT
scaledEvents_dataEXT_err = scaledEvents_dataEXT * math.sqrt(math.pow(nEvents_dataEXT_err / nEvents_dataEXT, 2) + math.pow(scale_dataEXT_err / scale_dataEXT, 2))

scaledEvents_overlayEXT = scaledEvents_overlays + scaledEvents_dataEXT
scaledEvents_overlayEXT_err = math.sqrt(math.pow(scaledEvents_overlays_err, 2) + math.pow(scaledEvents_dataEXT_err, 2))

ratio = nEvents_dataBNB / scaledEvents_overlayEXT
ratio_err = ratio * math.sqrt(math.pow(nEvents_dataBNB_err / nEvents_dataBNB, 2) + math.pow(scaledEvents_overlayEXT_err / scaledEvents_overlayEXT, 2))

# Print the output information
print "Beam on data"
print "nEvents      : {} +- {}".format(nEvents_dataBNB, nEvents_dataBNB_err)
print "POT          : {} +- {}".format(POT_dataBNB, POT_dataBNB_err)
print "triggers     : {} +- {}".format(triggers_dataBNB, triggers_dataBNB_err)
print

print "Overlays"
print "nEvents      : {} +- {}".format(nEvents_overlays, nEvents_overlays_err)
print "POT          : {} +- {}".format(POT_overlays, POT_overlays_err)
print "scaleFactor  : {} +- {}".format(scale_overlays, scale_overlays_err)
print "scaledEvents : {} +- {}".format(scaledEvents_overlays, scaledEvents_overlays_err)
print

print "Beam off data"
print "nEvents      : {} +- {}".format(nEvents_dataEXT, nEvents_dataEXT_err)
print "triggers     : {} +- {}".format(triggers_dataEXT, triggers_dataEXT_err)
print "scaleFactor  : {} +- {}".format(scale_dataEXT, scale_dataEXT_err)
print "scaledEvents : {} +- {}".format(scaledEvents_dataEXT, scaledEvents_dataEXT_err)
print

print "Overlay + beam off"
print "scaledEvents : {} +- {}".format(scaledEvents_overlayEXT, scaledEvents_overlayEXT_err)
print "ratio        : {} +- {}".format(ratio, ratio_err)
