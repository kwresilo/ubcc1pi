from subprocess import Popen, PIPE
import math

# Count the files in a given defintion
def count_definition_files(definition):
    p = Popen(['samweb', 'count-definition-files', definition], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output, err = p.communicate()
    rc = p.returncode

    if p.returncode != 0:
        raise Exception(err)

    return int(output)

# Get the definition name from a given data stream (bnb or extbnb) and epoch (e.g. C1)
def get_definition_name(stream, epoch):
    return "data_" + stream + "_mcc9.1_v08_00_00_25_reco2_" + epoch + "_nucc_reco2_high_lifetime"

# Print the config block for a given definition
def print_config_block(definition, dirname, maxfilesperjob = 1, useDataDriver = False):

    print ''
    print '<!-- ' + dirname + ' -->'

    numfiles = count_definition_files(definition)
    numjobs = int(math.ceil(float(numfiles) / float(maxfilesperjob)))

    print '<stage name="' + dirname + '">'
    print '    <!-- Settings -->'

    if useDataDriver:
        print '    <fcl>analysis_file_writer_driver_data.fcl</fcl>'
    else:
        print '    <fcl>analysis_file_writer_driver.fcl</fcl>'

    print '    <numjobs>' + str(numjobs) + '</numjobs>'
    print '    <maxfilesperjob>' + str(maxfilesperjob) + '</maxfilesperjob>'
    print '    <numevents>&numEvents;</numevents>'
    print ''
    print '    <!-- Input files -->'
    print '    <inputdef>' + definition + '</inputdef>'
    print ''
    print '    <!-- Output directories -->'
    print '    <outdir>/pnfs/uboone/scratch/users/&username;/&projectName;/' + dirname + '/</outdir>'
    print '    <logdir>/pnfs/uboone/scratch/users/&username;/&projectName;/' + dirname + '/</logdir>'
    print '    <workdir>/pnfs/uboone/resilient/users/&username;/work/&projectName;/' + dirname + '/</workdir>'
    print ''
    print '    <!-- Options -->'
    print '    <datatier>&dataTier;</datatier>'
    print '    <schema>&schema;</schema>'
    print '    <jobsub>&jobsub;</jobsub>'
    print '</stage>'

# Print the config block for a given stream and epoch and run number
def print_config_block_data(stream, run, epoch):

    definition = get_definition_name(stream, epoch)
    dirname = stream + '_run' + run + '-' + epoch

    print_config_block(definition, dirname, 10, True)

# ------------------------------------------------------------------------------------------------------------------------------------------

# Print the boilerplate
print '<?xml version="1.0"?>'
print ''
print '<!DOCTYPE project ['
print ''
print '<!-- Names -->'
print '<!ENTITY username        "asmith">'
print ''
print '<!-- LArSoft version -->'
print '<!ENTITY tag             "v08_00_00_25">'
print '<!ENTITY qual            "e17:prof">'
print ''
print '<!-- Job settings -->'
print '<!ENTITY numEvents       "9999999">'
print ''
print '<!-- Shorthands -->'
print '<!ENTITY projectName     "ubcc1pi">'
print '<!ENTITY dataTier        "ana">'
print '<!ENTITY schema          "root">'
print '<!ENTITY jobsub           "--append_condor_requirements=' + "'(TARGET.HAS_CVMFS_uboone_osgstorage_org==true)'" + '" -e XRD_LOADBALANCERTTL=7200 -e XRD_CONNECTIONRETRY=32 -e XRD_REQUESTTIMEOUT=3600 -e XRD_REDIRECTLIMIT=255>'
print ''
print ']>'
print ''
print '<project name="&projectName;">'
print ''
print '    <!-- =============================================================================================================================== -->'
print '    <!-- Set up the LArSoft version                                                                                                      -->'
print '    <!-- =============================================================================================================================== -->'
print '    <larsoft>'
print '        <tag>&tag;</tag>'
print '        <qual>&qual;</qual>'
print '        <local>/pnfs/uboone/resilient/users/&username;/tars/&projectName;/local.tar</local>'
print '    </larsoft>'
print '    '
print '    <!-- =============================================================================================================================== -->'
print '    <!-- Set up the project -->'
print '    <!-- =============================================================================================================================== -->'
print '    <group>uboone</group>'
print '    <os>SL7</os>'
print '    <resource>DEDICATED,OPPORTUNISTIC</resource>'
print ''
print '    <!-- To be overridden -->'
print '    <numjobs>0</numjobs>'
print '    <maxfilesperjob>0</maxfilesperjob>'
print '    <numevents>&numEvents;</numevents>'

# Monte-carlo
print_config_block("prodgenie_bnb_nu_uboone_overlay_ccinc_reweight_mcc9.1_run1_joelam", "overlays")
print_config_block("prodgenie_bnb_dirt_overlay_ccinc_withWeights_joelam", "dirt")

# Beam on data
print_config_block_data("bnb", "1", "C1")
print_config_block_data("bnb", "2", "D2")
print_config_block_data("bnb", "2", "E1")
print_config_block_data("bnb", "3", "F")
print_config_block_data("bnb", "3", "G1")

# Beam off data
print_config_block_data("extbnb", "1", "C1")
print_config_block_data("extbnb", "1", "C2")
print_config_block_data("extbnb", "2", "D2")
print_config_block_data("extbnb", "2", "E1")
print_config_block_data("extbnb", "2", "E2")
print_config_block_data("extbnb", "3", "F")
print_config_block_data("extbnb", "3", "G1")
print_config_block_data("extbnb", "3", "G2")
print_config_block_data("extbnb", "3", "G2a")

print '</project>'
